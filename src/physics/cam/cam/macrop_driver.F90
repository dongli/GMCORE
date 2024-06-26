module macrop_driver

  !-----------------------------------------------------------------------------
  ! Purpose:
  !
  ! Provides the CAM interface to the prognostic cloud macrophysics
  !
  ! Author: Andrew Gettelman, Cheryl Craig October 2010
  ! Origin: modified from stratiform.F90 elements
  !    (Boville 2002, Coleman 2004, Park 2009, Kay 2010)
  !-----------------------------------------------------------------------------

  use shr_kind_mod  , only: r8=>shr_kind_r8
  use spmd_utils    , only: masterproc
  use ppgrid        , only: pcols, pver, pverp
  use physconst     , only: latice, latvap
  use phys_control  , only: phys_getopts
  use constituents  , only: cnst_get_ind, pcnst
  use physics_buffer, only: physics_buffer_desc, pbuf_set_field, pbuf_get_field, pbuf_old_tim_idx
  use time_manager  , only: is_first_step
  use cldwat2m_macro, only: ini_macro
  use perf_mod      , only: t_startf, t_stopf
  use cam_logfile   , only: iulog
  use cam_abortutils, only: endrun
  use zm_conv_intr  , only: zmconv_microp

  implicit none

  private
  
  save

  public macrop_driver_readnl
  public macrop_driver_register
  public macrop_driver_init
  public macrop_driver_tend
  public liquid_macro_tend

  logical, public :: do_cldice  ! .true., park macrophysics is prognosing cldice
  logical, public :: do_cldliq  ! .true., park macrophysics is prognosing cldliq
  logical, public :: do_detrain ! .true., park macrophysics is detraining ice into stratiform

  ! 'cu_det_st' : If .true. (.false.), detrain cumulus liquid condensate into the pre-existing liquid stratus
  !               (environment) without (with) macrophysical evaporation. If there is no pre-esisting stratus,
  !               evaporate cumulus liquid condensate. This option only influences the treatment of cumulus
  !               liquid condensate, not cumulus ice condensate.
  logical, parameter :: cu_det_st = .false.

  ! Parameters used for selecting generalized critical RH for liquid and ice stratus
  integer :: rhminl_opt = 0
  integer :: rhmini_opt = 0


  character(16) :: shallow_scheme
  logical       :: use_shfrc      ! Local copy of flag from convect_shallow_use_shfrc

  
  integer ixcldliq    ! Cloud liquid amount index
  integer ixcldice    ! Cloud ice amount index
  integer ixnumliq    ! Cloud liquid number index
  integer ixnumice    ! Cloud ice water index
  integer qcwat_idx   !
  integer lcwat_idx   !
  integer iccwat_idx  !
  integer nlwat_idx   !
  integer niwat_idx   !
  integer tcwat_idx   !
  integer cc_t_idx    !
  integer cc_qv_idx   !
  integer cc_ql_idx   !
  integer cc_qi_idx   !
  integer cc_nl_idx   !
  integer cc_ni_idx   !
  integer cc_qlst_idx !
  integer cld_idx     ! Cloud fraction index
  integer ast_idx     ! Stratiform cloud fraction index
  integer aist_idx    ! Ice stratiform cloud fraction index
  integer alst_idx    ! Liquid stratiform cloud fraction index in physics buffer
  integer qist_idx    ! Ice stratiform in-cloud IWC
  integer qlst_idx    ! Liquid stratiform in-cloud LWC
  integer concld_idx  ! Convective cloud index
  integer fice_idx
  integer cmeliq_idx
  integer shfrc_idx

  integer :: dlfzm_idx    = -1 ! ZM detrained convective cloud water mixing ratio.
  integer :: difzm_idx    = -1 ! ZM detrained convective cloud ice mixing ratio.
  integer :: dnlfzm_idx   = -1 ! ZM detrained convective cloud water num concen.
  integer :: dnifzm_idx   = -1 ! ZM detrained convective cloud ice num concen.

  integer :: tke_idx      = -1 ! TKE defined at the model interfaces
  integer :: qtl_flx_idx  = -1 ! Overbar(w'qtl' where qtl = qv + ql) from the PBL scheme
  integer :: qti_flx_idx  = -1 ! Overbar(w'qti' where qti = qv + qi) from the PBL scheme
  integer :: cmfr_det_idx = -1 ! Detrained convective mass flux from UNICON
  integer :: qlr_det_idx  = -1 ! Detrained convective ql from UNICON
  integer :: qir_det_idx  = -1 ! Detrained convective qi from UNICON
  integer :: cmfmc_sh_idx = -1

contains

  subroutine macrop_driver_readnl(nlfile)

    use namelist_utils, only: find_group_name
    use units         , only: getunit, freeunit
    use mpishorthand

    character(*), intent(in) :: nlfile

    logical :: macro_park_do_cldice  = .true.   ! do_cldice = .true., park macrophysics is prognosing cldice
    logical :: macro_park_do_cldliq  = .true.   ! do_cldliq = .true., park macrophysics is prognosing cldliq
    logical :: macro_park_do_detrain = .true.   ! do_detrain = .true., park macrophysics is detraining ice into stratiform

    integer unitn, ierr
    character(*), parameter :: subname = 'macrop_driver_readnl'

    namelist /macro_park_nl/ macro_park_do_cldice, macro_park_do_cldliq, macro_park_do_detrain

    if (masterproc) then
      unitn = getunit()
      open(unitn, file=trim(nlfile), status='old')
      call find_group_name(unitn, 'macro_park_nl', status=ierr)
      if (ierr == 0) then
        read(unitn, macro_park_nl, iostat=ierr)
        if (ierr /= 0) then
          call endrun(subname // ':: ERROR reading namelist')
        end if
      end if
      close(unitn)
      call freeunit(unitn)

      do_cldice  = macro_park_do_cldice
      do_cldliq  = macro_park_do_cldliq
      do_detrain = macro_park_do_detrain
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(do_cldice , 1, mpilog, 0, mpicom)
    call mpibcast(do_cldliq , 1, mpilog, 0, mpicom)
    call mpibcast(do_detrain, 1, mpilog, 0, mpicom)
#endif

  end subroutine macrop_driver_readnl

  subroutine macrop_driver_register()

    !---------------------------------------------------------------------- !
    !                                                                       !
    ! Register the constituents (cloud liquid and cloud ice) and the fields !
    ! in the physics buffer.                                                !
    !                                                                       !
    !---------------------------------------------------------------------- !

    use physics_buffer, only: pbuf_add_field, dtype_r8, dyn_time_lvls

    call phys_getopts(shallow_scheme_out=shallow_scheme)

    call pbuf_add_field('AST'   , 'global' , dtype_r8, [pcols,pver,dyn_time_lvls], ast_idx   )
    call pbuf_add_field('AIST'  , 'global' , dtype_r8, [pcols,pver,dyn_time_lvls], aist_idx  )
    call pbuf_add_field('ALST'  , 'global' , dtype_r8, [pcols,pver,dyn_time_lvls], alst_idx  )
    call pbuf_add_field('QIST'  , 'global' , dtype_r8, [pcols,pver,dyn_time_lvls], qist_idx  )
    call pbuf_add_field('QLST'  , 'global' , dtype_r8, [pcols,pver,dyn_time_lvls], qlst_idx  )
    call pbuf_add_field('CLD'   , 'global' , dtype_r8, [pcols,pver,dyn_time_lvls], cld_idx   )
    call pbuf_add_field('CONCLD', 'global' , dtype_r8, [pcols,pver,dyn_time_lvls], concld_idx)
    call pbuf_add_field('QCWAT' , 'global' , dtype_r8, [pcols,pver,dyn_time_lvls], qcwat_idx )
    call pbuf_add_field('LCWAT' , 'global' , dtype_r8, [pcols,pver,dyn_time_lvls], lcwat_idx )
    call pbuf_add_field('ICCWAT', 'global' , dtype_r8, [pcols,pver,dyn_time_lvls], iccwat_idx)
    call pbuf_add_field('NLWAT' , 'global' , dtype_r8, [pcols,pver,dyn_time_lvls], nlwat_idx )
    call pbuf_add_field('NIWAT' , 'global' , dtype_r8, [pcols,pver,dyn_time_lvls], niwat_idx )
    call pbuf_add_field('TCWAT' , 'global' , dtype_r8, [pcols,pver,dyn_time_lvls], tcwat_idx )
    call pbuf_add_field('FICE'  , 'physpkg', dtype_r8, [pcols,pver              ], fice_idx  )
    call pbuf_add_field('CMELIQ', 'physpkg', dtype_r8, [pcols,pver              ], cmeliq_idx)

  end subroutine macrop_driver_register

  subroutine macrop_driver_init(pbuf2d)

    !-------------------------------------------- !
    !                                             !
    ! Initialize the cloud water parameterization !
    !                                             !
    !-------------------------------------------- !

    use physics_buffer , only: pbuf_get_index
    use cam_history    , only: addfld, add_default
    use convect_shallow, only: convect_shallow_use_shfrc

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    logical history_aerosol ! Output the MAM aerosol tendencies
    logical history_budget  ! Output tendencies and state variables for CAM4
                            ! temperature, water vapor, cloud ice and cloud
                            ! liquid budgets.
    integer history_budget_histfile_num ! Output history file number for budget fields
    integer istat
    character(*), parameter :: subname = 'macrop_driver_init'

    ! Initialization routine for cloud macrophysics
    if (shallow_scheme == 'UNICON') rhminl_opt = 1
    call ini_macro(rhminl_opt, rhmini_opt)

    call phys_getopts(history_aerosol_out            =history_aerosol      , &
                      history_budget_out             =history_budget       , &
                      history_budget_histfile_num_out=history_budget_histfile_num)

    ! Find out whether shfrc from convect_shallow will be used in cldfrc

    if (convect_shallow_use_shfrc()) then
      use_shfrc = .true.
      shfrc_idx = pbuf_get_index('shfrc')
    else
      use_shfrc = .false.
    endif

    call addfld('DPDLFLIQ' , ['lev'], 'A', 'kg/kg/s' , 'Detrained liquid water from deep convection'       )
    call addfld('DPDLFICE' , ['lev'], 'A', 'kg/kg/s' , 'Detrained ice from deep convection'                )
    call addfld('SHDLFLIQ' , ['lev'], 'A', 'kg/kg/s' , 'Detrained liquid water from shallow convection'    )
    call addfld('SHDLFICE' , ['lev'], 'A', 'kg/kg/s' , 'Detrained ice from shallow convection'             )
    call addfld('DPDLFT'   , ['lev'], 'A', 'K/s'     , 'T-tendency due to deep convective detrainment'     )
    call addfld('SHDLFT'   , ['lev'], 'A', 'K/s'     , 'T-tendency due to shallow convective detrainment'  )
    call addfld('ZMDLF'    , ['lev'], 'A', 'kg/kg/s' , 'Detrained liquid water from ZM convection'         )
    call addfld('MACPDT'   , ['lev'], 'A', 'W/kg'    , 'Heating tendency - Revised  macrophysics'          )
    call addfld('MACPDQ'   , ['lev'], 'A', 'kg/kg/s' , 'Q tendency - Revised macrophysics'                 )
    call addfld('MACPDLIQ' , ['lev'], 'A', 'kg/kg/s' , 'CLDLIQ tendency - Revised macrophysics'            )
    call addfld('MACPDICE' , ['lev'], 'A', 'kg/kg/s' , 'CLDICE tendency - Revised macrophysics'            )
    call addfld('CLDVAPADJ', ['lev'], 'A', 'kg/kg/s' , 'Q tendency associated with liq/ice adjustment - Revised macrophysics')
    call addfld('CLDLIQADJ', ['lev'], 'A', 'kg/kg/s' , 'CLDLIQ adjustment tendency - Revised macrophysics' )
    call addfld('CLDICEADJ', ['lev'], 'A', 'kg/kg/s' , 'CLDICE adjustment tendency - Revised macrophysics' )
    call addfld('CLDLIQDET', ['lev'], 'A', 'kg/kg/s' , 'Detrainment of conv cld liq into envrionment  - Revised macrophysics')
    call addfld('CLDICEDET', ['lev'], 'A', 'kg/kg/s' , 'Detrainment of conv cld ice into envrionment  - Revised macrophysics')
    call addfld('CLDLIQLIM', ['lev'], 'A', 'kg/kg/s' , 'CLDLIQ limiting tendency - Revised macrophysics'   )
    call addfld('CLDICELIM', ['lev'], 'A', 'kg/kg/s' , 'CLDICE limiting tendency - Revised macrophysics'   )
    call addfld('AST'      , ['lev'], 'A', '1'       , 'Stratus cloud fraction'                            )
    call addfld('LIQCLDF'  , ['lev'], 'A', '1'       , 'Stratus Liquid cloud fraction'                     )
    call addfld('ICECLDF'  , ['lev'], 'A', '1'       , 'Stratus ICE cloud fraction'                        )
    call addfld('CLDST'    , ['lev'], 'A', 'fraction', 'Stratus cloud fraction'                            )
    call addfld('CONCLD'   , ['lev'], 'A', 'fraction', 'Convective cloud cover'                            )
    call addfld('CLR_LIQ'  , ['lev'], 'A', 'fraction', 'Clear sky fraction for liquid stratus'             )
    call addfld('CLR_ICE'  , ['lev'], 'A', 'fraction', 'Clear sky fraction for ice stratus'                )
    call addfld('CLDLIQSTR', ['lev'], 'A', 'kg/kg'   , 'Stratiform CLDLIQ'                                 )
    call addfld('CLDICESTR', ['lev'], 'A', 'kg/kg'   , 'Stratiform CLDICE'                                 )
    call addfld('CLDLIQCON', ['lev'], 'A', 'kg/kg'   , 'Convective CLDLIQ'                                 )
    call addfld('CLDICECON', ['lev'], 'A', 'kg/kg'   , 'Convective CLDICE'                                 )
    call addfld('CLDSICE'  , ['lev'], 'A', 'kg/kg'   , 'CloudSat equivalent ice mass mixing ratio'         )
    call addfld('CMELIQ'   , ['lev'], 'A', 'kg/kg/s' , 'Rate of cond-evap of liq within the cloud'         )
    call addfld('TTENDICE' , ['lev'], 'A', 'K/s'     , 'T tendency from Ice Saturation Adjustment'         )
    call addfld('QVTENDICE', ['lev'], 'A', 'kg/kg/s' , 'Q tendency from Ice Saturation Adjustment'         )
    call addfld('QITENDICE', ['lev'], 'A', 'kg/kg/s' , 'CLDICE tendency from Ice Saturation Adjustment'    )
    call addfld('NITENDICE', ['lev'], 'A', 'kg/kg/s' , 'NUMICE tendency from Ice Saturation Adjustment'    )
    if (history_budget) then
      call add_default ('DPDLFLIQ ', history_budget_histfile_num, ' ')
      call add_default ('DPDLFICE ', history_budget_histfile_num, ' ')
      call add_default ('SHDLFLIQ ', history_budget_histfile_num, ' ')
      call add_default ('SHDLFICE ', history_budget_histfile_num, ' ')
      call add_default ('DPDLFT   ', history_budget_histfile_num, ' ')
      call add_default ('SHDLFT   ', history_budget_histfile_num, ' ')
      call add_default ('ZMDLF    ', history_budget_histfile_num, ' ')
      call add_default ('MACPDT   ', history_budget_histfile_num, ' ')
      call add_default ('MACPDQ   ', history_budget_histfile_num, ' ')
      call add_default ('MACPDLIQ ', history_budget_histfile_num, ' ')
      call add_default ('MACPDICE ', history_budget_histfile_num, ' ')
      call add_default ('CLDVAPADJ', history_budget_histfile_num, ' ')
      call add_default ('CLDLIQLIM', history_budget_histfile_num, ' ')
      call add_default ('CLDLIQDET', history_budget_histfile_num, ' ')
      call add_default ('CLDLIQADJ', history_budget_histfile_num, ' ')
      call add_default ('CLDICELIM', history_budget_histfile_num, ' ')
      call add_default ('CLDICEDET', history_budget_histfile_num, ' ')
      call add_default ('CLDICEADJ', history_budget_histfile_num, ' ')
      call add_default ('CMELIQ   ', history_budget_histfile_num, ' ')
    end if

    ! Get constituent indices
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('NUMLIQ', ixnumliq)
    call cnst_get_ind('NUMICE', ixnumice)

    ! Get physics buffer indices
    cc_t_idx     = pbuf_get_index('CC_T'    )
    cc_qv_idx    = pbuf_get_index('CC_qv'   )
    cc_ql_idx    = pbuf_get_index('CC_ql'   )
    cc_qi_idx    = pbuf_get_index('CC_qi'   )
    cc_nl_idx    = pbuf_get_index('CC_nl'   )
    cc_ni_idx    = pbuf_get_index('CC_ni'   )
    cc_qlst_idx  = pbuf_get_index('CC_qlst' )
    cmfmc_sh_idx = pbuf_get_index('CMFMC_SH')

    if (zmconv_microp) then
      dlfzm_idx  = pbuf_get_index('DLFZM' )
      difzm_idx  = pbuf_get_index('DIFZM' )
      dnlfzm_idx = pbuf_get_index('DNLFZM')
      dnifzm_idx = pbuf_get_index('DNIFZM')
    end if

    if (rhminl_opt > 0 .or. rhmini_opt > 0) then
      cmfr_det_idx = pbuf_get_index('cmfr_det', istat)
      if (istat < 0) call endrun(subname // ': macrop option requires cmfr_det in pbuf')
      if (rhminl_opt > 0) then
        qlr_det_idx  = pbuf_get_index('qlr_det', istat)
        if (istat < 0) call endrun(subname // ': macrop option requires qlr_det in pbuf')
      end if
      if (rhmini_opt > 0) then
        qir_det_idx  = pbuf_get_index('qir_det', istat)
        if (istat < 0) call endrun(subname // ': macrop option requires qir_det in pbuf')
      end if
    end if

    if (rhminl_opt == 2 .or. rhmini_opt == 2) then
      tke_idx = pbuf_get_index('tke')
      if (rhminl_opt == 2) then
        qtl_flx_idx = pbuf_get_index('qtl_flx', istat)
        if (istat < 0) call endrun(subname // ': macrop option requires qtl_flx in pbuf')
      end if
      if (rhmini_opt == 2) then
        qti_flx_idx = pbuf_get_index('qti_flx', istat)
        if (istat < 0) call endrun(subname // ': macrop option requires qti_flx in pbuf')
      end if
    end if

    ! Init pbuf fields.  Note that the fields CLD, CONCLD, QCWAT, LCWAT,
    ! ICCWAT, and TCWAT are initialized in phys_inidat.
    if (is_first_step()) then
      call pbuf_set_field(pbuf2d, ast_idx  ,  0.0_r8)
      call pbuf_set_field(pbuf2d, aist_idx ,  0.0_r8)
      call pbuf_set_field(pbuf2d, alst_idx ,  0.0_r8)
      call pbuf_set_field(pbuf2d, qist_idx ,  0.0_r8)
      call pbuf_set_field(pbuf2d, qlst_idx ,  0.0_r8)
      call pbuf_set_field(pbuf2d, nlwat_idx,  0.0_r8)
      call pbuf_set_field(pbuf2d, niwat_idx,  0.0_r8)
    end if

    ! The following are physpkg, so they need to be init every time
    call pbuf_set_field(pbuf2d, fice_idx  , 0.0_r8)
    call pbuf_set_field(pbuf2d, cmeliq_idx, 0.0_r8)

  end subroutine macrop_driver_init

  subroutine macrop_driver_tend(state, ptend, dtime, landfrac, &
                                ocnfrac, snowh, dlf, dlf2, cmfmc, ts, &
                                sst, zdu, pbuf, det_s, det_ice)

    !-------------------------------------------------------- !
    !                                                         !
    ! Purpose:                                                !
    !                                                         !
    ! Interface to detrain, cloud fraction and                !
    !     cloud macrophysics subroutines                      !
    !                                                         !
    ! Author: A. Gettelman, C. Craig, Oct 2010                !
    ! based on stratiform_tend by D.B. Coleman 4/2010         !
    !                                                         !
    !-------------------------------------------------------- !

    use cloud_fraction, only: cldfrc, cldfrc_fice
    use physics_types , only: physics_state, physics_ptend
    use physics_types , only: physics_ptend_init, physics_update
    use physics_types , only: physics_ptend_sum,  physics_state_copy
    use physics_types , only: physics_state_dealloc
    use cam_history   , only: outfld
    use constituents  , only: cnst_get_ind, pcnst
    use cldwat2m_macro, only: mmacro_pcond
    use physconst     , only: cpair, tmelt, gravit
    use time_manager  , only: get_nstep
    use ref_pres      , only: top_lev => trop_cloud_top_lev

    type(physics_state), intent(in)    :: state       ! State variables
    type(physics_ptend), intent(out)   :: ptend       ! macrophysics parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)     ! Physics buffer

    real(r8), intent(in) :: dtime                     ! Timestep
    real(r8), intent(in) :: landfrac(pcols)           ! Land fraction (fraction)
    real(r8), intent(in) :: ocnfrac (pcols)           ! Ocean fraction (fraction)
    real(r8), intent(in) :: snowh(pcols)              ! Snow depth over land, water equivalent (m)
    real(r8), intent(in) :: dlf(pcols,pver)           ! Detrained water from convection schemes
    real(r8), intent(in) :: dlf2(pcols,pver)          ! Detrained water from shallow convection scheme
    real(r8), intent(in) :: cmfmc(pcols,pverp)        ! Deep + Shallow Convective mass flux [ kg /s/m^2 ]
                                                 
    real(r8), intent(in) :: ts(pcols)                 ! Surface temperature
    real(r8), intent(in) :: sst(pcols)                ! Sea surface temperature
    real(r8), intent(in) :: zdu(pcols,pver)           ! Detrainment rate from deep convection

    ! These two variables are needed for energy check
    real(r8), intent(out) :: det_s(pcols)             ! Integral of detrained static energy from ice
    real(r8), intent(out) :: det_ice(pcols)           ! Integral of detrained ice for energy check

    type(physics_state) state_loc                     ! Local copy of the state variable
    type(physics_ptend) ptend_loc                     ! Local parameterization tendencies

    integer i, k, lchnk, ncol

    integer itim_old
    real(r8), pointer, dimension(:,:) :: qcwat        ! Cloud water old q
    real(r8), pointer, dimension(:,:) :: tcwat        ! Cloud water old temperature
    real(r8), pointer, dimension(:,:) :: lcwat        ! Cloud liquid water old q
    real(r8), pointer, dimension(:,:) :: iccwat       ! Cloud ice water old q
    real(r8), pointer, dimension(:,:) :: nlwat        ! Cloud liquid droplet number condentration. old.
    real(r8), pointer, dimension(:,:) :: niwat        ! Cloud ice    droplet number condentration. old.
    real(r8), pointer, dimension(:,:) :: cc_T         ! Grid-mean microphysical tendency
    real(r8), pointer, dimension(:,:) :: cc_qv        ! Grid-mean microphysical tendency
    real(r8), pointer, dimension(:,:) :: cc_ql        ! Grid-mean microphysical tendency
    real(r8), pointer, dimension(:,:) :: cc_qi        ! Grid-mean microphysical tendency
    real(r8), pointer, dimension(:,:) :: cc_nl        ! Grid-mean microphysical tendency
    real(r8), pointer, dimension(:,:) :: cc_ni        ! Grid-mean microphysical tendency
    real(r8), pointer, dimension(:,:) :: cc_qlst      ! In-liquid stratus microphysical tendency
    real(r8), pointer, dimension(:,:) :: cld          ! Total cloud fraction
    real(r8), pointer, dimension(:,:) :: ast          ! Relative humidity cloud fraction
    real(r8), pointer, dimension(:,:) :: aist         ! Physical ice stratus fraction
    real(r8), pointer, dimension(:,:) :: alst         ! Physical liquid stratus fraction
    real(r8), pointer, dimension(:,:) :: qist         ! Physical in-cloud IWC
    real(r8), pointer, dimension(:,:) :: qlst         ! Physical in-cloud LWC
    real(r8), pointer, dimension(:,:) :: concld       ! Convective cloud fraction
    real(r8), pointer, dimension(:,:) :: shfrc        ! Cloud fraction from shallow convection scheme
    real(r8), pointer, dimension(:,:) :: cmfmc_sh     ! Shallow convective mass flux (pcols,pverp) [ kg/s/m^2 ]
    real(r8), pointer, dimension(:,:) :: cmeliq
    real(r8), pointer, dimension(:,:) :: tke
    real(r8), pointer, dimension(:,:) :: qtl_flx
    real(r8), pointer, dimension(:,:) :: qti_flx
    real(r8), pointer, dimension(:,:) :: cmfr_det
    real(r8), pointer, dimension(:,:) :: qlr_det
    real(r8), pointer, dimension(:,:) :: qir_det
    ! Convective cloud to the physics buffer for purposes of ql contrib. to radn.
    real(r8), pointer, dimension(:,:) :: fice_ql      ! Cloud ice/water partitioning ratio.

    ! ZM microphysics
    real(r8), pointer :: dlfzm (:,:) ! ZM detrained convective cloud water mixing ratio.
    real(r8), pointer :: difzm (:,:) ! ZM detrained convective cloud ice mixing ratio.
    real(r8), pointer :: dnlfzm(:,:) ! ZM detrained convective cloud water num concen.
    real(r8), pointer :: dnifzm(:,:) ! ZM detrained convective cloud ice num concen.

    real(r8) latsub

    ! Tendencies for ice saturation adjustment
    real(r8) stend  (pcols,pver)
    real(r8) qvtend (pcols,pver)
    real(r8) qitend (pcols,pver)
    real(r8) initend(pcols,pver)

    ! Local variables for cldfrc

    real(r8) cldst  (pcols,pver)                     ! Stratus cloud fraction
    real(r8) rhcloud(pcols,pver)                     ! Relative humidity cloud (last timestep)
    real(r8) clc    (pcols)                          ! Column convective cloud amount
    real(r8) rhu00  (pcols,pver)                     ! RH threshold for cloud
    real(r8) icecldf(pcols,pver)                     ! Ice cloud fraction
    real(r8) liqcldf(pcols,pver)                     ! Liquid cloud fraction (combined into cloud)
    real(r8) relhum (pcols,pver)                     ! RH, output to determine drh/da

    ! Local variables for macrophysics

    real(r8) rdtime                                  ! 1./dtime
    real(r8) qtend   (pcols,pver)                    ! Moisture tendencies
    real(r8) ttend   (pcols,pver)                    ! Temperature tendencies
    real(r8) ltend   (pcols,pver)                    ! Cloud liquid water tendencies
    real(r8) fice    (pcols,pver)                    ! Fractional ice content within cloud
    real(r8) fsnow   (pcols,pver)                    ! Fractional snow production
    real(r8) homoo   (pcols,pver)
    real(r8) qcreso  (pcols,pver)
    real(r8) prcio   (pcols,pver)
    real(r8) praio   (pcols,pver)
    real(r8) qireso  (pcols,pver)
    real(r8) ftem    (pcols,pver)
    real(r8) pracso  (pcols,pver)
    real(r8) dpdlfliq(pcols,pver)
    real(r8) dpdlfice(pcols,pver)
    real(r8) shdlfliq(pcols,pver)
    real(r8) shdlfice(pcols,pver)
    real(r8) dpdlft  (pcols,pver)
    real(r8) shdlft  (pcols,pver)

    real(r8) dum1
    real(r8) qc(pcols,pver)
    real(r8) qi(pcols,pver)
    real(r8) nc(pcols,pver)
    real(r8) ni(pcols,pver)

    logical  lq(pcnst)

    ! Output from mmacro_pcond

    real(r8) tlat (pcols,pver)
    real(r8) qvlat(pcols,pver)
    real(r8) qcten(pcols,pver)
    real(r8) qiten(pcols,pver)
    real(r8) ncten(pcols,pver)
    real(r8) niten(pcols,pver)

    ! Output from mmacro_pcond

    real(r8) qvadj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (vapor)
    real(r8) qladj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (liquid)
    real(r8) qiadj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (ice)
    real(r8) qllim(pcols,pver)                       ! Macro-physics tendency from "instratus_condensate" call (liquid)
    real(r8) qilim(pcols,pver)                       ! Macro-physics tendency from "instratus_condensate" call (ice)

    ! For revised macophysics, mmacro_pcond

    real(r8) itend     (pcols,pver)
    real(r8) lmitend   (pcols,pver)
    real(r8) zeros     (pcols,pver)
    real(r8) t_inout   (pcols,pver)
    real(r8) qv_inout  (pcols,pver)
    real(r8) ql_inout  (pcols,pver)
    real(r8) qi_inout  (pcols,pver)
    real(r8) concld_old(pcols,pver)

    ! Note that below 'clr_old' is defined using 'alst_old' not 'ast_old' for full consistency with the
    ! liquid condensation process which is using 'alst' not 'ast'.
    ! For microconsistency use 'concld_old', since 'alst_old' was computed using 'concld_old'.
    ! Since convective updraft fractional area is small, it does not matter whether 'concld' or 'concld_old' is used.
    ! Note also that 'clri_old' is defined using 'ast_old' since current microphysics is operating on 'ast_old'
    real(r8) clrw_old(pcols,pver) ! (1 - concld_old - alst_old)
    real(r8) clri_old(pcols,pver) ! (1 - concld_old -  ast_old)

    real(r8) nl_inout(pcols,pver)
    real(r8) ni_inout(pcols,pver)

    real(r8) nltend(pcols,pver)
    real(r8) nitend(pcols,pver)

    ! For detraining cumulus condensate into the 'stratus' without evaporation
    ! This is for use in mmacro_pcond

    real(r8) dlf_T (pcols,pver)
    real(r8) dlf_qv(pcols,pver)
    real(r8) dlf_ql(pcols,pver)
    real(r8) dlf_qi(pcols,pver)
    real(r8) dlf_nl(pcols,pver)
    real(r8) dlf_ni(pcols,pver)

    ! Local variables for CFMIP calculations
    real(r8) mr_lsliq(pcols,pver)  ! mixing_ratio_large_scale_cloud_liquid (kg/kg)
    real(r8) mr_lsice(pcols,pver)  ! mixing_ratio_large_scale_cloud_ice (kg/kg)
    real(r8) mr_ccliq(pcols,pver)  ! mixing_ratio_convective_cloud_liquid (kg/kg)
    real(r8) mr_ccice(pcols,pver)  ! mixing_ratio_convective_cloud_ice (kg/kg)

    ! CloudSat equivalent ice mass mixing ratio (kg/kg)
    real(r8) cldsice(pcols,pver)

    lchnk = state%lchnk
    ncol  = state%ncol

    call physics_state_copy(state, state_loc)

    ! Associate pointers with physics buffer fields

    itim_old = pbuf_old_tim_idx()

    call pbuf_get_field(pbuf, qcwat_idx  , qcwat  , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, tcwat_idx  , tcwat  , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, lcwat_idx  , lcwat  , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, iccwat_idx , iccwat , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, nlwat_idx  , nlwat  , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, niwat_idx  , niwat  , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, cc_t_idx   , cc_t   , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, cc_qv_idx  , cc_qv  , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, cc_ql_idx  , cc_ql  , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, cc_qi_idx  , cc_qi  , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, cc_nl_idx  , cc_nl  , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, cc_ni_idx  , cc_ni  , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, cc_qlst_idx, cc_qlst, start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, cld_idx    , cld    , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, concld_idx , concld , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, ast_idx    , ast    , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, aist_idx   , aist   , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, alst_idx   , alst   , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, qist_idx   , qist   , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, qlst_idx   , qlst   , start=[1,1,itim_old], kount=[pcols,pver,1])
    call pbuf_get_field(pbuf, cmeliq_idx , cmeliq)

    ! For purposes of convective ql.

    call pbuf_get_field(pbuf, fice_idx,     fice_ql )

    call pbuf_get_field(pbuf, cmfmc_sh_idx, cmfmc_sh)

    ! Check that qcwat and tcwat were initialized; if not then do it now.
    if (qcwat(1,1) == huge(1.0_r8)) then
      qcwat(:ncol,:) = state%q(:ncol,:,1)
    end if
    if (tcwat(1,1) == huge(1.0_r8)) then
      tcwat(:ncol,:) = state%t(:ncol,:)
    end if

    ! Initialize convective detrainment tendency

    dlf_T  = 0
    dlf_qv = 0
    dlf_ql = 0
    dlf_qi = 0
    dlf_nl = 0
    dlf_ni = 0

    ! ------------------------------------- !
    ! From here, process computation begins !
    ! ------------------------------------- !

    ! ----------------------------------------------------------------------------- !
    ! Detrainment of convective condensate into the environment or stratiform cloud !
    ! ----------------------------------------------------------------------------- !

    lq(:)        = .false.
    lq(ixcldliq) = .true.
    lq(ixcldice) = .true.
    lq(ixnumliq) = .true.
    lq(ixnumice) = .true.
    call physics_ptend_init(ptend_loc, state%psetcols, 'pcwdetrain', ls=.true., lq=lq)

    ! Procedures :
    ! (1) Partition detrained convective cloud water into liquid and ice based on T.
    !     This also involves heating.
    !     If convection scheme can handle this internally, this step is not necssary.
    ! (2) Assuming a certain effective droplet radius, computes number concentration
    !     of detrained convective cloud liquid and ice.
    ! (3) If 'cu_det_st = .true' ('false'), detrain convective cloud 'liquid' into
    !     the pre-existing 'liquid' stratus ( mean environment ).  The former does
    !     not involve any macrophysical evaporation while the latter does. This is
    !     a kind of 'targetted' deposition. Then, force in-stratus LWC to be bounded
    !     by qcst_min and qcst_max in mmacro_pcond.
    ! (4) In contrast to liquid, convective ice is detrained into the environment
    !     and involved in the sublimation. Similar bounds as liquid stratus are imposed.
    ! This is the key procesure generating upper-level cirrus clouds.
    ! The unit of dlf : [ kg/kg/s ]

    if (zmconv_microp) then
      call pbuf_get_field(pbuf, dlfzm_idx , dlfzm )
      call pbuf_get_field(pbuf, difzm_idx , difzm )
      call pbuf_get_field(pbuf, dnlfzm_idx, dnlfzm)
      call pbuf_get_field(pbuf, dnifzm_idx, dnifzm)
    end if

    det_s    = 0
    det_ice  = 0
    dpdlfliq = 0
    dpdlfice = 0
    shdlfliq = 0
    shdlfice = 0
    dpdlft   = 0
    shdlft   = 0

    do k = top_lev, pver
      do i = 1, state_loc%ncol
        if (state_loc%t(i,k) > 268.15_r8) then
          dum1 = 0
        else if (state_loc%t(i,k) < 238.15_r8) then
          dum1 = 1
        else
          dum1 = (268.15_r8 - state_loc%t(i,k)) / 30.0_r8
        end if

        ! If detrainment was done elsewhere, still update the variables used for output
        ! assuming that the temperature split between liquid and ice is the same as assumed
        ! here.
        if (zmconv_microp) then
          ptend_loc%q(i,k,ixcldliq) = dlfzm(i,k) + dlf2(i,k) * (1 - dum1)
          ptend_loc%q(i,k,ixcldice) = difzm(i,k) + dlf2(i,k) * dum1

          ptend_loc%q(i,k,ixnumliq) = dnlfzm(i,k) + 3 * (dlf2(i,k) * (1 - dum1)) &
                                                  / (4 * 3.14_r8 * 10.0e-6_r8**3 * 997.0_r8)   ! Shallow Convection
          ptend_loc%q(i,k,ixnumice) = dnifzm(i,k) + 3 * (dlf2(i,k) * dum1) &
                                                  / (4 * 3.14_r8 * 50.0e-6_r8**3 * 500.0_r8)   ! Shallow Convection
          ptend_loc%s(i,k)          = dlf2(i,k) * dum1 * latice
        else if (do_detrain) then
          ptend_loc%q(i,k,ixcldliq) = dlf(i,k) * (1 - dum1)
          ptend_loc%q(i,k,ixcldice) = dlf(i,k) * dum1
          ptend_loc%q(i,k,ixnumliq) = 3 * (max(0.0_r8, (dlf(i,k) - dlf2(i,k))) * (1 - dum1)) / &
                                      (4 * 3.14_r8 *  8.0e-6_r8**3 * 997.0_r8) + & ! Deep    Convection
                                      3 * (                        dlf2(i,k)   * (1 - dum1)) / &
                                      (4 * 3.14_r8 * 10.0e-6_r8**3 * 997.0_r8)     ! Shallow Convection
          ptend_loc%q(i,k,ixnumice) = 3 * (max(0.0_r8, (dlf(i,k) - dlf2(i,k))) * dum1) / &
                                      (4 * 3.14_r8 * 25.0e-6_r8**3 * 500.0_r8) + & ! Deep    Convection
                                      3 * (                        dlf2(i,k)   * dum1) / &
                                      (4 * 3.14_r8 * 50.0e-6_r8**3 * 500.0_r8)     ! Shallow Convection
          ptend_loc%s(i,k)          = dlf(i,k) * dum1 * latice
        else
          ptend_loc%q(i,k,ixcldliq) = 0
          ptend_loc%q(i,k,ixcldice) = 0
          ptend_loc%q(i,k,ixnumliq) = 0
          ptend_loc%q(i,k,ixnumice) = 0
          ptend_loc%s(i,k)          = 0
        end if

        ! Only rliq is saved from deep convection, which is the reserved liquid.  We need to keep
        !   track of the integrals of ice and static energy that is effected from conversion to ice
        !   so that the energy checker doesn't complain.
        det_s  (i) = det_s  (i) + ptend_loc%s(i,k) * state_loc%pdel(i,k) / gravit
        det_ice(i) = det_ice(i) - ptend_loc%q(i,k,ixcldice) * state_loc%pdel(i,k) / gravit

        ! Targetted detrainment of convective liquid water either directly into the
        ! existing liquid stratus or into the environment.
        if (cu_det_st) then
          dlf_T (i,k) = ptend_loc%s(i,k) / cpair
          dlf_qv(i,k) = 0
          dlf_ql(i,k) = ptend_loc%q(i,k,ixcldliq)
          dlf_qi(i,k) = ptend_loc%q(i,k,ixcldice)
          dlf_nl(i,k) = ptend_loc%q(i,k,ixnumliq)
          dlf_ni(i,k) = ptend_loc%q(i,k,ixnumice)
          ptend_loc%q(i,k,ixcldliq) = 0
          ptend_loc%q(i,k,ixcldice) = 0
          ptend_loc%q(i,k,ixnumliq) = 0
          ptend_loc%q(i,k,ixnumice) = 0
          ptend_loc%s(i,k)          = 0
          dpdlfliq(i,k)             = 0
          dpdlfice(i,k)             = 0
          shdlfliq(i,k)             = 0
          shdlfice(i,k)             = 0
          dpdlft  (i,k)             = 0
          shdlft  (i,k)             = 0
        else
          if (zmconv_microp) then
            dpdlfliq(i,k) = dlfzm(i,k)
            dpdlfice(i,k) = difzm(i,k)
            dpdlft  (i,k) = 0
          else
            dpdlfliq(i,k) = (dlf(i,k) - dlf2(i,k)) * (1 - dum1)
            dpdlfice(i,k) = (dlf(i,k) - dlf2(i,k)) * dum1
            dpdlft  (i,k) = (dlf(i,k) - dlf2(i,k)) * dum1 * latice / cpair
          end if

          shdlfliq(i,k) = dlf2(i,k) * (1 - dum1)
          shdlfice(i,k) = dlf2(i,k) * dum1
          shdlft  (i,k) = dlf2(i,k) * dum1 * latice / cpair
        end if
      end do
    end do

    ! Divide by density of water
    det_ice(:ncol) = det_ice(:ncol) / 1000.0_r8

    ! Add the detrainment tendency to the output tendency
    call physics_ptend_init(ptend, state%psetcols, 'macrop')
    call physics_ptend_sum(ptend_loc, ptend, ncol)

    ! update local copy of state with the detrainment tendency
    ! ptend_loc is reset to zero by this call
    call physics_update(state_loc, ptend_loc, dtime)

    ! -------------------------------------- !
    ! Computation of Various Cloud Fractions !
    ! -------------------------------------- !

    ! ----------------------------------------------------------------------------- !
    ! Treatment of cloud fraction in CAM4 and CAM5 differs                          !
    ! (1) CAM4                                                                      !
    !     . Cumulus AMT = Deep    Cumulus AMT ( empirical fcn of mass flux ) +      !
    !                     Shallow Cumulus AMT ( empirical fcn of mass flux )        !
    !     . Stratus AMT = max( RH stratus AMT, Stability Stratus AMT )              !
    !     . Cumulus and Stratus are 'minimally' overlapped without hierarchy.       !
    !     . Cumulus LWC,IWC is assumed to be the same as Stratus LWC,IWC            !
    ! (2) CAM5                                                                      !
    !     . Cumulus AMT = Deep    Cumulus AMT ( empirical fcn of mass flux ) +      !
    !                     Shallow Cumulus AMT ( internally fcn of mass flux and w ) !
    !     . Stratus AMT = fcn of environmental-mean RH ( no Stability Stratus )     !
    !     . Cumulus and Stratus are non-overlapped with higher priority on Cumulus  !
    !     . Cumulus ( both Deep and Shallow ) has its own LWC and IWC.              !
    ! ----------------------------------------------------------------------------- !

    concld_old(:ncol,top_lev:pver) = concld(:ncol,top_lev:pver)

    nullify(tke, qtl_flx, qti_flx, cmfr_det, qlr_det, qir_det)
    if (tke_idx      > 0) call pbuf_get_field(pbuf, tke_idx, tke)
    if (qtl_flx_idx  > 0) call pbuf_get_field(pbuf, qtl_flx_idx,  qtl_flx)
    if (qti_flx_idx  > 0) call pbuf_get_field(pbuf, qti_flx_idx,  qti_flx)
    if (cmfr_det_idx > 0) call pbuf_get_field(pbuf, cmfr_det_idx, cmfr_det)
    if (qlr_det_idx  > 0) call pbuf_get_field(pbuf, qlr_det_idx,  qlr_det)
    if (qir_det_idx  > 0) call pbuf_get_field(pbuf, qir_det_idx,  qir_det)

    clrw_old(:ncol,:top_lev-1) = 0
    clri_old(:ncol,:top_lev-1) = 0
    do k = top_lev, pver
      do i = 1, ncol
        clrw_old(i,k) = max(0.0_r8, min(1.0_r8, 1.0_r8 - concld(i,k) - alst(i,k)))
        clri_old(i,k) = max(0.0_r8, min(1.0_r8, 1.0_r8 - concld(i,k) -  ast(i,k)))
      end do
    end do

    if (use_shfrc) then
      call pbuf_get_field(pbuf, shfrc_idx, shfrc)
    else
      allocate(shfrc(pcols,pver))
      shfrc = 0
    end if

    ! CAM5 only uses 'concld' output from the below subroutine.
    ! Stratus ('ast' = max(alst,aist)) and total cloud fraction ('cld = ast + concld')
    ! will be computed using this updated 'concld' in the stratiform macrophysics
    ! scheme (mmacro_pcond) later below.

    call cldfrc(                              &
      lchnk       =lchnk                    , &
      ncol        =ncol                     , &
      pbuf        =pbuf                     , &
      pmid        =state_loc%pmid           , &
      temp        =state_loc%t              , &
      q           =state_loc%q(:,:,1)       , &
      omga        =state_loc%omega          , &
      phis        =state_loc%phis           , &
      shfrc       =shfrc                    , &
      use_shfrc   =use_shfrc                , &
      cloud       =cld                      , &
      rhcloud     =rhcloud                  , &
      clc         =clc                      , &
      pdel        =state_loc%pdel           , &
      cmfmc       =cmfmc                    , &
      cmfmc2      =cmfmc_sh                 , &
      landfrac    =landfrac                 , &
      snowh       =snowh                    , &
      concld      =concld                   , &
      cldst       =cldst                    , &
      ts          =ts                       , &
      sst         =sst                      , &
      ps          =state_loc%pint(:,pverp)  , &
      zdu         =zdu                      , &
      ocnfrac     =ocnfrac                  , &
      rhu00       =rhu00                    , &
      cldice      =state_loc%q(:,:,ixcldice), &
      icecldf     =icecldf                  , &
      liqcldf     =liqcldf                  , &
      relhum      =relhum                   , &
      dindex      =0                        )

    ! ---------------------------------------------- !
    ! Stratiform Cloud Macrophysics and Microphysics !
    ! ---------------------------------------------- !

    lchnk  = state_loc%lchnk
    ncol   = state_loc%ncol
    rdtime = 1.0_r8 / dtime

    ! Define fractional amount of stratus condensate and precipitation in ice phase.
    ! This uses a ramp ( -30 ~ -10 for fice, -5 ~ 0 for fsnow ).
    ! The ramp within convective cloud may be different

    call cldfrc_fice(ncol, state_loc%t, fice, fsnow)

    lq(:)        = .false.
    lq(1)        = .true.
    lq(ixcldice) = .true.
    lq(ixcldliq) = .true.
    lq(ixnumliq) = .true.
    lq(ixnumice) = .true.
 
    ! Initialize local physics_ptend object again
    call physics_ptend_init(ptend_loc, state%psetcols, 'macro_park', ls=.true., lq=lq )

    ! --------------------------------- !
    ! Liquid Macrop_Driver Macrophysics !
    ! --------------------------------- !

    zeros(:ncol,top_lev:pver)  = 0
    qc(:ncol,top_lev:pver) = state_loc%q(:ncol,top_lev:pver,ixcldliq)
    qi(:ncol,top_lev:pver) = state_loc%q(:ncol,top_lev:pver,ixcldice)
    nc(:ncol,top_lev:pver) = state_loc%q(:ncol,top_lev:pver,ixnumliq)
    ni(:ncol,top_lev:pver) = state_loc%q(:ncol,top_lev:pver,ixnumice)

    ! In CAM5, 'microphysical forcing' ( cc_... ) and 'the other advective forcings' ( ttend, ... )
    ! are separately provided into the prognostic microp_driver macrophysics scheme. This is an
    ! attempt to resolve in-cloud and out-cloud forcings.

    if (get_nstep() <= 1) then
      tcwat  (:ncol,top_lev:pver) = state_loc%t(:ncol,top_lev:pver)
      qcwat  (:ncol,top_lev:pver) = state_loc%q(:ncol,top_lev:pver,1)
      lcwat  (:ncol,top_lev:pver) = qc(:ncol,top_lev:pver) + qi(:ncol,top_lev:pver)
      iccwat (:ncol,top_lev:pver) = qi(:ncol,top_lev:pver)
      nlwat  (:ncol,top_lev:pver) = nc(:ncol,top_lev:pver)
      niwat  (:ncol,top_lev:pver) = ni(:ncol,top_lev:pver)
      ttend  (:ncol,:) = 0
      qtend  (:ncol,:) = 0
      ltend  (:ncol,:) = 0
      itend  (:ncol,:) = 0
      nltend (:ncol,:) = 0
      nitend (:ncol,:) = 0
      cc_T   (:ncol,:) = 0
      cc_qv  (:ncol,:) = 0
      cc_ql  (:ncol,:) = 0
      cc_qi  (:ncol,:) = 0
      cc_nl  (:ncol,:) = 0
      cc_ni  (:ncol,:) = 0
      cc_qlst(:ncol,:) = 0
    else
      ttend(:ncol,top_lev:pver)  = (state_loc%t(:ncol,top_lev:pver)   -  tcwat(:ncol,top_lev:pver)) * rdtime &
                                 - cc_t(:ncol,top_lev:pver)
      qtend(:ncol,top_lev:pver)  = (state_loc%q(:ncol,top_lev:pver,1) -  qcwat(:ncol,top_lev:pver)) * rdtime &
                                 - cc_qv(:ncol,top_lev:pver)
      ltend(:ncol,top_lev:pver)  = (qc(:ncol,top_lev:pver) + qi(:ncol,top_lev:pver) - lcwat(:ncol,top_lev:pver)) * rdtime &
                                 - (cc_ql(:ncol,top_lev:pver) + cc_qi(:ncol,top_lev:pver))
      itend(:ncol,top_lev:pver)  = (qi(:ncol,top_lev:pver)   - iccwat(:ncol,top_lev:pver)) * rdtime &
                                 - cc_qi(:ncol,top_lev:pver)
      nltend(:ncol,top_lev:pver) = (nc(:ncol,top_lev:pver)   - nlwat(:ncol,top_lev:pver)) * rdtime &
                                 - cc_nl(:ncol,top_lev:pver)
      nitend(:ncol,top_lev:pver) = (ni(:ncol,top_lev:pver)   - niwat(:ncol,top_lev:pver)) * rdtime &
                                 - cc_ni(:ncol,top_lev:pver)
    end if
    lmitend (:ncol,top_lev:pver) = ltend (:ncol,top_lev:pver) - itend(:ncol,top_lev:pver)
    t_inout (:ncol,top_lev:pver) = tcwat (:ncol,top_lev:pver)
    qv_inout(:ncol,top_lev:pver) = qcwat (:ncol,top_lev:pver)
    ql_inout(:ncol,top_lev:pver) = lcwat (:ncol,top_lev:pver) - iccwat(:ncol,top_lev:pver)
    qi_inout(:ncol,top_lev:pver) = iccwat(:ncol,top_lev:pver)
    nl_inout(:ncol,top_lev:pver) = nlwat (:ncol,top_lev:pver)
    ni_inout(:ncol,top_lev:pver) = niwat (:ncol,top_lev:pver)

    ! Liquid Microp_Driver Macrophysics.
    ! The main roles of this subroutines are
    ! (1) compute net condensation rate of stratiform liquid ( cmeliq )
    ! (2) compute liquid stratus and ice stratus fractions.
    ! Note 'ttend...' are advective tendencies except microphysical process while
    !      'CC...'    are microphysical tendencies.

    call mmacro_pcond( &
      lchnk       =lchnk          , &
      ncol        =ncol           , &
      dt          =dtime          , &
      p           =state_loc%pmid , &
      dp          =state_loc%pdel , &
      t0          =t_inout        , &
      qv0         =qv_inout       , &
      ql0         =ql_inout       , &
      qi0         =qi_inout       , &
      nl0         =nl_inout       , &
      ni0         =ni_inout       , &
      a_t         =ttend          , &
      a_qv        =qtend          , &
      a_ql        =lmitend        , &
      a_qi        =itend          , &
      a_nl        =nltend         , &
      a_ni        =nitend         , &
      c_t         =cc_t           , &
      c_qv        =cc_qv          , &
      c_ql        =cc_ql          , &
      c_qi        =cc_qi          , &
      c_nl        =cc_nl          , &
      c_ni        =cc_ni          , &
      c_qlst      =cc_qlst        , &
      d_t         =dlf_t          , &
      d_qv        =dlf_qv         , &
      d_ql        =dlf_ql         , &
      d_qi        =dlf_qi         , &
      d_nl        =dlf_nl         , &
      d_ni        =dlf_ni         , &
      a_cud       =concld_old     , &
      a_cu0       =concld         , &
      clrw_old    =clrw_old       , &
      clri_old    =clri_old       , &
      landfrac    =landfrac       , &
      snowh       =snowh          , &
      tke         =tke            , &
      qtl_flx     =qtl_flx        , &
      qti_flx     =qti_flx        , &
      cmfr_det    =cmfr_det       , &
      qlr_det     =qlr_det        , &
      qir_det     =qir_det        , &
      s_tendout   =tlat           , &
      qv_tendout  =qvlat          , &
      ql_tendout  =qcten          , &
      qi_tendout  =qiten          , &
      nl_tendout  =ncten          , &
      ni_tendout  =niten          , &
      qme         =cmeliq         , &
      qvadj       =qvadj          , &
      qladj       =qladj          , &
      qiadj       =qiadj          , &
      qllim       =qllim          , &
      qilim       =qilim          , &
      cld         =cld            , &
      al_st_star  =alst           , &
      ai_st_star  =aist           , &
      ql_st_star  =qlst           , &
      qi_st_star  =qist           , &
      do_cldice   =do_cldice      )

    ! Copy of concld/fice to put in physics buffer
    ! Below are used only for convective cloud.

    fice_ql(:ncol,:top_lev-1)   = 0
    fice_ql(:ncol,top_lev:pver) = fice(:ncol,top_lev:pver)

    ! Compute net stratus fraction using maximum over-lapping assumption
    ast(:ncol,:top_lev-1)   = 0
    ast(:ncol,top_lev:pver) = max(alst(:ncol,top_lev:pver), aist(:ncol,top_lev:pver))

    do k = top_lev, pver
      do i = 1, ncol
        ptend_loc%s(i,k)          =  tlat(i,k)
        ptend_loc%q(i,k,1)        = qvlat(i,k)
        ptend_loc%q(i,k,ixcldliq) = qcten(i,k)
        ptend_loc%q(i,k,ixcldice) = qiten(i,k)
        ptend_loc%q(i,k,ixnumliq) = ncten(i,k)
        ptend_loc%q(i,k,ixnumice) = niten(i,k)

        ! Check to make sure that the macrophysics code is respecting the flags that control
        ! whether cldwat should be prognosing cloud ice and cloud liquid or not.
        if (.not. do_cldice .and. qiten(i,k) /= 0) then
          call endrun('macrop_driver: ERROR - ' // &
            'Cldwat is configured not to prognose cloud ice, but mmacro_pcond has ice mass tendencies.')
        end if
        if (.not. do_cldice .and. niten(i,k) /= 0) then
          call endrun('macrop_driver: ERROR - ' // &
            'Cldwat is configured not to prognose cloud ice, but mmacro_pcond has ice number tendencies.')
        end if

        if (.not. do_cldliq .and. qcten(i,k) /= 0) then
          call endrun('macrop_driver: ERROR - ' // &
            'Cldwat is configured not to prognose cloud liquid, but mmacro_pcond has liquid mass tendencies.')
        end if
        if (.not. do_cldliq .and. ncten(i,k) /= 0) then
          call endrun('macrop_driver: ERROR - ' // &
            'Cldwat is configured not to prognose cloud liquid, but mmacro_pcond has liquid number tendencies.')
        end if
      end do
    end do

    ! Update the output tendencies with the mmacro_pcond tendencies
    call physics_ptend_sum(ptend_loc, ptend, ncol)

    ! state_loc is the equlibrium state after macrophysics
    call physics_update(state_loc, ptend_loc, dtime)

    ! Calculations and outfld calls for CLDLIQSTR, CLDICESTR, CLDLIQCON, CLDICECON for CFMIP

    ! Initialize local variables
    mr_ccliq = 0   ! not seen by radiation, so setting to 0
    mr_ccice = 0   ! not seen by radiation, so setting to 0
    mr_lsliq = 0
    mr_lsice = 0

    do k = top_lev, pver
      do i = 1, ncol
        if (cld(i,k) > 0) then
          mr_lsliq(i,k) = state_loc%q(i,k,ixcldliq)
          mr_lsice(i,k) = state_loc%q(i,k,ixcldice)
        else
          mr_lsliq(i,k) = 0
          mr_lsice(i,k) = 0
        end if
      end do
    end do

    ! ------------------------------------------------- !
    ! Save equilibrium state variables for macrophysics !
    ! at the next time step                             !
    ! ------------------------------------------------- !
    cldsice = 0
    do k = top_lev, pver
      tcwat  (:ncol,k) = state_loc%t(:ncol,k)
      qcwat  (:ncol,k) = state_loc%q(:ncol,k,1)
      lcwat  (:ncol,k) = state_loc%q(:ncol,k,ixcldliq) + state_loc%q(:ncol,k,ixcldice)
      iccwat (:ncol,k) = state_loc%q(:ncol,k,ixcldice)
      nlwat  (:ncol,k) = state_loc%q(:ncol,k,ixnumliq)
      niwat  (:ncol,k) = state_loc%q(:ncol,k,ixnumice)
      cldsice(:ncol,k) = lcwat(:ncol,k) * min(1.0_r8, max(0.0_r8, (tmelt - tcwat(:ncol,k)) / 20.0_r8))
    end do

    call physics_state_dealloc(state_loc)

  end subroutine macrop_driver_tend


  ! Saturation adjustment for liquid
  !
  ! With CLUBB, we are seeing relative humidity with respect to water
  ! greater than 1. This should not be happening and is not what the
  ! microphsyics expects from the macrophysics. As a work around while
  ! this issue is investigated in CLUBB, this routine will enfornce a
  ! maximum RHliq of 1 everywhere in the atmosphere. Any excess water will
  ! be converted into cloud drops.
  elemental subroutine liquid_macro_tend(npccn, t, p, qv, qc, nc, xxlv, deltat, stend, qvtend, qctend, nctend)

    use wv_sat_methods, only: wv_sat_qsat_ice, wv_sat_qsat_water
    use micro_mg_utils, only: rhow
    use physconst     , only: rair
    use cldfrc2m      , only: rhmini_const, rhmaxi_const

    real(r8), intent(in)  :: npccn  !Activated number of cloud condensation nuclei
    real(r8), intent(in)  :: t      !temperature (k)
    real(r8), intent(in)  :: p      !pressure (pa)
    real(r8), intent(in)  :: qv     !water vapor mixing ratio
    real(r8), intent(in)  :: qc     !liquid mixing ratio
    real(r8), intent(in)  :: nc     !liquid number concentration
    real(r8), intent(in)  :: xxlv   !latent heat of vaporization
    real(r8), intent(in)  :: deltat !timestep
    real(r8), intent(out) :: stend  ! 'temperature' tendency
    real(r8), intent(out) :: qvtend !vapor tendency
    real(r8), intent(out) :: qctend !liquid mass tendency
    real(r8), intent(out) :: nctend !liquid number tendency

    real(r8) ESL
    real(r8) QSL

    stend  = 0
    qvtend = 0
    qctend = 0
    nctend = 0

    ! Calculate qsatl from t, p, q
    call wv_sat_qsat_water(t, p, ESL, QSL)

    ! Don't allow supersaturation with respect to liquid.
    if (qv > QSL) then
      qctend = (qv - QSL) / deltat
      qvtend = 0.0_r8 - qctend
      stend  = qctend * xxlv  ! Moist static energy tend...[J/kg/s] !

      ! If drops  exists (more than 1 L-1) and there is condensation,
      ! do not add to number (= growth), otherwise  add 6um drops.
      !
      ! This is somewhat arbitrary, but ensures that some reasonable droplet
      ! size is create to remove the excess water. This could be enhanced to
      ! look at npccn, but ideally this entire routine should go away.
      if (nc * p / rair / t < 1e3_r8 .and. (qc + qctend * deltat) > 1e-18_r8) then
        nctend = nctend + 3 * qctend / (4 * 3.14_r8 * 6.0e-6_r8**3 * rhow)
      end if
    end if

  end subroutine liquid_macro_tend

end module macrop_driver
