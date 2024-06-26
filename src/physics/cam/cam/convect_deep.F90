
module convect_deep

  !---------------------------------------------------------------------------------
  ! Purpose:
  !
  ! CAM interface to several deep convection interfaces. Currently includes:
  !    Zhang-McFarlane (default)
  !    Kerry Emanuel 
  !
  !
  ! Author: D.B. Coleman, Sep 2004
  !
  !---------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid      , only: pver, pcols, pverp
  use cam_logfile , only: iulog

  implicit none
  private
  save

  public convect_deep_register
  public convect_deep_init
  public convect_deep_tend
  public convect_deep_tend_2
  public deep_scheme_does_scav_trans
   
  character(16) deep_scheme

  integer :: icwmrdp_idx     = 0 
  integer :: rprddp_idx      = 0 
  integer :: nevapr_dpcu_idx = 0 
  integer :: cldtop_idx      = 0 
  integer :: cldbot_idx      = 0 
  integer :: cld_idx         = 0 
  integer :: fracis_idx      = 0 
  integer :: pblh_idx        = 0 
  integer :: tpert_idx       = 0 
  integer :: prec_dp_idx     = 0
  integer :: snow_dp_idx     = 0
  integer :: ttend_dp_idx    = 0

contains 

  pure logical function deep_scheme_does_scav_trans() result(res)
    !
    ! Function called by tphysbc to determine if it needs to do scavenging and convective transport
    ! or if those have been done by the deep convection scheme. Each scheme could have its own
    ! identical query function for a less-knowledgable interface but for now, we know that KE 
    ! does scavenging & transport, and ZM doesn't
    !

    res = deep_scheme == 'KE'

  end function deep_scheme_does_scav_trans

  subroutine convect_deep_register

    use physics_buffer, only: pbuf_add_field, dtype_r8
    use zm_conv_intr  , only: zm_conv_register
    use phys_control  , only: phys_getopts, use_gw_convect_dp

    integer idx

    call phys_getopts(deep_scheme_out=deep_scheme)

    select case (deep_scheme)
    case('ZM') !    Zhang-McFarlane (default)
      call zm_conv_register
    case('off', 'UNICON') ! Off needs to setup the following fields
      call pbuf_add_field('ICWMRDP'    , 'physpkg', dtype_r8, [pcols,pver], icwmrdp_idx)
      call pbuf_add_field('RPRDDP'     , 'physpkg', dtype_r8, [pcols,pver], rprddp_idx)
      call pbuf_add_field('NEVAPR_DPCU', 'physpkg', dtype_r8, [pcols,pver], nevapr_dpcu_idx)
      call pbuf_add_field('PREC_DP'    , 'physpkg', dtype_r8, [pcols]     , prec_dp_idx)
      call pbuf_add_field('SNOW_DP'    , 'physpkg', dtype_r8, [pcols]     , snow_dp_idx)
    end select

    ! If gravity waves from deep convection are on, output this field.
    if (use_gw_convect_dp .and. deep_scheme == 'ZM') then
      call pbuf_add_field('TTEND_DP','physpkg',dtype_r8,[pcols,pver],ttend_dp_idx)
    end if

  end subroutine convect_deep_register

  subroutine convect_deep_init(pref_edge)

    use cam_history   , only: addfld                          
    use pmgrid        , only: plevp
    use spmd_utils    , only: masterproc
    use zm_conv_intr  , only: zm_conv_init
    use cam_abortutils, only: endrun
    use physics_buffer, only: physics_buffer_desc, pbuf_get_index

    real(r8), intent(in) :: pref_edge(plevp) ! Reference pressures at interfaces

    select case (deep_scheme)
    case('off')
      if (masterproc) write(iulog, *) 'convect_deep: no deep convection selected'
    case('CLUBB_SGS')
      if (masterproc) write(iulog, *) 'convect_deep: CLUBB_SGS selected'
    case('ZM')
      if (masterproc) write(iulog, *) 'convect_deep initializing Zhang-McFarlane convection'
      call zm_conv_init(pref_edge)
    case ('UNICON')
      if (masterproc) write(iulog, *)'convect_deep: deep convection done by UNICON'
    case('SPCAM')
      if (masterproc) write(iulog, *)'convect_deep: deep convection done by SPCAM'
      return
    case default
      if (masterproc) write(iulog, *)'WARNING: convect_deep: no deep convection scheme. May fail.'
    end select

    icwmrdp_idx     = pbuf_get_index('ICWMRDP')
    rprddp_idx      = pbuf_get_index('RPRDDP')
    nevapr_dpcu_idx = pbuf_get_index('NEVAPR_DPCU')
    prec_dp_idx     = pbuf_get_index('PREC_DP')
    snow_dp_idx     = pbuf_get_index('SNOW_DP')
    cldtop_idx      = pbuf_get_index('CLDTOP')
    cldbot_idx      = pbuf_get_index('CLDBOT')
    cld_idx         = pbuf_get_index('CLD')
    fracis_idx      = pbuf_get_index('FRACIS')
    pblh_idx        = pbuf_get_index('pblh')
    tpert_idx       = pbuf_get_index('tpert')

    call addfld('ICWMRDP', ['lev'], 'A', 'kg/kg', 'Deep Convection in-cloud water mixing ratio ')

  end subroutine convect_deep_init

  subroutine convect_deep_tend(mcon, cme, pflx, zdu, rliq, rice, ztodt, state, ptend, landfrac, pbuf)

    use physics_types , only: physics_state, physics_ptend, physics_tend, physics_ptend_init
    use cam_history   , only: outfld
    use constituents  , only: pcnst
    use zm_conv_intr  , only: zm_conv_tend
    use cam_history   , only: outfld
    use physconst     , only: cpair
    use physics_buffer, only: physics_buffer_desc, pbuf_get_field

    type(physics_state), intent(in ) :: state   ! Physics state variables
    type(physics_ptend), intent(out) :: ptend   ! Individual parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(in ) :: ztodt              ! 2 delta t (model time increment)
    real(r8), intent(in ) :: landfrac(pcols)    ! Land fraction
    real(r8), intent(out) :: mcon(pcols,pverp)  ! Convective mass flux--m sub c
    real(r8), intent(out) :: pflx(pcols,pverp)  ! Scattered precip flux at each level
    real(r8), intent(out) :: cme(pcols,pver)    ! Cmf condensation - evaporation
    real(r8), intent(out) :: zdu(pcols,pver)    ! Detraining mass flux
    real(r8), intent(out) :: rliq(pcols)        ! Reserved liquid (not yet in cldliq) for energy integrals
    real(r8), intent(out) :: rice(pcols)        ! Reserved ice (not yet in cldice) for energy integrals

    real(r8), pointer :: prec   (:)             ! Total precipitation
    real(r8), pointer :: snow   (:)             ! Snow from ZM convection 
    real(r8), pointer :: jctop  (:)
    real(r8), pointer :: jcbot  (:)
    real(r8), pointer :: cld    (:,:,:)
    real(r8), pointer :: ql     (:,:)           ! wg grid slice of cloud liquid water.
    real(r8), pointer :: rprd   (:,:)           ! rain production rate
    real(r8), pointer :: fracis (:,:,:)         ! fraction of transported species that are insoluble
    real(r8), pointer :: evapcdp(:,:)           ! Evaporation of deep convective precipitation
    real(r8), pointer :: pblh   (:)             ! Planetary boundary layer height
    real(r8), pointer :: tpert  (:)             ! Thermal temperature excess
    real(r8), pointer :: ttend_dp(:,:)          ! Temperature tendency from deep convection (pbuf pointer).
    real(r8) zero(pcols,pver)
    integer i, k

    call pbuf_get_field(pbuf, cldtop_idx , jctop)
    call pbuf_get_field(pbuf, cldbot_idx , jcbot)
    call pbuf_get_field(pbuf, icwmrdp_idx, ql   )

    select case (deep_scheme)
    case('off', 'UNICON', 'CLUBB_SGS') ! in UNICON case the run method is called from convect_shallow_tend
      zero = 0     
      mcon = 0
      pflx = 0
      cme  = 0
      zdu  = 0
      rliq = 0
      rice = 0

      call physics_ptend_init(ptend, state%psetcols, 'convect_deep')

      call pbuf_get_field(pbuf, cld_idx        , cld    , start=[1,1]  , kount=[pcols,pver]) 
      call pbuf_get_field(pbuf, rprddp_idx     , rprd   )
      call pbuf_get_field(pbuf, fracis_idx     , fracis , start=[1,1,1], kount=[pcols,pver,pcnst])
      call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp)
      call pbuf_get_field(pbuf, prec_dp_idx    , prec   )
      call pbuf_get_field(pbuf, snow_dp_idx    , snow   )

      prec    = 0
      snow    = 0
      jctop   = pver
      jcbot   = 1._r8
      cld     = 0
      ql      = 0
      rprd    = 0
      fracis  = 0
      evapcdp = 0
    case ('ZM')
      call pbuf_get_field(pbuf, pblh_idx,  pblh)
      call pbuf_get_field(pbuf, tpert_idx, tpert)
      call zm_conv_tend(pblh, mcon, cme, tpert, pflx, zdu, rliq, rice, ztodt, jctop, jcbot, &
                        state, ptend, landfrac, pbuf)
    end select

    ! If we added temperature tendency to pbuf, set it now.
    if (ttend_dp_idx > 0) then
      call pbuf_get_field(pbuf, ttend_dp_idx, ttend_dp)
      ttend_dp(:state%ncol,:pver) = ptend%s(:state%ncol,:pver) / cpair
    end if

    call outfld('ICWMRDP ', ql, pcols, state%lchnk)

  end subroutine convect_deep_tend

  subroutine convect_deep_tend_2(state,  ptend,  ztodt, pbuf)

    use physics_types , only: physics_state, physics_ptend, physics_ptend_init
    use physics_buffer, only: physics_buffer_desc
    use constituents  , only: pcnst
    use zm_conv_intr  , only: zm_conv_tend_2

    type(physics_state), intent(in ) :: state          ! Physics state variables
    type(physics_ptend), intent(out) :: ptend          ! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(in) :: ztodt                      ! 2 delta t (model time increment)

    if (deep_scheme == 'ZM') then
      call zm_conv_tend_2(state, ptend, ztodt, pbuf) 
    else
      call physics_ptend_init(ptend, state%psetcols, 'convect_deep')
    end if

  end subroutine convect_deep_tend_2

end module convect_deep
