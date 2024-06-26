! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module cam_physics_driver_mod

  ! CAM modules
  use shr_kind_mod      , only: r8 => SHR_KIND_R8
  use shr_flux_mod      , only: flux_atmOcn
  use shr_log_mod       , only: shr_log_Unit
  use shr_orb_mod       , only: SHR_ORB_UNDEF_INT, SHR_ORB_UNDEF_REAL, shr_orb_params
  use spmd_utils        , only: spmdinit
  use check_energy      , only: check_energy_timestep_init
  use dyn_grid          , only: dyn_grid_init, dyn_grid_final
  use dyn_comp          , only: dyn_init, dyn_final, dyn_import_t, dyn_export_t
  use phys_grid         , only: phys_grid_init, local_dp_map
  use physics_types     , only: physics_state, physics_tend, set_state_pdry
  use physics_buffer    , only: physics_buffer_desc, pbuf_get_chunk
  use physpkg           , only: phys_register, phys_init, phys_run1, phys_run2, phys_final
  use physconst         , only: cappa, cpair, gravit
  use ppgrid            , only: begchunk, endchunk, pcols, pver, pverp
  use qneg_module       , only: qneg3
  use constituents      , only: pcnst, cnst_name, cnst_longname, cnst_get_type_byind, qmin
  use chem_surfvals     , only: chem_surfvals_init
  use air_composition   , only: air_composition_init
  use camsrfexch        , only: cam_out_t, cam_in_t, hub2atm_alloc, atm2hub_alloc, cam_export
  use runtime_opts      , only: read_namelist
  use shr_cal_mod       , only: shr_cal_gregorian
  use shr_pio_mod       , only: shr_pio_init1, shr_pio_init2
  use cam_pio_utils     , only: init_pio_subsystem
  use cam_instance      , only: cam_instance_init
  use cam_initfiles     , only: cam_initfiles_open
  use cam_control_mod   , only: cam_ctrl_init, cam_ctrl_set_orbit
  use cam_history       , only: intht
  use orbit             , only: zenith
  use time_manager      , only: timemgr_init, get_curr_calday, get_curr_date, advance_timestep
  ! GMCORE modules
  use flogger
  use string            , only: to_int
  use formula_mod       , only: dry_mixing_ratio
  use namelist_mod      , only: dt_phys, dt_adv, restart, cam_namelist_path, case_name, case_desc, use_aqua_planet, filter_ptend
  use process_mod       , only: proc, process_barrier
  use time_mod          , only: start_time, end_time, curr_time
  use tracer_mod        , only: tracer_add, tracers
  use albedo_mod        , only: albedo_ocnice
  use physics_types_mod , only: physics_use_wet_tracers
  use cam_physics_types_mod
  use cam_physics_objects_mod
  use cam_physics_output_mod
  use aquaplanet_test_mod

  implicit none

  private

  public cam_physics_init_stage2
  public cam_physics_final
  public cam_physics_run_stage1
  public cam_physics_sfc_flux
  public cam_physics_run_stage2
  public cam_physics_d2p
  public cam_physics_p2d
  public cam_physics_add_output
  public cam_physics_output
  public objects

  type(physics_state), pointer :: phys_state(:) => null()
  type(physics_tend), pointer :: phys_tend(:) => null()
  type(physics_buffer_desc), pointer :: pbuf2d(:,:) => null()

  type(dyn_import_t) dyn_in
  type(dyn_export_t) dyn_out
  type(cam_in_t), pointer :: cam_in(:) => null()
  type(cam_out_t), pointer :: cam_out(:) => null()

  integer, parameter :: atm_id = 1
  integer, parameter :: ncomps = 1
  integer, parameter :: comp_id(ncomps) = [atm_id]
  character(3), parameter :: comp_name(ncomps) = ['atm']
  logical comp_iamin(ncomps)
  integer comp_comm(ncomps)
  integer comp_comm_iam(ncomps)

  ! Local map between CAM physics chunked columns and GMCORE dynamics columns.
  integer, allocatable :: icol_map(:,:)

contains

  subroutine cam_physics_init_stage2(namelist_path, mesh, dt_adv, dt_phys)

    character(*), intent(in) :: namelist_path
    type(physics_mesh_type), intent(in) :: mesh(:)
    real(r8), intent(in) :: dt_adv
    real(r8), intent(in) :: dt_phys

    integer ncol, c, i, icol, m
    real(r8) :: eccen  = SHR_ORB_UNDEF_REAL
    real(r8) :: obliq  = SHR_ORB_UNDEF_REAL
    real(r8) :: mvelp  = SHR_ORB_UNDEF_REAL
    real(r8) :: obliqr = SHR_ORB_UNDEF_REAL
    real(r8) :: lambm0 = SHR_ORB_UNDEF_REAL
    real(r8) :: mvelpp = SHR_ORB_UNDEF_REAL

    call shr_orb_params( &
      iyear_AD=2000    , &
      eccen=eccen      , &
      obliq=obliq      , &
      mvelp=mvelp      , &
      obliqr=obliqr    , &
      lambm0=lambm0    , &
      mvelpp=mvelpp    , &
      log_print=.true. )

    comp_comm(1) = proc%comm_model
    comp_comm_iam(1) = proc%id_model
    comp_iamin(1) = .true.

    call spmdinit(proc%comm_model)
    call cam_instance_init(atm_id, 'atm', 1, '_1')
    call shr_pio_init1(ncomps, merge(namelist_path, cam_namelist_path, cam_namelist_path=='N/A'), proc%comm_model)
    call shr_pio_init2(comp_id, comp_name, comp_iamin, comp_comm, comp_comm_iam)
    call init_pio_subsystem()
    call cam_ctrl_init(                &
      caseid_in=case_name            , &
      ctitle_in=case_desc            , &
      initial_run_in=.not.restart    , &
      restart_run_in=restart         , &
      branch_run_in=.false.          , &
      post_assim_in=.false.          , &
      aqua_planet_in=use_aqua_planet , &
      brnch_retain_casename_in=.true.)
    call cam_ctrl_set_orbit(eccen_in=eccen, obliqr_in=obliqr, lambm0_in=lambm0, mvelpp_in=mvelpp)
    call timemgr_init(                                                        &
      dtime_in=int(dt_phys)                                                 , &
      calendar_in=shr_cal_gregorian                                         , &
      start_ymd=to_int(start_time%format('%Y%m%d'))                         , &
      start_tod=start_time%hour*3600+start_time%minute*60+start_time%second , &
      ref_ymd=20000101                                                      , &
      ref_tod=0                                                             , &
      stop_ymd=to_int(end_time%format('%Y%m%d'))                            , &
      stop_tod=end_time%hour*3600+end_time%minute*60+end_time%second        , &
      curr_ymd=to_int(curr_time%format('%Y%m%d'))                           , &
      curr_tod=curr_time%hour*3600+curr_time%minute*60+curr_time%second     , &
      perpetual_run=.false.                                                 , &
      perpetual_ymd=20000101                                                , &
      initial_run=.not. restart                                             )
    call read_namelist(merge(namelist_path, cam_namelist_path, cam_namelist_path=='N/A'))
    call cam_initfiles_open()
    call dyn_grid_init()
    call phys_grid_init()
    call phys_register()

    if (allocated(physics_use_wet_tracers)) deallocate(physics_use_wet_tracers)
    allocate(physics_use_wet_tracers(pcnst))
    do m = 1, pcnst
      call tracer_add('cam_cnst', dt_adv, cnst_name(m), cnst_longname(m), 'kg kg-1', type=0)
      physics_use_wet_tracers(m) = cnst_get_type_byind(m) == 'wet'
    end do

    call cam_physics_objects_init(mesh)
    call dyn_init(dyn_in, dyn_out)
    call chem_surfvals_init()
    call air_composition_init()

    if (.not. restart) then
      call hub2atm_alloc(cam_in)
      call atm2hub_alloc(cam_out)
    end if

    call phys_init(phys_state, phys_tend, pbuf2d, cam_in, cam_out)

    if (use_aqua_planet) then
      do c = begchunk, endchunk
        ncol = phys_state(c)%ncol
        call aquaplanet_test_set_bc( &
          phys_state(c)%lon (:ncol), &
          phys_state(c)%lat (:ncol), &
          cam_in(c)%sst     (:ncol), &
          cam_in(c)%landfrac(:ncol), &
          cam_in(c)%ocnfrac (:ncol), &
          cam_in(c)%icefrac (:ncol))
      end do
    end if

    ! Create column map.
    allocate(icol_map(pcols,begchunk:endchunk))
    icol = 1
    do c = begchunk, endchunk
      do i = 1, pcols
        icol_map(i,c) = icol
        icol = icol + 1
      end do
    end do

  end subroutine cam_physics_init_stage2

  subroutine cam_physics_run_stage1()

    integer ncol, c, i
    real(r8) coszrs(pcols)
    real(r8) calday
    logical, save :: first_call = .true.
    type(physics_buffer_desc), pointer :: pbuf(:)

    if (first_call) then
      do c = begchunk, endchunk
        pbuf => pbuf2d(:,c)
        call cam_export(phys_state(c), cam_out(c), pbuf)
      end do
      call cam_physics_sfc_flux()
      first_call = .false.
    end if

    calday = get_curr_calday()
    do c = begchunk, endchunk
      ncol = phys_state(c)%ncol
      call zenith(calday, phys_state(c)%lat, phys_state(c)%lon, coszrs, ncol)
      call albedo_ocnice(           &
        ncol=ncol                 , &
        lndfrac=cam_in(c)%landfrac, &
        ocnfrac=cam_in(c)%ocnfrac , &
        icefrac=cam_in(c)%icefrac , &
        coszrs=coszrs             , &
        asdir=cam_in(c)%asdir     , &
        asdif=cam_in(c)%asdif     , &
        aldir=cam_in(c)%aldir     , &
        aldif=cam_in(c)%aldif     )
    end do

    ! Copy into GMCORE physics state.
    associate (pstate => objects(1)%state)
    do c = begchunk, endchunk
      ncol = phys_state(c)%ncol
      do i = 1, ncol
        pstate%asdir(icol_map(i,c)) = cam_in(c)%asdir(i)
        pstate%asdif(icol_map(i,c)) = cam_in(c)%asdif(i)
        pstate%aldir(icol_map(i,c)) = cam_in(c)%aldir(i)
        pstate%aldif(icol_map(i,c)) = cam_in(c)%aldif(i)
      end do
    end do
    end associate
    call phys_run1(dt_phys, phys_state, phys_tend, pbuf2d, cam_in, cam_out)

  end subroutine cam_physics_run_stage1

  subroutine cam_physics_sfc_flux()

    integer ncol, c, i
    real(r8) zeros(pcols)

    zeros = 0
    do c = begchunk, endchunk
      ncol = cam_in(c)%ncol
      call flux_atmOcn(                 &
        logunit=shr_log_Unit          , &
        nMax=ncol                     , &
        mask=int(cam_in(c)%ocnfrac)   , &
        ocn_surface_flux_scheme=0     , &
        add_gusts=.true.              , &
        zbot=cam_out(c)%zbot          , &
        ubot=cam_out(c)%ubot          , &
        vbot=cam_out(c)%vbot          , &
        thbot=cam_out(c)%thbot        , &
        qbot=cam_out(c)%qbot(:ncol,1) , &
        s16O=zeros                    , &
        sHDO=zeros                    , &
        s18O=zeros                    , &
        r16O=zeros                    , &
        rHDO=zeros                    , &
        r18O=zeros                    , &
        rbot=cam_out(c)%rho           , &
        tbot=cam_out(c)%tbot          , &
        rainc=cam_out(c)%precc        , &
        us=zeros                      , &
        vs=zeros                      , &
        pslv=cam_out(c)%psl           , &
        ts=cam_in(c)%sst              , &
        seq_flux_atmocn_minwind=0.5_r8, &
        sen=cam_in(c)%shf             , &
        lat=cam_in(c)%lhf             , &
        lwup=cam_in(c)%lwup           , &
        evap=cam_in(c)%cflx(:ncol,1)  , &
        evap_16O=zeros                , &
        evap_HDO=zeros                , &
        evap_18O=zeros                , &
        taux=cam_in(c)%wsx            , &
        tauy=cam_in(c)%wsy            , &
        tref=cam_in(c)%tref           , &
        qref=cam_in(c)%qref           , &
        duu10n=cam_in(c)%u10          , &
        ugust_out=cam_in(c)%ugust_out , &
        u10res=cam_in(c)%u10_with_gust)
      cam_in(c)%u10 (:ncol  ) = sqrt(cam_in(c)%u10(:ncol))
      cam_in(c)%wsx (:ncol  ) = -cam_in(c)%wsx (:ncol  )
      cam_in(c)%wsy (:ncol  ) = -cam_in(c)%wsy (:ncol  )
      cam_in(c)%lwup(:ncol  ) = -cam_in(c)%lwup(:ncol  )
      cam_in(c)%shf (:ncol  ) = -cam_in(c)%shf (:ncol  )
      cam_in(c)%lhf (:ncol  ) = -cam_in(c)%lhf (:ncol  )
      cam_in(c)%cflx(:ncol,1) = -cam_in(c)%cflx(:ncol,1)
    end do

    ! Copy into GMCORE physics state.
    associate (pstate => objects(1)%state)
    do c = begchunk, endchunk
      do i = 1, cam_in(c)%ncol
        pstate%wsp10(icol_map(i,c)  ) = cam_in(c)%u10 (i)
        pstate%qflx (icol_map(i,c),1) = cam_in(c)%cflx(i,1)
      end do
    end do
    end associate

  end subroutine cam_physics_sfc_flux

  subroutine cam_physics_run_stage2()

    call phys_run2(phys_state, dt_phys, phys_tend, pbuf2d,  cam_out, cam_in)
    call advance_timestep()

  end subroutine cam_physics_run_stage2

  subroutine cam_physics_final()

    call phys_final(phys_state, phys_tend , pbuf2d)
    call dyn_grid_final()
    call dyn_final()
    call cam_physics_objects_final()

    if (allocated(icol_map)) deallocate(icol_map)

  end subroutine cam_physics_final

  subroutine cam_physics_d2p()

    integer c, i, k, m
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    associate (mesh => objects(1)%mesh, pstate => objects(1)%state)

    if (local_dp_map) then
      ! Copy data into CAM physics state object.
      ! Assume only one block in dycore.
      do c = begchunk, endchunk
        do i = 1, phys_state(c)%ncol
          phys_state(c)%ps  (i) = pstate%ps(icol_map(i,c))
          phys_state(c)%phis(i) = pstate%zs(icol_map(i,c)) * gravit
        end do
        do k = 1, mesh%nlev
          do i = 1, phys_state(c)%ncol
            phys_state(c)%u      (i,k) = pstate %u     (icol_map(i,c),k)
            phys_state(c)%v      (i,k) = pstate %v     (icol_map(i,c),k)
            phys_state(c)%t      (i,k) = pstate %t     (icol_map(i,c),k)
            phys_state(c)%exner  (i,k) = (pstate%ps    (icol_map(i,c)) / pstate%p(icol_map(i,c),k))**cappa
            phys_state(c)%omega  (i,k) = pstate %omg   (icol_map(i,c),k)
            phys_state(c)%pmid   (i,k) = pstate %p     (icol_map(i,c),k)
            phys_state(c)%pdel   (i,k) = pstate %dp    (icol_map(i,c),k)
            phys_state(c)%pdeldry(i,k) = pstate %dp_dry(icol_map(i,c),k)
            phys_state(c)%rpdel  (i,k) = 1.0_r8 / phys_state(c)%pdel(i,k)
            phys_state(c)%lnpmid (i,k) = log(phys_state(c)%pmid(i,k))
            phys_state(c)%zm     (i,k) = pstate %z     (icol_map(i,c),k)
          end do
        end do
        do k = 1, mesh%nlev + 1
          do i = 1, phys_state(c)%ncol
            phys_state(c)%pint   (i,k) = pstate%p_lev(icol_map(i,c),k)
            phys_state(c)%zi     (i,k) = pstate%z_lev(icol_map(i,c),k)
            phys_state(c)%lnpint (i,k) = log(phys_state(c)%pint(i,k))
          end do
        end do
        do m = 1, pcnst
          do k = 1, mesh%nlev
            do i = 1, phys_state(c)%ncol
              phys_state(c)%q(i,k,m) = pstate%q(icol_map(i,c),k,m)
            end do
          end do
        end do
      end do
    else
      call log_error('cam_physics_d2p: Distributed physics columns are not supported yet!', __FILE__, __LINE__)
    end if

    do c = begchunk, endchunk
      ! Compute initial dry static energy, include surface geopotential
      do k = 1, mesh%nlev
        do i = 1, phys_state(c)%ncol
          phys_state(c)%s(i,k) = cpair  * phys_state(c)%t   (i,k) &
                               + gravit * phys_state(c)%zm  (i,k) &
                               +          phys_state(c)%phis(i)
        end do
      end do

      call set_state_pdry(phys_state(c), pdeld_calc=.false.)

      call qneg3('cam_physics_d2p', c, phys_state(c)%ncol, pcols, mesh%nlev, 1, pcnst, qmin, phys_state(c)%q)

      ! Compute energy and water integrals of input state
      pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
      call check_energy_timestep_init(phys_state(c), phys_tend(c), pbuf_chnk)
    end do

    end associate

  end subroutine cam_physics_d2p

  subroutine cam_physics_p2d()

    integer c, i, j, k, m

    associate (mesh => objects(1)%mesh, pstate => objects(1)%state, ptend => objects(1)%tend)

    if (local_dp_map) then
      do c = begchunk, endchunk
        do k = 1, mesh%nlev
          do i = 1, phys_state(c)%ncol
            ptend%dudt(icol_map(i,c),k) = phys_tend(c)%dudt(i,k)
            ptend%dvdt(icol_map(i,c),k) = phys_tend(c)%dvdt(i,k)
            ptend%dtdt(icol_map(i,c),k) = phys_tend(c)%dtdt(i,k)
          end do
        end do
        ptend%updated_u = .true.
        ptend%updated_v = .true.
        ptend%updated_t = .true.
        do m = 1, pcnst
          do k = 1, mesh%nlev
            do i = 1, phys_state(c)%ncol
              ptend%dqdt(icol_map(i,c),k,m) = (phys_state(c)%q(i,k,m) - pstate%q(icol_map(i,c),k,m)) / dt_phys
              pstate%q(icol_map(i,c),k,m) = phys_state(c)%q(i,k,m)
            end do
          end do
          ptend%updated_q(m) = .true.
        end do
      end do
    else
      call log_error('cam_physics_p2d: Distributed physics columns are not supported yet!', __FILE__, __LINE__)
    end if

    end associate

  end subroutine cam_physics_p2d

end module cam_physics_driver_mod
