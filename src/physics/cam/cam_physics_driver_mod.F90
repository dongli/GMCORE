module cam_physics_driver_mod

  use shr_kind_mod    , only: r8 => SHR_KIND_R8
  use shr_flux_mod    , only: shr_flux_atmOcn
  use shr_orb_mod     , only: SHR_ORB_UNDEF_INT, SHR_ORB_UNDEF_REAL, shr_orb_params
  use spmd_utils      , only: spmdinit
  use check_energy    , only: check_energy_timestep_init
  use dyn_grid        , only: dyn_grid_init, dyn_grid_final
  use dyn_comp        , only: dyn_init, dyn_final, dyn_import_t, dyn_export_t
  use phys_grid       , only: phys_grid_init, local_dp_map, get_lon_all_p, get_lat_all_p
  use physics_types   , only: physics_state, physics_tend, set_state_pdry
  use physics_buffer  , only: physics_buffer_desc, pbuf_get_chunk
  use physpkg         , only: phys_register, phys_init, phys_run1, phys_run2, phys_final
  use physconst       , only: cappa, cpair, gravit
  use ppgrid          , only: begchunk, endchunk, pcols, pver, pverp
  use qneg_module     , only: qneg3
  use constituents    , only: pcnst, cnst_name, cnst_longname, cnst_get_type_byind, qmin
  use chem_surfvals   , only: chem_surfvals_init
  use air_composition , only: air_composition_init
  use camsrfexch      , only: cam_out_t, cam_in_t, hub2atm_alloc, atm2hub_alloc, cam_export
  use runtime_opts    , only: read_namelist
  use shr_cal_mod     , only: shr_cal_gregorian
  use shr_pio_mod     , only: shr_pio_init1, shr_pio_init2
  use cam_pio_utils   , only: init_pio_subsystem
  use cam_instance    , only: cam_instance_init
  use cam_initfiles   , only: cam_initfiles_open
  use cam_control_mod , only: cam_ctrl_init, cam_ctrl_set_orbit
  use cam_history     , only: intht
  use orbit           , only: zenith
  use time_manager    , only: timemgr_init, get_curr_calday, get_curr_date, advance_timestep

  use flogger
  use string            , only: to_int
  use formula_mod       , only: wet_mixing_ratio
  use namelist_mod      , only: dt_phys, dt_adv, restart, cam_namelist_path, case_name, case_desc, use_aqua_planet, filter_ptend
  use process_mod       , only: proc, process_barrier
  use time_mod          , only: start_time, end_time, curr_time
  use tracer_mod        , only: tracer_add, tracers
  use albedo_mod        , only: albedo_ocnice
  use physics_types_mod , only: physics_use_wet_tracers
  use cam_physics_types_mod
  use cam_physics_objects_mod
  use aquaplanet_test_mod

  implicit none

  private

  public cam_physics_init
  public cam_physics_final
  public cam_physics_run1
  public cam_physics_sfc_flux
  public cam_physics_run2
  public cam_physics_d2p
  public cam_physics_p2d
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

contains

  subroutine cam_physics_init(namelist_path, mesh, dt_adv, dt_phys)

    character(*), intent(in) :: namelist_path
    type(physics_mesh_type), intent(in) :: mesh(:)
    real(r8), intent(in) :: dt_adv
    real(r8), intent(in) :: dt_phys

    integer ncol, c, m
    real(r8) :: eccen = SHR_ORB_UNDEF_REAL
    real(r8) :: obliq = SHR_ORB_UNDEF_REAL
    real(r8) :: mvelp = SHR_ORB_UNDEF_REAL
    real(r8) :: obliqr = SHR_ORB_UNDEF_REAL
    real(r8) :: lambm0 = SHR_ORB_UNDEF_REAL
    real(r8) :: mvelpp = SHR_ORB_UNDEF_REAL

    call shr_orb_params(iyear_AD=2000, eccen=eccen, obliq=obliq, mvelp=mvelp, obliqr=obliqr, lambm0=lambm0, mvelpp=mvelpp, log_print=.true.)

    comp_comm(1) = proc%comm_model
    comp_comm_iam(1) = proc%id_model
    comp_iamin(1) = .true.

    call spmdinit(proc%comm_model)
    call cam_instance_init(atm_id, 'atm', 1, '_1')
    call shr_pio_init1(ncomps, trim(merge(namelist_path, cam_namelist_path, cam_namelist_path == 'N/A')), proc%comm_model)
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
    call timemgr_init( &
      dtime_in=int(dt_phys), &
      calendar_in=shr_cal_gregorian, &
      start_ymd=to_int(start_time%format('%Y%m%d')), &
      start_tod=start_time%hour*3600+start_time%minute*60+start_time%second, &
      ref_ymd=20000101, &
      ref_tod=0, &
      stop_ymd=to_int(end_time%format('%Y%m%d')), &
      stop_tod=end_time%hour*3600+end_time%minute*60+end_time%second, &
      curr_ymd=to_int(curr_time%format('%Y%m%d')), &
      curr_tod=curr_time%hour*3600+curr_time%minute*60+curr_time%second, &
      perpetual_run=.false., &
      perpetual_ymd=20000101, &
      initial_run=.not. restart &
    )
    call read_namelist(trim(merge(namelist_path, cam_namelist_path, cam_namelist_path == 'N/A')))
    call cam_initfiles_open()
    call dyn_grid_init()
    call phys_grid_init()
    call phys_register()
    call dyn_init(dyn_in, dyn_out)
    call chem_surfvals_init()
    call air_composition_init()

    if (allocated(physics_use_wet_tracers)) deallocate(physics_use_wet_tracers)
    allocate(physics_use_wet_tracers(pcnst))
    do m = 1, pcnst
      call tracer_add('cam_cnst', dt_adv, cnst_name(m), cnst_longname(m), 'kg kg-1', type=0)
      physics_use_wet_tracers(m) = cnst_get_type_byind(m) == 'wet'
    end do

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

    ! call intht('gmcore')

  end subroutine cam_physics_init

  subroutine cam_physics_run1()

    integer ncol, c
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
    call phys_run1(dt_phys, phys_state, phys_tend, pbuf2d, cam_in, cam_out)

  end subroutine cam_physics_run1

  subroutine cam_physics_sfc_flux()

    integer ncol, c
    real(r8) zeros(pcols)

    zeros = 0
    do c = begchunk, endchunk
      ncol = cam_in(c)%ncol
      call shr_flux_atmOcn(             &
        nMax=ncol                     , &
        mask=int(cam_in(c)%ocnfrac)   , &
        ocn_surface_flux_scheme=0     , &
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
        us=zeros                      , &
        vs=zeros                      , &
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
        duu10n=cam_in(c)%u10          )
      cam_in(c)%u10 (:ncol  ) = sqrt(cam_in(c)%u10(:ncol))
      cam_in(c)%lwup(:ncol  ) = -cam_in(c)%lwup(:ncol)
      cam_in(c)%shf (:ncol  ) = -cam_in(c)%shf(:ncol)
      cam_in(c)%lhf (:ncol  ) = -cam_in(c)%lhf(:ncol)
      cam_in(c)%cflx(:ncol,1) = -cam_in(c)%cflx(:ncol,1)
    end do

  end subroutine cam_physics_sfc_flux

  subroutine cam_physics_run2()

    call phys_run2(phys_state, dt_phys, phys_tend, pbuf2d,  cam_out, cam_in)
    call advance_timestep()

  end subroutine cam_physics_run2

  subroutine cam_physics_final()

    call phys_final(phys_state, phys_tend , pbuf2d)
    call dyn_grid_final()
    call dyn_final()

  end subroutine cam_physics_final

  subroutine cam_physics_d2p()

    integer ncol, c, i, k, m
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    if (local_dp_map) then
      ! TODO: Check if begchunk is 1 and endchunk is mesh%ncol.
      do c = begchunk, endchunk
        associate (mesh  => objects(c)%mesh , &
                   state => objects(c)%state)
        do i = 1, mesh%ncol
          phys_state(c)%ps(i) = state%ps(i)
          phys_state(c)%phis(i) = state%zs(i) * gravit
        end do
        do k = 1, mesh%nlev
          do i = 1, mesh%ncol
            phys_state(c)%u      (i,k) = state%u(i,k)
            phys_state(c)%v      (i,k) = state%v(i,k)
            phys_state(c)%t      (i,k) = state%t(i,k)
            phys_state(c)%exner  (i,k) = (state%ps(i) / state%p(i,k))**cappa
            phys_state(c)%omega  (i,k) = state%omg(i,k)
            phys_state(c)%pmid   (i,k) = state%p(i,k)
            phys_state(c)%pdel   (i,k) = state%dp(i,k)
            phys_state(c)%pdeldry(i,k) = state%dp_dry(i,k)
            phys_state(c)%rpdel  (i,k) = 1.0_r8 / phys_state(c)%pdel(i,k)
            phys_state(c)%lnpmid (i,k) = log(phys_state(c)%pmid(i,k))
            phys_state(c)%zm     (i,k) = state%z(i,k)
          end do
        end do
        do k = 1, mesh%nlev + 1
          do i = 1, mesh%ncol
            phys_state(c)%pint   (i,k) = state%p_lev(i,k)
            phys_state(c)%lnpint (i,k) = log(phys_state(c)%pint(i,k))
            phys_state(c)%zi     (i,k) = state%z_lev(i,k)
          end do
        end do
        do m = 1, pcnst
          do k = 1, mesh%nlev
            do i = 1, mesh%ncol
              phys_state(c)%q(i,k,m) = state%q(i,k,m)
            end do
          end do
        end do
        end associate
      end do
    else
      call log_error('cam_physics_d2p: Distributed physics columns are not supported yet!', __FILE__, __LINE__)
    end if

    do c = begchunk, endchunk
      ncol = phys_state(c)%ncol
      ! Compute initial dry static energy, include surface geopotential
      do k = 1, pver
        do i = 1, ncol
          phys_state(c)%s(i,k) = cpair * phys_state(c)%t(i,k) &
                               + gravit * phys_state(c)%zm(i,k) + phys_state(c)%phis(i)
        end do
      end do

      call set_state_pdry(phys_state(c), pdeld_calc=.false.)

      call qneg3('cam_physics_d2p', c, ncol, pcols, pver, 1, pcnst, qmin, phys_state(c)%q)

      ! Compute energy and water integrals of input state
      pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
      call check_energy_timestep_init(phys_state(c), phys_tend(c), pbuf_chnk)
    end do

  end subroutine cam_physics_d2p

  subroutine cam_physics_p2d()

    integer ncol, c, i, j, k, m
    integer ilon(pcols), jlat(pcols)

    if (local_dp_map) then
      do c = begchunk, endchunk
        associate (mesh  => objects(c)%mesh , &
                   state => objects(c)%state, &
                   tend  => objects(c)%tend , &
                   q     => tracers(c)%q%d  , &
                   qm    => tracers(c)%qm%d )
        call get_lon_all_p(c, mesh%ncol, ilon)
        call get_lat_all_p(c, mesh%ncol, jlat)
        do k = 1, mesh%nlev
          do i = 1, mesh%ncol
            tend%dudt(i,k) = phys_tend(c)%dudt(i,k)
            tend%dvdt(i,k) = phys_tend(c)%dvdt(i,k)
            tend%dtdt(i,k) = phys_tend(c)%dtdt(i,k)
          end do
        end do
        tend%updated_u = .true.
        tend%updated_v = .true.
        tend%updated_t = .true.
        do m = 1, pcnst
          if (physics_use_wet_tracers(m)) then
            do k = 1, pver
              do i = 1, ncol
                tend%dqdt(i,k,m) = (phys_state(c)%q(i,k,m) - wet_mixing_ratio(q(ilon(i),jlat(i),k,m), qm(ilon(i),jlat(i),k))) / dt_phys
              end do
            end do
          else
            do k = 1, pver
              do i = 1, ncol
                tend%dqdt(i,k,m) = (phys_state(c)%q(i,k,m) - q(ilon(i),jlat(i),k,m)) / dt_phys
              end do
            end do
          end if
          tend%updated_q(m) = .true.
        end do
        end associate
      end do
    else
      call log_error('cam_physics_p2d: Distributed physics columns are not supported yet!', __FILE__, __LINE__)
    end if

  end subroutine cam_physics_p2d

end module cam_physics_driver_mod
