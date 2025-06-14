! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module namelist_mod

  use string
  use flogger
  use const_mod, only: r8, const_init, time_scale, earth_day_seconds, mars_sol_seconds

  implicit none

  character(30)   :: planet               = 'earth'
  logical         :: use_aqua_planet      = .false.

  integer         :: start_time(5)        = 0
  integer         :: end_time(5)          = 0
  real(r8)        :: run_hours            = 0
  real(r8)        :: run_days             = 0
  real(r8)        :: run_years            = 0
  real(r8)        :: run_my               = 0
  real(r8)        :: run_sol              = 0

  ! Individual time step sizes in target planet time system
  real(r8)        :: dt_dyn               = 0
  real(r8)        :: dt_adv               = 0
  real(r8)        :: dt_phys              = 0
  ! Physics-dynamics coupling type
  ! 1 : Update dynamics and tracers after dynamics.
  ! 2 : Update dynamics and tracers after advection.
  ! 3 : Update dynamics and tracers after physics.
  ! 14: Update modified potential temperature in RK substeps and others like 1.
  integer         :: pdc_type             = 1

  character(256)  :: case_desc            = 'N/A'
  character(256)  :: case_name            = 'N/A'
  character(30 )  :: test_case            = 'N/A'
  character(30)   :: initial_time         = '2000-01-01T00:00:00'
  character(30 )  :: history_interval(1)  = 'N/A'
  character(30 )  :: restart_interval     = 'N/A'
  character(30 )  :: print_interval       = '1 hours'
  character(256)  :: initial_file         = 'N/A'
  character(256)  :: restart_file         = 'N/A'
  character(256)  :: topo_file            = 'N/A'
  character(30 )  :: topo_type            = 'etopo1' ! etopo1, gmted, mola32
  character(256)  :: bkg_file             = 'N/A'
  character(30 )  :: bkg_type             = 'era5'

  integer         :: nlon
  integer         :: nlat
  integer         :: nlev                 = 1

  logical         :: baroclinic           = .false.
  logical         :: hydrostatic          = .true.
  logical         :: nonhydrostatic       = .false.
  logical         :: ideal_dry_core       = .false.
  logical         :: advection            = .false.
  logical         :: restart              = .false.

  logical         :: prepare_regrid_gz    = .true.
  logical         :: init_hydrostatic_gz  = .false.

  character(30)   :: physics_suite        = 'N/A'
  character(30)   :: mp_scheme            = 'N/A'
  character(30)   :: pbl_scheme           = 'N/A'
  character(256)  :: cam_namelist_path    = 'N/A'
  logical         :: filter_ptend         = .false.

  character(256)  :: gmcore_data_dir      = 'N/A'

  integer         :: nproc_io             = 0         ! Number of processes for IO (asynchronized output)
  integer         :: nproc_x(20)          = 0
  integer         :: nproc_y(20)          = 0
  integer         :: lon_hw               = 3
  integer         :: lat_hw               = 3
  character(30)   :: proc_layout          = 'lon>lat' ! or 'lat>lon'
  logical         :: use_async_io         = .false.   ! Use asynchronized output
  integer         :: proc_io_stride       = 0         ! ID stride for IO processes

  character(30)   :: tangent_wgt_scheme   = 'classic'

  real(r8)        :: implicit_w_wgt       = 0.55

  character(30)   :: vert_coord_scheme    = 'hybrid'
  character(30)   :: vert_coord_template  = 'N/A'
  character(30)   :: refer_state_scheme   = 'wrf'
  real(r8)        :: ptop                 = 2.194e2_r8
  real(r8)        :: hybrid_coord_p0      = 1.0e5_r8

  ! Parameters for generating hybrid levels from WRF.
  real(r8)        :: tiso                 = 300.0_r8  ! Isothermal temperature (K)
  real(r8)        :: dzbot                = 10.0_r8   ! Bottom layer thickness (m)
  real(r8)        :: dzmax                = 4000.0_r8 ! Maximum layer thickness (m)
  real(r8)        :: dzstretch_s          = 1.3_r8    ! Stretching factor for surface layers
  real(r8)        :: dzstretch_u          = 1.1_r8    ! Stretching factor for upper layers
  real(r8)        :: eta_b                = 0.2_r8    ! Transition eta for hybrid levels

  ! Parameters for generating hybrid levels from NCEP.
  real(r8)        :: hybrid_coord_ncep_psig   = 0
  real(r8)        :: hybrid_coord_ncep_ppre   = 0
  real(r8)        :: hybrid_coord_ncep_dpbot  = 0
  real(r8)        :: hybrid_coord_ncep_dpsig  = 0
  real(r8)        :: hybrid_coord_ncep_dppre  = 0
  real(r8)        :: hybrid_coord_ncep_dptop  = 0

  integer         :: ke_scheme            = 2
  real(r8)        :: ke_cell_wgt          = 0.5_r8

  character(30)   :: pv_adv_scheme        = 'weno'   ! midpoint, upwind, weno
  logical         :: pv_pole_stokes       = .true.
  integer         :: upwind_order_pv      = 5
  real(r8)        :: upwind_wgt_pv        = 1
  integer         :: weno_order_pv        = 5

  character(8)    :: pgf_scheme           = ''       ! lin97, ptb

  character(30)   :: bg_adv_scheme        = 'ffsl'
  character(30)   :: pt_adv_scheme        = 'upwind'
  character(30)   :: nh_adv_scheme        = 'upwind'
  character(8)    :: limiter_type         = 'mono'
  character(8)    :: ffsl_flux_type       = 'ppm'
  character(8)    :: tvd_limiter_type     = 'van_leer'

  character(8)    :: zonal_tridiag_solver = 'spk'   ! mkl, spk

  integer         :: weno_order           = 5       ! 3, 5
  integer         :: weno_order_h         = 5       ! 3, 5
  integer         :: weno_order_v         = 5       ! 3, 5
  integer         :: upwind_order         = 5       ! 0, 1, 3, 5
  integer         :: upwind_order_h       = -1      ! 0, 1, 3, 5
  integer         :: upwind_order_v       = -1      ! 0, 1, 3, 5
  real(r8)        :: upwind_wgt           = 1

  character(30)   :: time_scheme          = 'wrfrk3'
  logical         :: save_dyn_calc        = .true.

  ! Filter settings
  real(r8)        :: filter_wave_speed    = 300.0_r8
  real(r8)        :: filter_coef_a        = 3.5_r8
  real(r8)        :: filter_coef_b        = 0.5_r8
  real(r8)        :: filter_coef_c        = 0.3_r8
  real(r8)        :: filter_gauss_sigma   = 8.0_r8
  real(r8)        :: filter_min_width     = 0.0_r8

  ! Damping settings
  logical         :: use_topo_smooth      = .false.
  logical         :: use_zs_polar_filter  = .false.
  real(r8)        :: topo_max_slope       = 0.12_r8
  integer         :: topo_smooth_order    = 2
  real(r8)        :: topo_smooth_coef     = 1.0e4_r8
  integer         :: topo_smooth_cycles   = 500
  logical         :: use_div_damp         = .false.
  integer         :: div_damp_cycles      = 1
  integer         :: div_damp_order       = 2
  real(r8)        :: div_damp_top         = 1
  integer         :: div_damp_k0          = 10
  real(r8)        :: div_damp_pole        = 100
  real(r8)        :: div_damp_lat0        = 80
  real(r8)        :: div_damp_coef2       = 1.0_r8 / 128.0_r8
  real(r8)        :: div_damp_coef4       = 0.01_r8
  logical         :: use_vor_damp         = .false.
  integer         :: vor_damp_cycles      = 1
  integer         :: vor_damp_order       = 2
  real(r8)        :: vor_damp_coef2       = 0.05_r8
  real(r8)        :: vor_damp_lat0        = 80
  logical         :: use_rayleigh_damp_w  = .false.
  real(r8)        :: rayleigh_damp_w_coef = 0.2       ! s-1
  real(r8)        :: rayleigh_damp_top    = 10.0d3    ! m
  logical         :: use_p_damp           = .false.
  real(r8)        :: p_damp_coef          = 0.12_r8
  logical         :: use_smag_damp        = .false.
  integer         :: smag_damp_cycles     = 1
  real(r8)        :: smag_damp_coef       = 0.015
  logical         :: use_laplace_damp     = .false.
  integer         :: laplace_damp_order   = 4
  real(r8)        :: laplace_damp_coef    = 1.0e12_r8
  logical         :: use_sponge_layer     = .false.
  integer         :: sponge_layer_k0      = 6
  real(r8)        :: sponge_layer_coef    = 1.0e6_r8

  ! Input settings
  integer         :: input_ngroups        = 0

  ! Output settings
#if (REAL_KIND == 4)
  character(8)    :: output_i0_dtype      = 'r4'
#elif (REAL_KIND == 8)
  character(8)    :: output_i0_dtype      = 'r8'
#endif
  logical         :: output_h0            = .true.
  logical         :: append_h0            = .true.
  character(8)    :: output_h0_dtype      = 'r4'
  logical         :: output_h1            = .false.
  logical         :: output_h2            = .false.
  character(30)   :: output_h0_new_file   = ''
  character(8)    :: output_h0_vars(100)  = ''
  integer         :: output_ngroups       = 0

  integer         :: output_nlev          = 0
  real(r8)        :: output_plev_hPa(100) = 0

  namelist /gmcore_control/     &
    planet                    , &
    use_aqua_planet           , &
    case_name                 , &
    test_case                 , &
    case_desc                 , &
    nlon                      , &
    nlat                      , &
    nlev                      , &
    nonhydrostatic            , &
    ideal_dry_core            , &
    advection                 , &
    nproc_io                  , &
    nproc_x                   , &
    nproc_y                   , &
    lon_hw                    , &
    lat_hw                    , &
    proc_layout               , &
    use_async_io              , &
    proc_io_stride            , &
    initial_time              , &
    start_time                , &
    end_time                  , &
    dt_dyn                    , &
    dt_adv                    , &
    dt_phys                   , &
    pdc_type                  , &
    run_years                 , &
    run_my                    , &
    run_sol                   , &
    run_hours                 , &
    run_days                  , &
    history_interval          , &
    restart_interval          , &
    print_interval            , &
    initial_file              , &
    output_i0_dtype           , &
    restart_file              , &
    restart                   , &
    prepare_regrid_gz         , &
    init_hydrostatic_gz       , &
    topo_file                 , &
    topo_type                 , &
    bkg_file                  , &
    bkg_type                  , &
    tangent_wgt_scheme        , &
    implicit_w_wgt            , &
    vert_coord_scheme         , &
    vert_coord_template       , &
    refer_state_scheme        , &
    ptop                      , &
    hybrid_coord_p0           , &
    tiso                      , &
    dzbot                     , &
    dzmax                     , &
    dzstretch_s               , &
    dzstretch_u               , &
    eta_b                     , &
    hybrid_coord_ncep_psig    , &
    hybrid_coord_ncep_ppre    , &
    hybrid_coord_ncep_dpbot   , &
    hybrid_coord_ncep_dpsig   , &
    hybrid_coord_ncep_dppre   , &
    hybrid_coord_ncep_dptop   , &
    ke_scheme                 , &
    ke_cell_wgt               , &
    pv_adv_scheme             , &
    pv_pole_stokes            , &
    upwind_order_pv           , &
    upwind_wgt_pv             , &
    weno_order_pv             , &
    pgf_scheme                , &
    bg_adv_scheme             , &
    pt_adv_scheme             , &
    nh_adv_scheme             , &
    limiter_type              , &
    ffsl_flux_type            , &
    tvd_limiter_type          , &
    zonal_tridiag_solver      , &
    weno_order                , &
    weno_order_h              , &
    weno_order_v              , &
    upwind_order              , &
    upwind_order_h            , &
    upwind_order_v            , &
    upwind_wgt                , &
    time_scheme               , &
    save_dyn_calc             , &
    filter_wave_speed         , &
    filter_coef_a             , &
    filter_coef_b             , &
    filter_coef_c             , &
    filter_gauss_sigma        , &
    filter_min_width          , &
    physics_suite             , &
    mp_scheme                 , &
    pbl_scheme                , &
    cam_namelist_path         , &
    filter_ptend              , &
    gmcore_data_dir           , &
    use_topo_smooth           , &
    use_zs_polar_filter       , &
    topo_max_slope            , &
    topo_smooth_order         , &
    topo_smooth_coef          , &
    topo_smooth_cycles        , &
    use_div_damp              , &
    div_damp_cycles           , &
    div_damp_order            , &
    div_damp_coef2            , &
    div_damp_coef4            , &
    div_damp_k0               , &
    div_damp_top              , &
    div_damp_pole             , &
    div_damp_lat0             , &
    use_vor_damp              , &
    vor_damp_cycles           , &
    vor_damp_order            , &
    vor_damp_coef2            , &
    vor_damp_lat0             , &
    use_rayleigh_damp_w       , &
    rayleigh_damp_w_coef      , &
    rayleigh_damp_top         , &
    use_p_damp                , &
    p_damp_coef               , &
    use_smag_damp             , &
    smag_damp_cycles          , &
    smag_damp_coef            , &
    use_laplace_damp          , &
    laplace_damp_order        , &
    laplace_damp_coef         , &
    use_sponge_layer          , &
    sponge_layer_k0           , &
    sponge_layer_coef         , &
    input_ngroups             , &
    output_h0                 , &
    append_h0                 , &
    output_h0_dtype           , &
    output_h1                 , &
    output_h2                 , &
    output_h0_new_file        , &
    output_h0_vars            , &
    output_ngroups            , &
    output_nlev              , &
    output_plev_hPa

contains

  subroutine parse_namelist(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path, status='old')
    read(10, nml=gmcore_control)
    close(10)

    if (use_async_io .and. nproc_io < 1) then
      nproc_io = 1
    end if

    if (use_async_io) then
      input_ngroups = nproc_io
      output_ngroups = nproc_io
    end if

    if (.not. restart) then
      append_h0 = .false.
    end if

    ! Here we set baroclinic according to levels.
    baroclinic = nlev > 1
    if (.not. baroclinic) then
      hydrostatic    = .false.
      nonhydrostatic = .false.
      ke_scheme      = 1
    else
      hydrostatic = .not. nonhydrostatic
    end if

    if (advection) then
      hydrostatic    = .false.
      baroclinic     = .false.
      nonhydrostatic = .false.
    end if

    ! Set default for nonhydrostatic temporally.
    if (pgf_scheme == '') then
      if (nonhydrostatic) then
        pgf_scheme   = 'ptb'
      else
        pgf_scheme   = 'lin97'
      end if
    end if

    if (upwind_order_h == -1) upwind_order_h = upwind_order
    if (upwind_order_v == -1) upwind_order_v = upwind_order

    if (dt_dyn  == 0) dt_dyn  = dt_adv
    if (dt_adv  == 0) dt_adv  = dt_dyn
    if (dt_phys == 0) dt_phys = dt_adv

    ! Convert time step sizes to Earth time system
    select case (planet)
    case ('mars')
      time_scale = mars_sol_seconds / earth_day_seconds
    end select
    dt_dyn  = dt_dyn  * time_scale
    dt_adv  = dt_adv  * time_scale
    dt_phys = dt_phys * time_scale

    if (physics_suite == 'N/A') then
      pdc_type       = 0
    end if

    if (.not. use_div_damp) then
      div_damp_order = 0
    end if

    if (.not. use_vor_damp) then
      vor_damp_order = 0
    end if

  end subroutine parse_namelist

  subroutine print_namelist()

      write(*, *) '=================== GMCORE Parameters ==================='
      write(*, *) 'case_name           = ', trim(case_name)
      write(*, *) 'nlon                = ', to_str(nlon)
      write(*, *) 'nlat                = ', to_str(nlat)
      write(*, *) 'nlev                = ', to_str(nlev)
      write(*, *) 'physics_suite       = ', trim(physics_suite)
      write(*, *) 'mp_scheme           = ', trim(mp_scheme)
      write(*, *) 'pbl_scheme          = ', trim(pbl_scheme)
      write(*, *) 'hydrostatic         = ', to_str(hydrostatic)
      write(*, *) 'nonhydrostatic      = ', to_str(nonhydrostatic)
      write(*, *) 'ideal_dry_core      = ', to_str(ideal_dry_core)
      write(*, *) 'vert_coord_scheme   = ', trim(vert_coord_scheme)
      write(*, *) 'vert_coord_template = ', trim(vert_coord_template)
      write(*, *) 'ptop                = ', to_str(ptop, 4)
      write(*, *) 'hybrid_coord_p0     = ', to_str(hybrid_coord_p0, 2)
      write(*, *) 'dt_dyn              = ', to_str(dt_dyn , 2)
      write(*, *) 'dt_adv              = ', to_str(dt_adv , 2)
      write(*, *) 'dt_phys             = ', to_str(dt_phys, 2)
      write(*, *) 'time_scheme         = ', trim(time_scheme)
      write(*, *) 'save_dyn_calc       = ', to_str(save_dyn_calc)
      write(*, *) 'pdc_type            = ', to_str(pdc_type)
      write(*, *) 'filter_wave_speed   = ', filter_wave_speed
      write(*, *) 'filter_coef_a       = ', filter_coef_a
      write(*, *) 'filter_coef_b       = ', filter_coef_b
      write(*, *) 'filter_coef_c       = ', filter_coef_c
      write(*, *) 'filter_gauss_sigma  = ', filter_gauss_sigma
      write(*, *) 'filter_min_width    = ', filter_min_width
      write(*, *) 'filter_ptend        = ', to_str(filter_ptend)
      write(*, *) 'pgf_scheme          = ', trim(pgf_scheme)
      write(*, *) 'bg_adv_scheme       = ', trim(bg_adv_scheme)
      write(*, *) 'pt_adv_scheme       = ', trim(pt_adv_scheme)
      write(*, *) 'nh_adv_scheme       = ', trim(nh_adv_scheme)
      write(*, *) 'limiter_type        = ', trim(limiter_type)
    if (pt_adv_scheme == 'ffsl') then
      write(*, *) 'ffsl_flux_type      = ', trim(ffsl_flux_type)
    end if
      write(*, *) 'ke_scheme           = ', to_str(ke_scheme)
    if (ke_scheme == 2) then
      write(*, *) 'ke_cell_wgt         = ', to_str(ke_cell_wgt, 2)
    end if
      write(*, *) 'pv_adv_scheme       = ', trim(pv_adv_scheme)
      write(*, *) 'pv_pole_stokes      = ', to_str(pv_pole_stokes)
    if (pv_adv_scheme == 'upwind') then
      write(*, *) 'upwind_order_pv     = ', to_str(upwind_order_pv)
      write(*, *) 'upwind_wgt_pv       = ', to_str(upwind_wgt_pv, 2)
    else if (pv_adv_scheme == 'weno') then
      write(*, *) 'weno_order_pv       = ', to_str(weno_order_pv)
    end if
    if (pt_adv_scheme == 'upwind' .or. nh_adv_scheme == 'upwind') then
      write(*, *) 'upwind_order_h      = ', to_str(upwind_order_h)
      write(*, *) 'upwind_order_v      = ', to_str(upwind_order_v)
      write(*, *) 'upwind_wgt          = ', to_str(upwind_wgt, 4)
    else if (pt_adv_scheme == 'weno' .or. nh_adv_scheme == 'weno') then
      write(*, *) 'weno_order_h        = ', to_str(weno_order_h)
      write(*, *) 'weno_order_v        = ', to_str(weno_order_v)
    end if
      write(*, *) 'use_topo_smooth     = ', to_str(use_topo_smooth)
    if (use_topo_smooth) then
      write(*, *) 'topo_smooth_cycles  = ', to_str(topo_smooth_cycles)
    end if
      write(*, *) 'use_div_damp        = ', to_str(use_div_damp)
    if (use_div_damp) then
      write(*, *) 'div_damp_cycles     = ', to_str(div_damp_cycles)
      write(*, *) 'div_damp_order      = ', to_str(div_damp_order)
      write(*, *) 'div_damp_coef2      = ', div_damp_coef2
      write(*, *) 'div_damp_coef4      = ', div_damp_coef4
      write(*, *) 'div_damp_top        = ', to_str(div_damp_top, 3)
      write(*, *) 'div_damp_pole       = ', to_str(div_damp_pole, 3)
      write(*, *) 'div_damp_lat0       = ', to_str(div_damp_lat0, 3)
    end if
      write(*, *) 'use_vor_damp        = ', to_str(use_vor_damp)
    if (use_vor_damp) then
      write(*, *) 'vor_damp_cycles     = ', to_str(vor_damp_cycles)
      write(*, *) 'vor_damp_order      = ', to_str(vor_damp_order)
      write(*, *) 'vor_damp_coef2      = ', vor_damp_coef2
      write(*, *) 'vor_damp_lat0       = ', to_str(vor_damp_lat0, 3)
    end if
    if (nonhydrostatic) then
      write(*, *) 'implicit_w_wgt      = ', to_str(implicit_w_wgt, 3)
      write(*, *) 'use_rayleigh_damp_w = ', to_str(use_rayleigh_damp_w)
    end if
    if (use_rayleigh_damp_w) then
      write(*, *) 'rayleigh_damp_w_coef= ', rayleigh_damp_w_coef
      write(*, *) 'rayleigh_damp_top   = ', to_str(rayleigh_damp_top   , 2)
    end if
      write(*, *) 'use_p_damp          = ', to_str(use_p_damp)
    if (use_p_damp) then
      write(*, *) 'p_damp_coef         = ', p_damp_coef
    end if
      write(*, *) 'use_smag_damp       = ', to_str(use_smag_damp)
    if (use_smag_damp) then
      write(*, *) 'smag_damp_cycles    = ', to_str(smag_damp_cycles)
      write(*, *) 'smag_damp_coef      = ', smag_damp_coef
    end if
      write(*, *) 'use_laplace_damp    = ', to_str(use_laplace_damp)
    if (use_laplace_damp) then
      write(*, *) 'laplace_damp_order  = ', to_str(laplace_damp_order)
      write(*, *) 'laplace_damp_coef   = ', laplace_damp_coef
    end if
      write(*, *) 'use_sponge_layer    = ', to_str(use_sponge_layer)
    if (use_sponge_layer) then
      write(*, *) 'sponge_layer_k0     = ', to_str(sponge_layer_k0)
      write(*, *) 'sponge_layer_coef   = ', sponge_layer_coef
    end if
      write(*, *) '========================================================='

  end subroutine print_namelist

  subroutine namelist_add_atts(tag)

    use fiona

    character(*), intent(in) :: tag

    call fiona_add_att(tag, 'planet', planet)
    call fiona_add_att(tag, 'nonhydrostatic', nonhydrostatic)
    call fiona_add_att(tag, 'case_name', case_name)
    call fiona_add_att(tag, 'time_scheme', time_scheme)
    call fiona_add_att(tag, 'dt_dyn', dt_dyn / time_scale)
    call fiona_add_att(tag, 'dt_adv', dt_adv / time_scale)
    call fiona_add_att(tag, 'dt_phys', dt_phys / time_scale)
    call fiona_add_att(tag, 'pt_adv_scheme', pt_adv_scheme)
    call fiona_add_att(tag, 'nh_adv_scheme', nh_adv_scheme)
    call fiona_add_att(tag, 'physics_suite', physics_suite)
    call fiona_add_att(tag, 'filter_coef_a', filter_coef_a)
    call fiona_add_att(tag, 'filter_coef_b', filter_coef_b)
    call fiona_add_att(tag, 'filter_coef_c', filter_coef_c)
    call fiona_add_att(tag, 'use_div_damp', merge(1, 0, use_div_damp))
    if (use_div_damp) then
      call fiona_add_att(tag, 'div_damp_coef2', div_damp_coef2)
      call fiona_add_att(tag, 'div_damp_top', div_damp_top)
      call fiona_add_att(tag, 'div_damp_k0', div_damp_k0)
      call fiona_add_att(tag, 'div_damp_pole', div_damp_pole)
      call fiona_add_att(tag, 'div_damp_lat0', div_damp_lat0)
    end if
    call fiona_add_att(tag, 'use_vor_damp', merge(1, 0, use_vor_damp))
    if (use_vor_damp) then
      call fiona_add_att(tag, 'vor_damp_coef2', vor_damp_coef2)
      call fiona_add_att(tag, 'vor_damp_lat0', vor_damp_lat0)
    end if
    call fiona_add_att(tag, 'use_smag_damp', merge(1, 0, use_smag_damp))
    if (use_smag_damp) then
      call fiona_add_att(tag, 'smag_damp_coef', smag_damp_coef)
    end if
    call fiona_add_att(tag, 'use_laplace_damp', merge(1, 0, use_laplace_damp))
    if (use_laplace_damp) then
      call fiona_add_att(tag, 'laplace_damp_order', laplace_damp_order)
      call fiona_add_att(tag, 'laplace_damp_coef', laplace_damp_coef)
    end if
    call fiona_add_att(tag, 'use_sponge_layer', merge(1, 0, use_sponge_layer))
    if (use_sponge_layer) then
      call fiona_add_att(tag, 'sponge_layer_k0', sponge_layer_k0)
      call fiona_add_att(tag, 'sponge_layer_coef', sponge_layer_coef)
    end if

  end subroutine namelist_add_atts

end module namelist_mod
