! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module dynamics_types_mod

  use container
  use const_mod
  use namelist_mod
  use latlon_field_types_mod
  use tracer_types_mod

  implicit none

  private

  public dstate_type
  public dtend_type
  public static_type
  public aux_array_type

  ! NOTE:
  !   Variables with '_lon', '_lat' and '_lev' are on the half grids on the corresponding direction.
  type dstate_type
    integer :: id = 0
    type(array_type) fields
    type(latlon_field3d_type), pointer :: u       => null()
    type(latlon_field3d_type), pointer :: v       => null()
    type(latlon_field3d_type), pointer :: u_lon   => null()
    type(latlon_field3d_type), pointer :: v_lat   => null()
    type(latlon_field3d_type), pointer :: we_lev  => null()
    type(latlon_field3d_type), pointer :: gz      => null()
    type(latlon_field3d_type), pointer :: gz_lev  => null()
    type(latlon_field3d_type), pointer :: dmg     => null()
    type(latlon_field3d_type), pointer :: dmg_lev => null()
    type(latlon_field3d_type), pointer :: pt      => null()
    type(latlon_field3d_type), pointer :: t       => null()
    type(latlon_field3d_type), pointer :: tv      => null()
    type(latlon_field3d_type), pointer :: mg      => null()
    type(latlon_field3d_type), pointer :: mg_lev  => null()
    type(latlon_field2d_type), pointer :: mgs     => null()
    type(latlon_field3d_type), pointer :: ph      => null()
    type(latlon_field3d_type), pointer :: ph_lev  => null()
    type(latlon_field2d_type), pointer :: phs     => null()
    type(latlon_field3d_type), pointer :: rhod    => null()
    ! Nonhydrostatic variable
    type(latlon_field3d_type), pointer :: we      => null()
    type(latlon_field3d_type), pointer :: w       => null()
    type(latlon_field3d_type), pointer :: w_lev   => null()
    type(latlon_field3d_type), pointer :: p       => null()
    type(latlon_field3d_type), pointer :: p_lev   => null()
    type(latlon_field2d_type), pointer :: ps      => null()
    ! Total diagnostics
    real(r8) tm
    real(r8) te, te_ke, te_ie, te_pe
    real(r8) tpe
  contains
    procedure :: init  => dstate_init
    procedure :: clear => dstate_clear
    procedure :: c2a   => dstate_c2a
    procedure :: a2c   => dstate_a2c
    final dstate_final
  end type dstate_type

  type dtend_type
    type(array_type) fields
    type(latlon_field3d_type), pointer :: du    => null()
    type(latlon_field3d_type), pointer :: dv    => null()
    type(latlon_field3d_type), pointer :: dgz   => null()
    type(latlon_field3d_type), pointer :: dpt   => null()
    type(latlon_field2d_type), pointer :: dmgs  => null()
#ifdef OUTPUT_H1_DTEND
    type(latlon_field3d_type), pointer :: dudt_coriolis => null()
    type(latlon_field3d_type), pointer :: dvdt_coriolis => null()
    type(latlon_field3d_type), pointer :: dudt_wedudeta => null()
    type(latlon_field3d_type), pointer :: dvdt_wedvdeta => null()
    type(latlon_field3d_type), pointer :: dudt_dkedx    => null()
    type(latlon_field3d_type), pointer :: dvdt_dkedy    => null()
    type(latlon_field3d_type), pointer :: dudt_pgf      => null()
    type(latlon_field3d_type), pointer :: dvdt_pgf      => null()
#endif
    logical :: update_u   = .false.
    logical :: update_v   = .false.
    logical :: update_gz  = .false.
    logical :: update_pt  = .false.
    logical :: update_mgs = .false.
  contains
    procedure :: init        => dtend_init
    procedure :: reset_flags => dtend_reset_flags
    procedure :: clear       => dtend_clear
    final dtend_final
  end type dtend_type

  type static_type
    type(array_type) fields
    type(latlon_field2d_type), pointer :: landmask    => null()
    ! Topography
    type(latlon_field2d_type), pointer :: gzs         => null()
    type(latlon_field2d_type), pointer :: zs_std      => null()
    type(latlon_field2d_type), pointer :: dzsdx       => null()
    type(latlon_field2d_type), pointer :: dzsdy       => null()
    ! Reference surface pressure
    type(latlon_field2d_type), pointer :: ref_ps      => null()
    type(latlon_field2d_type), pointer :: ref_ps_smth => null()
    type(latlon_field2d_type), pointer :: ref_ps_perb => null()
  contains
    procedure :: init_stage1 => static_init_stage1
    procedure :: init_stage2 => static_init_stage2
    procedure :: clear       => static_clear
    final static_final
  end type static_type

  type aux_array_type
    type(array_type) fields
    ! Smagorinsky damping variables
    type(latlon_field3d_type), pointer :: smag_t      => null() ! tension strain
    type(latlon_field3d_type), pointer :: smag_s      => null() ! shear strain on vertex
    type(latlon_field3d_type), pointer :: kmh         => null() ! nonlinear diffusion coef
    type(latlon_field3d_type), pointer :: kmh_lon     => null() ! nonlinear diffusion coef on zonal edge
    type(latlon_field3d_type), pointer :: kmh_lat     => null() ! nonlinear diffusion coef on meridional edge
    ! Other variables
    type(latlon_field3d_type), pointer :: v_lon       => null() ! Meridional wind speed at lon edge (m s-1)
    type(latlon_field3d_type), pointer :: u_lat       => null() ! Zonal wind speed at lat edge (m s-1)
    type(latlon_field3d_type), pointer :: ke          => null() ! Kinetic energy
    type(latlon_field3d_type), pointer :: pv_lon      => null() ! Potential vorticity on zonal edge
    type(latlon_field3d_type), pointer :: pv_lat      => null() ! Potential vorticity on merdional edge
    type(latlon_field3d_type), pointer :: dmg_lon     => null() ! Mass on zonal edge
    type(latlon_field3d_type), pointer :: dmg_lat     => null() ! Mass on merdional edge
    type(latlon_field3d_type), pointer :: dmg_vtx     => null() ! Mass on vertex
    type(latlon_field3d_type), pointer :: pkh_lev     => null() ! Exner pressure on half levels
    type(latlon_field3d_type), pointer :: we_lev_lon  => null() ! Vertical coordinate speed multiplied by 𝛛π/𝛛η on zonal edge
    type(latlon_field3d_type), pointer :: we_lev_lat  => null() ! Vertical coordinate speed multiplied by 𝛛π/𝛛η on merdional edge
    type(latlon_field3d_type), pointer :: ptfx        => null() ! Potential temperature on the zonal edge
    type(latlon_field3d_type), pointer :: ptfy        => null() ! Potential temperature on the merdional edge
    type(latlon_field3d_type), pointer :: ptfz        => null() ! Potential temperature on the vertical edge
    type(latlon_field3d_type), pointer :: mfx_lon     => null() ! Normal mass flux on zonal edge
    type(latlon_field3d_type), pointer :: mfy_lat     => null() ! Normal mass flux on merdional edge
    type(latlon_field3d_type), pointer :: mfx_lat     => null() ! Tangient mass flux on zonal edge
    type(latlon_field3d_type), pointer :: mfy_lon     => null() ! Tangient mass flux on merdional edge
    type(latlon_field3d_type), pointer :: vor         => null() ! Vorticity (s-1)
    type(latlon_field3d_type), pointer :: pv          => null() ! Potential vorticity
    type(latlon_field3d_type), pointer :: div         => null() ! Divergence (s-1)
    type(latlon_field3d_type), pointer :: div2        => null() ! Laplacian of divergence (s-1)
    type(latlon_field3d_type), pointer :: dmf         => null() ! Mass flux divergence on full level (Pa s-1)
    type(latlon_field3d_type), pointer :: dmf_lev     => null() ! Mass flux divergence on half level (Pa s-1)
    type(latlon_field3d_type), pointer :: omg         => null() ! Vertical pressure velocity (Pa s-1)
    ! Tendencies from physics
    type(latlon_field3d_type), pointer :: dudt_phys   => null()
    type(latlon_field3d_type), pointer :: dvdt_phys   => null()
    type(latlon_field3d_type), pointer :: dptdt_phys  => null()
    type(latlon_field4d_type), pointer :: dqdt_phys   => null()
    ! Perturbed quantities for calculating HPGF
    type(latlon_field3d_type), pointer :: p_ptb       => null()
    type(latlon_field3d_type), pointer :: gz_ptb      => null()
    type(latlon_field3d_type), pointer :: dp_ptb      => null()
    type(latlon_field3d_type), pointer :: ad_ptb      => null()
    ! Nonhydrostatic variables
    type(latlon_field3d_type), pointer :: u_lev_lon   => null()
    type(latlon_field3d_type), pointer :: v_lev_lat   => null()
    type(latlon_field3d_type), pointer :: mfx_lev_lon => null()
    type(latlon_field3d_type), pointer :: mfy_lev_lat => null()
    type(latlon_field3d_type), pointer :: adv_w_lev   => null()
    type(latlon_field3d_type), pointer :: adv_gz_lev  => null()
  contains
    procedure :: init      => aux_array_init
    procedure :: init_phys => aux_array_init_phys
    procedure :: clear     => aux_array_clear
    final aux_array_final
  end type aux_array_type

  interface add_var
    module procedure add_var2d
    module procedure add_var3d
    module procedure add_var4d
  end interface add_var

contains

  subroutine add_var2d(fields, name, long_name, units, loc, mesh, halo, output, restart, ptr, link)

    type(array_type), intent(inout) :: fields
    character(*), intent(in) :: name
    character(*), intent(in) :: long_name
    character(*), intent(in) :: units
    character(*), intent(in) :: loc
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)
    character(*), intent(in) :: output
    logical, intent(in) :: restart
    type(latlon_field2d_type), intent(inout), pointer :: ptr
    type(latlon_field2d_type), intent(in), optional :: link

    if (associated(ptr)) deallocate(ptr)
    allocate(ptr)
    call fields%append_ptr(ptr)
    call ptr%init(name, long_name, units, loc, mesh, halo, link=link, output=output, restart=restart)

  end subroutine add_var2d

  subroutine add_var3d(fields, name, long_name, units, loc, mesh, halo, output, restart, ptr, halo_cross_pole, link)

    type(array_type), intent(inout) :: fields
    character(*), intent(in) :: name
    character(*), intent(in) :: long_name
    character(*), intent(in) :: units
    character(*), intent(in) :: loc
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)
    character(*), intent(in) :: output
    logical, intent(in) :: restart
    type(latlon_field3d_type), intent(inout), pointer :: ptr
    logical, intent(in), optional :: halo_cross_pole
    type(latlon_field3d_type), intent(in), optional :: link

    if (associated(ptr)) deallocate(ptr)
    allocate(ptr)
    call fields%append_ptr(ptr)
    call ptr%init(name, long_name, units, loc, mesh, halo, halo_cross_pole=halo_cross_pole, link=link, output=output, restart=restart)

  end subroutine add_var3d

  subroutine add_var4d(fields, name, long_name, units, loc, mesh, dim4_name, dim4_size, var4_names, halo, output, restart, ptr, halo_cross_pole)

    type(array_type), intent(inout) :: fields
    character(*), intent(in) :: name
    character(*), intent(in) :: long_name
    character(*), intent(in) :: units
    character(*), intent(in) :: loc
    type(latlon_mesh_type), intent(in) :: mesh
    character(*), intent(in) :: dim4_name
    integer, intent(in) :: dim4_size
    character(*), intent(in) :: var4_names(:)
    type(latlon_halo_type), intent(in) :: halo(:)
    character(*), intent(in) :: output
    logical, intent(in) :: restart
    type(latlon_field4d_type), intent(inout), pointer :: ptr
    logical, intent(in), optional :: halo_cross_pole

    if (associated(ptr)) deallocate(ptr)
    allocate(ptr)
    call fields%append_ptr(ptr)
    call ptr%init(name, long_name, units, loc, mesh, halo, halo_cross_pole=halo_cross_pole, &
      dim4_name=dim4_name, dim4_size=dim4_size, var4_names=var4_names, output=output, restart=restart)

  end subroutine add_var4d

  subroutine dstate_init(this, filter_mesh, filter_halo, mesh, halo)

    class(dstate_type), intent(inout), target :: this
    type(latlon_mesh_type), intent(in) :: filter_mesh
    type(latlon_halo_type), intent(in) :: filter_halo(:)
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)

    character(field_name_len     ) name
    character(field_long_name_len) long_name
    character(field_units_len    ) units

    call this%clear()

    this%fields = array(30)

    if (.not. advection) then
      call add_var(this%fields                                               , &
        name            ='u'                                                 , &
        long_name       ='Zonal wind component'                              , &
        units           ='m s-1'                                             , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.false.                                             , &
        ptr             =this%u                                              )
    end if
    if (.not. advection) then
      call add_var(this%fields                                               , &
        name            ='v'                                                 , &
        long_name       ='Meridional wind component'                         , &
        units           ='m s-1'                                             , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.false.                                             , &
        ptr             =this%v                                              )
    end if
    call add_var(this%fields                                                 , &
      name              ='u_lon'                                             , &
      long_name         ='Zonal wind component'                              , &
      units             ='m s-1'                                             , &
      loc               ='lon'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.true.                                              , &
      ptr               =this%u_lon                                          )
    call add_var(this%fields                                                 , &
      name              ='v_lat'                                             , &
      long_name         ='Meridional wind component'                         , &
      units             ='m s-1'                                             , &
      loc               ='lat'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.true.                                              , &
      ptr               =this%v_lat                                          )
    if (baroclinic .or. advection) then
      call add_var(this%fields                                               , &
        name            ='we_lev'                                            , &
        long_name       ='Vertical mass flux'                                , &
        units           ='Pa s-1'                                            , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%we_lev                                         )
    end if
    call add_var(this%fields                                                 , &
      name              ='gz'                                                , &
      long_name         ='Geopotential'                                      , &
      units             ='m2 s-2'                                            , &
      loc               ='cell'                                              , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h0'                                                , &
      restart           =.not. baroclinic                                    , & ! FIXME: Revise this.
      ptr               =this%gz                                             )
    if (nonhydrostatic) then
      call add_var(this%fields                                               , &
        name            ='gz_lev'                                            , &
        long_name       ='Geopotential'                                      , &
        units           ='m2 s-2'                                            , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.true.                                              , &
        ptr             =this%gz_lev                                         , &
        halo_cross_pole =.true.                                              )
    else
      call add_var(this%fields                                               , &
        name            ='gz_lev'                                            , &
        long_name       ='Geopotential'                                      , &
        units           ='m2 s-2'                                            , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%gz_lev                                         )
    end if
    call add_var(this%fields                                                 , &
      name              ='dmg'                                               , &
      long_name         ='Dry-air weight'                                    , &
      units             ='Pa'                                                , &
      loc               ='cell'                                              , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      ptr               =this%dmg                                            )
    if (baroclinic .or. advection) then
      call add_var(this%fields                                               , &
        name            ='dmg_lev'                                           , &
        long_name       ='Dry-air weight'                                    , &
        units           ='Pa'                                                , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%dmg_lev                                        )
    end if
    if (baroclinic .or. advection) then
      call add_var(this%fields                                               , &
        name            ='pt'                                                , &
        long_name       ='Modified potential temperature'                    , &
        units           ='K'                                                 , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.true.                                              , &
        ptr             =this%pt                                             , &
        halo_cross_pole =.true.                                              )
    end if
    if (baroclinic .or. advection) then
      call add_var(this%fields                                               , &
        name            ='t'                                                 , &
        long_name       ='Temperature'                                       , &
        units           ='K'                                                 , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.false.                                             , &
        ptr             =this%t                                              )
    end if
    if (baroclinic .or. advection) then
      call add_var(this%fields                                               , &
        name            ='tv'                                                , &
        long_name       ='Virtual temperature'                               , &
        units           ='K'                                                 , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%tv                                             )
    end if
    call add_var(this%fields                                                 , &
      name              ='mg'                                                , &
      long_name         ='Dry-air weight'                                    , &
      units             ='Pa'                                                , &
      loc               ='cell'                                              , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      ptr               =this%mg                                             )
    if (baroclinic .or. advection) then
      call add_var(this%fields                                               , &
        name            ='mg_lev'                                            , &
        long_name       ='Dry-air weight'                                    , &
        units           ='Pa'                                                , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%mg_lev                                         )
    end if
    if (baroclinic .or. advection) then
      call add_var(this%fields                                               , &
        name            ='mgs'                                               , &
        long_name       ='Dry-air weight on surface'                         , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.true.                                              , &
        ptr             =this%mgs                                            )
    end if
    if (baroclinic .or. advection) then
      call add_var(this%fields                                               , &
        name            ='ph'                                                , &
        long_name       ='Hydrostatic pressure'                              , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%ph                                             )
    end if
    if (baroclinic .or. advection) then
      call add_var(this%fields                                               , &
        name            ='ph_lev'                                            , &
        long_name       ='Hydrostatic pressure'                              , &
        units           ='Pa'                                                , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%ph_lev                                         )
    end if
    if (baroclinic) then
      call add_var(this%fields                                               , &
        name            ='phs'                                               , &
        long_name       ='Hydrostatic pressure on surface'                   , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.true.                                              , &
        ptr             =this%phs                                            )
      call this%phs%link(this%ph_lev, mesh%half_nlev)
    end if
    if (baroclinic) then
      call add_var(this%fields                                               , &
        name            ='rhod'                                              , &
        long_name       ='Dry-air density'                                   , &
        units           ='kg m-3'                                            , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.false.                                             , &
        ptr             =this%rhod                                           )
    end if
    if (nonhydrostatic) then
      call add_var(this%fields                                               , &
        name            ='we'                                                , &
        long_name       ='Vertical mass flux'                                , &
        units           ='Pa s-1'                                            , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%we                                             )
    end if
    if (nonhydrostatic) then
      call add_var(this%fields                                               , &
        name            ='w'                                                 , &
        long_name       ='Vertical wind speed'                               , &
        units           ='m s-1'                                             , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%w                                              )
    end if
    if (nonhydrostatic) then
      call add_var(this%fields                                               , &
        name            ='w_lev'                                             , &
        long_name       ='Vertical wind speed'                               , &
        units           ='m s-1'                                             , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.true.                                              , &
        ptr             =this%w_lev                                          , &
        halo_cross_pole =.true.                                              )
    end if
    if (hydrostatic) then
      call add_var(this%fields                                               , &
        name            ='p'                                                 , &
        long_name       ='Pressure'                                          , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%p                                              , &
        link            =this%ph                                             )
    else if (nonhydrostatic) then
      call add_var(this%fields                                               , &
        name            ='p'                                                 , &
        long_name       ='Pressure'                                          , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%p                                              )
    end if
    if (hydrostatic) then
      call add_var(this%fields                                               , &
        name            ='p_lev'                                             , &
        long_name       ='Pressure'                                          , &
        units           ='Pa'                                                , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%p_lev                                          , &
        link            =this%ph_lev                                         )
    else if (nonhydrostatic) then
      call add_var(this%fields                                               , &
        name            ='p_lev'                                             , &
        long_name       ='Pressure'                                          , &
        units           ='Pa'                                                , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%p_lev                                          )
    end if
    if (hydrostatic) then
      call add_var(this%fields                                               , &
        name            ='ps'                                                , &
        long_name       ='Surface pressure'                                  , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%ps                                             , &
        link            =this%phs                                            )
    else if (nonhydrostatic) then
      call add_var(this%fields                                               , &
        name            ='ps'                                                , &
        long_name       ='Surface pressure'                                  , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%ps                                             )
    end if

  end subroutine dstate_init

  subroutine dstate_clear(this)

    class(dstate_type), intent(inout) :: this

    class(*), pointer :: field
    integer i

    do i = 1, this%fields%size
      field => this%fields%value_at(i)
      select type (field)
      type is (latlon_field2d_type)
        call field%clear()
      type is (latlon_field3d_type)
        call field%clear()
      type is (latlon_field4d_type)
        call field%clear()
      end select
      deallocate(field)
    end do
    call this%fields%clear()

  end subroutine dstate_clear

  subroutine dstate_final(this)

    type(dstate_type), intent(inout) :: this

    call this%clear()

  end subroutine dstate_final

  subroutine dstate_a2c(this)

    class(dstate_type), intent(inout) :: this

    integer i, j, k

    associate (mesh  => this%u%mesh, &
               u     => this%u     , & ! in
               v     => this%v     , & ! in
               u_lon => this%u_lon , & ! out
               v_lat => this%v_lat )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          u_lon%d(i,j,k) = 0.5_r8 * (u%d(i,j,k) + u%d(i+1,j,k))
        end do
      end do
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          v_lat%d(i,j,k) = 0.5_r8 * (v%d(i,j,k) + v%d(i,j+1,k))
        end do
      end do
    end do
    end associate

  end subroutine dstate_a2c

  subroutine dstate_c2a(this)

    class(dstate_type), intent(inout) :: this

    integer i, j, k

    associate (mesh  => this%u%mesh, &
               u     => this%u     , & ! out
               v     => this%v     , & ! out
               u_lon => this%u_lon , & ! in
               v_lat => this%v_lat )   ! in
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          u%d(i,j,k) = 0.5_r8 * (u_lon%d(i,j,k) + u_lon%d(i-1,j,k))
          v%d(i,j,k) = 0.5_r8 * (v_lat%d(i,j,k) + v_lat%d(i,j-1,k))
        end do
      end do
    end do
    end associate

  end subroutine dstate_c2a

  subroutine dtend_init(this, filter_mesh, filter_halo, mesh, halo)

    class(dtend_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: filter_mesh
    type(latlon_halo_type), intent(in) :: filter_halo(:)
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)

    character(field_name_len     ) name
    character(field_long_name_len) long_name
    character(field_units_len    ) units

    call this%clear()

    this%fields = array(30)

    call add_var(this%fields                                                 , &
      name              ='dudt'                                              , &
      long_name         ='Dynamic tendency of u'                             , &
      units             ='m s-2'                                             , &
      loc               ='lon'                                               , &
      mesh              =filter_mesh                                         , &
      halo              =filter_halo                                         , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      ptr               =this%du                                             )
    call add_var(this%fields                                                 , &
      name              ='dvdt'                                              , &
      long_name         ='Dynamic tendency of v'                             , &
      units             ='m s-2'                                             , &
      loc               ='lat'                                               , &
      mesh              =filter_mesh                                         , &
      halo              =filter_halo                                         , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      ptr               =this%dv)
    if (baroclinic) then
      call add_var(this%fields                                               , &
        name            ='dptdt'                                             , &
        long_name       ='Dynamic tendency of pt'                            , &
        units           ='K s-1'                                             , &
        loc             ='cell'                                              , &
        mesh            =filter_mesh                                         , &
        halo            =filter_halo                                         , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%dpt                                            )
      call add_var(this%fields                                               , &
        name            ='dmgsdt'                                            , &
        long_name       ='Dynamic tendency of mgs'                           , &
        units           ='Pa s-1'                                            , &
        loc             ='cell'                                              , &
        mesh            =filter_mesh                                         , &
        halo            =filter_halo                                         , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%dmgs                                           )
    end if

    if (nonhydrostatic .or. .not. baroclinic) then
      call add_var(this%fields                                               , &
        name            ='dgzdt'                                             , &
        long_name       ='Dynamic tendency of gz'                            , &
        units           ='m2 s-2'                                            , &
        loc             ='cell'                                              , &
        mesh            =filter_mesh                                         , &
        halo            =filter_halo                                         , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%dgz                                            )
    end if

#ifdef OUTPUT_H1_DTEND
    call add_var(this%fields                                                 , &
      name              ='dudt_coriolis'                                     , &
      long_name         ='Dynamic tendency of u due to Coriolis force'       , &
      units             ='m s-2'                                             , &
      loc               ='lon'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      ptr               =this%dudt_coriolis                                  )
    call add_var(this%fields                                                 , &
      name              ='dvdt_coriolis'                                     , &
      long_name         ='Dynamic tendency of v due to Coriolis force'       , &
      units             ='m s-2'                                             , &
      loc               ='lat'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      ptr               =this%dvdt_coriolis                                  )
    call add_var(this%fields                                                 , &
      name              ='dudt_dkedx'                                        , &
      long_name         ='Dynamic tendency of u due to kinetic gradient'     , &
      units             ='m s-2'                                             , &
      loc               ='lon'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      ptr               =this%dudt_dkedx                                     )
    call add_var(this%fields                                                 , &
      name              ='dvdt_dkedy'                                        , &
      long_name         ='Dynamic tendency of v due to kinetic gradient'     , &
      units             ='m s-2'                                             , &
      loc               ='lat'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      ptr               =this%dvdt_dkedy                                     )
    call add_var(this%fields                                                 , &
      name              ='dudt_pgf'                                          , &
      long_name         ='Dynamic tendency of u due to PGF'                  , &
      units             ='m s-2'                                             , &
      loc               ='lon'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      ptr               =this%dudt_pgf                                       )
    call add_var(this%fields                                                 , &
      name              ='dvdt_pgf'                                          , &
      long_name         ='Dynamic tendency of v due to PGF'                  , &
      units             ='m s-2'                                             , &
      loc               ='lat'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      ptr               =this%dvdt_pgf                                       )
    if (baroclinic) then
      call add_var(this%fields                                               , &
        name            ='dudt_wedudeta'                                     , &
        long_name       ='Dynamic tendency of u due to vertical advection'   , &
        units           ='m s-2'                                             , &
        loc             ='lon'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%dudt_wedudeta                                  )
      call add_var(this%fields                                               , &
        name            ='dvdt_wedvdeta'                                     , &
        long_name       ='Dynamic tendency of v due to vertical advection'   , &
        units           ='m s-2'                                             , &
        loc             ='lat'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%dvdt_wedvdeta                                  )
    end if
#endif

  end subroutine dtend_init

  subroutine dtend_reset_flags(this)

    class(dtend_type), intent(inout) :: this

    this%du%d = 0
    this%dv%d = 0

    this%update_u   = .false.
    this%update_v   = .false.
    this%update_gz  = .false.
    this%update_pt  = .false.
    this%update_mgs = .false.

  end subroutine dtend_reset_flags

  subroutine dtend_clear(this)

    class(dtend_type), intent(inout) :: this

    class(*), pointer :: field
    integer i

    do i = 1, this%fields%size
      field => this%fields%value_at(i)
      select type (field)
      type is (latlon_field2d_type)
        call field%clear()
      type is (latlon_field3d_type)
        call field%clear()
      type is (latlon_field4d_type)
        call field%clear()
      end select
      deallocate(field)
    end do
    call this%fields%clear()

  end subroutine dtend_clear

  subroutine dtend_final(this)

    type(dtend_type), intent(inout) :: this

    call this%clear()

  end subroutine dtend_final

  subroutine static_init_stage1(this, filter_mesh, filter_halo, mesh, halo)

    class(static_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: filter_mesh
    type(latlon_halo_type), intent(in) :: filter_halo(:)
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)

    call this%clear()

    this%fields = array(30)

    call add_var(this%fields                                                 , &
      name              ='gzs'                                               , &
      long_name         ='Surface geopotential'                              , &
      units             ='m2 s-2'                                            , &
      loc               ='cell'                                              , &
      mesh              =filter_mesh                                         , &
      halo              =filter_halo                                         , &
      output            ='h0'                                                , &
      restart           =.true.                                              , &
      ptr               =this%gzs                                            )
    if (.not. advection) then
      call add_var(this%fields                                               , &
        name            ='landmask'                                          , &
        long_name       ='Land mask'                                         , &
        units           ='1'                                                 , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.true.                                              , &
        ptr             =this%landmask                                       )
      call add_var(this%fields                                               , &
        name            ='zs_std'                                            , &
        long_name       ='Subgrid variance of zs'                            , &
        units           ='m2 s-2'                                            , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.true.                                              , &
        ptr             =this%zs_std                                         )
      call add_var(this%fields                                               , &
        name            ='dzsdx'                                             , &
        long_name       ='Zonal gradient of zs'                              , &
        units           ='1'                                                 , &
        loc             ='lon'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.true.                                              , &
        ptr             =this%dzsdx)
      call add_var(this%fields                                               , &
        name            ='dzsdy'                                             , &
        long_name       ='Meridional gradient of zs'                         , &
        units           ='1'                                                 , &
        loc             ='lat'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.true.                                              , &
        ptr             =this%dzsdy                                          )
    end if
    call add_var(this%fields                                                 , &
      name              ='ref_ps'                                            , &
      long_name         ='Reference surface pressure'                        , &
      units             ='Pa'                                                , &
      loc               ='cell'                                              , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.true.                                              , &
      ptr               =this%ref_ps                                         )
    call add_var(this%fields                                                 , &
      name              ='ref_ps_smth'                                       , &
      long_name         ='Smoothed reference surface pressure'               , &
      units             ='Pa'                                                , &
      loc               ='cell'                                              , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.true.                                              , &
      ptr               =this%ref_ps_smth                                    )
    call add_var(this%fields                                                 , &
      name              ='ref_ps_perb'                                       , &
      long_name         ='Perturbation of reference surface pressure'        , &
      units             ='Pa'                                                , &
      loc               ='cell'                                              , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.true.                                              , &
      ptr               =this%ref_ps_perb                                    )

  end subroutine static_init_stage1

  subroutine static_init_stage2(this, mesh)

    class(static_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: mesh

  end subroutine static_init_stage2

  subroutine static_clear(this)

    class(static_type), intent(inout) :: this

    class(*), pointer :: field
    integer i

    do i = 1, this%fields%size
      field => this%fields%value_at(i)
      select type (field)
      type is (latlon_field2d_type)
        call field%clear()
      type is (latlon_field3d_type)
        call field%clear()
      type is (latlon_field4d_type)
        call field%clear()
      end select
      deallocate(field)
    end do
    call this%fields%clear()

  end subroutine static_clear

  subroutine static_final(this)

    type(static_type), intent(inout) :: this

    call this%clear()

  end subroutine static_final

  subroutine aux_array_init(this, filter_mesh, filter_halo, mesh, halo)

    class(aux_array_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: filter_mesh
    type(latlon_halo_type), intent(in) :: filter_halo(:)
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)

    call this%clear()

    this%fields = array(50)

    if (use_smag_damp) then
      call add_var(this%fields                                               , &
        name            ='smag_t'                                            , &
        long_name       ='Tension of horizontal wind for Smagorinsky damping', &
        units           ='s-2'                                               , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%smag_t                                         )
      call add_var(this%fields                                               , &
        name            ='smag_s'                                            , &
        long_name       ='Shear of horizontal wind for Smagorinsky damping'  , &
        units           ='s-2'                                               , &
        loc             ='vtx'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%smag_s                                         )
      call add_var(this%fields                                               , &
        name            ='kmh'                                               , &
        long_name       ='Horizontal eddy viscosity for Smagorinsky damping' , &
        units           ='s-1'                                               , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%kmh                                            )
      call add_var(this%fields                                               , &
        name            ='kmh_lon'                                           , &
        long_name       ='Horizontal eddy viscosity for Smagorinsky damping' , &
        units           ='s-1'                                               , &
        loc             ='lon'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%kmh_lon                                        )
      call add_var(this%fields                                               , &
        name            ='kmh_lat'                                           , &
        long_name       ='Horizontal eddy viscosity for Smagorinsky damping' , &
        units           ='s-1'                                               , &
        loc             ='lat'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%kmh_lat                                        )
    end if
    if (.not. advection) then
      call add_var(this%fields                                               , &
        name            ='v_lon'                                             , &
        long_name       ='Zonal wind component'                              , &
        units           ='m s-1'                                             , &
        loc             ='lon'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%v_lon                                          )
      call add_var(this%fields                                               , &
        name            ='u_lat'                                             , &
        long_name       ='Meridional wind component'                         , &
        units           ='m s-1'                                             , &
        loc             ='lat'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%u_lat                                          )
      call add_var(this%fields                                               , &
        name            ='ke'                                                , &
        long_name       ='Kinetic energy'                                    , &
        units           ='m2 s-2'                                            , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%ke                                             )
      call add_var(this%fields                                               , &
        name            ='pv_lon'                                            , &
        long_name       ='Potential vorticity'                               , &
        units           ='Pa-1 s-1'                                          , &
        loc             ='lon'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%pv_lon                                         )
      call add_var(this%fields                                               , &
        name            ='pv_lat'                                            , &
        long_name       ='Potential vorticity'                               , &
        units           ='Pa-1 s-1'                                          , &
        loc             ='lat'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%pv_lat                                         )
    end if
    call add_var(this%fields                                                 , &
      name              ='dmg_lon'                                           , &
      long_name         ='Dry-air weight'                                    , &
      units             ='Pa'                                                , &
      loc               ='lon'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      ptr               =this%dmg_lon                                        )
    call add_var(this%fields                                                 , &
      name              ='dmg_lat'                                           , &
      long_name         ='Dry-air weight'                                    , &
      units             ='Pa'                                                , &
      loc               ='lat'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      ptr               =this%dmg_lat                                        )
    call add_var(this%fields                                                 , &
      name              ='dmg_vtx'                                           , &
      long_name         ='Dry-air weight'                                    , &
      units             ='Pa'                                                , &
      loc               ='vtx'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      ptr               =this%dmg_vtx                                        )
    if (baroclinic) then
      call add_var(this%fields                                               , &
        name            ='pkh_lev'                                           , &
        long_name       ='Hydrostatic pressure under Kappa exponent'         , &
        units           ='Pa'                                                , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%pkh_lev                                        )
      call add_var(this%fields                                               , &
        name            ='we_lev_lon'                                        , &
        long_name       ='Vertical mass flux'                                , &
        units           ='Pa s-1'                                            , &
        loc             ='lev_lon'                                           , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%we_lev_lon                                     )
      call add_var(this%fields                                               , &
        name            ='we_lev_lat'                                        , &
        long_name       ='Vertical mass flux'                                , &
        units           ='Pa s-1'                                            , &
        loc             ='lev_lat'                                           , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%we_lev_lat                                     )
      call add_var(this%fields                                               , &
        name            ='ptfx'                                              , &
        long_name       ='Zonal flux of pt'                                  , &
        units           ='K m s-1'                                           , &
        loc             ='lon'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%ptfx                                           )
      call add_var(this%fields                                               , &
        name            ='ptfy'                                              , &
        long_name       ='Meridional flux of pt'                             , &
        units           ='K m s-1'                                           , &
        loc             ='lat'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%ptfy                                           )
      call add_var(this%fields                                               , &
        name            ='ptfz'                                              , &
        long_name       ='Vertical flux of pt'                               , &
        units           ='K m s-1'                                           , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%ptfz                                           )
    end if
    call add_var(this%fields                                                 , &
      name              ='mfx_lon'                                           , &
      long_name         ='Zonal mass flux'                                   , &
      units             ='Pa m s-1'                                          , &
      loc               ='lon'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      ptr               =this%mfx_lon                                        )
    call add_var(this%fields                                                 , &
      name              ='mfy_lat'                                           , &
      long_name         ='Meridional mass flux'                              , &
      units             ='Pa m s-1'                                          , &
      loc               ='lat'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      ptr               =this%mfy_lat)
    call add_var(this%fields                                                 , &
      name              ='mfx_lat'                                           , &
      long_name         ='Zonal mass flux'                                   , &
      units             ='Pa m s-1'                                          , &
      loc               ='lat'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      ptr               =this%mfx_lat                                        )
    call add_var(this%fields                                                 , &
      name              ='mfy_lon'                                           , &
      long_name         ='Meridional mass flux'                              , &
      units             ='Pa m s-1'                                          , &
      loc               ='lon'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      ptr               =this%mfy_lon                                        )
    if (.not. advection) then
      call add_var(this%fields                                               , &
        name            ='vor'                                               , &
        long_name       ='Relative vorticity'                                , &
        units           ='s-1'                                               , &
        loc             ='vtx'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.false.                                             , &
        ptr             =this%vor                                            )
    end if
    call add_var(this%fields                                                 , &
      name              ='pv'                                                , &
      long_name         ='Potential vorticity'                               , &
      units             ='Pa-1 s-1'                                          , &
      loc               ='vtx'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      ptr               =this%pv                                             , &
      halo_cross_pole   =.true.                                              )
    if (.not. advection) then
      call add_var(this%fields                                               , &
        name            ='div'                                               , &
        long_name       ='Divergence'                                        , &
        units           ='s-1'                                               , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.false.                                             , &
        ptr             =this%div                                            )
    end if
    if (use_div_damp .and. div_damp_order == 4) then
      call add_var(this%fields                                               , &
        name            ='div2'                                              , &
        long_name       ='Gradient of divergence'                            , &
        units           ='m-1 s-1'                                           , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.false.                                             , &
        ptr             =this%div2                                           )
    end if
    call add_var(this%fields                                                 , &
      name              ='dmf'                                               , &
      long_name         ='Mass flux divergence'                              , &
      units             ='Pa s-1'                                            , &
      loc               ='cell'                                              , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      ptr               =this%dmf                                            )
    if (nonhydrostatic) then
      call add_var(this%fields                                               , &
        name            ='dmf_lev'                                           , &
        long_name       ='Mass flux divergence'                              , &
        units           ='Pa s-1'                                            , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%dmf_lev                                        )
    end if
    if (baroclinic) then
      call add_var(this%fields                                               , &
        name            ='omg'                                               , &
        long_name       ='Omega'                                             , &
        units           ='Pa s-1'                                            , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%omg                                            )
    end if
    if (pgf_scheme == 'ptb') then
      call add_var(this%fields                                               , &
        name            ='p_ptb'                                             , &
        long_name       ='Perturbation of pressure'                          , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%p_ptb                                          )
      call add_var(this%fields                                               , &
        name            ='gz_ptb'                                            , &
        long_name       ='Perturbation of geopotential'                      , &
        units           ='m2 s-2'                                            , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%gz_ptb                                         )
      call add_var(this%fields                                               , &
        name            ='dp_ptb'                                            , &
        long_name       ='Perturbation of dry-air weight'                    , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%dp_ptb                                         )
      call add_var(this%fields                                               , &
        name            ='ad_ptb'                                            , &
        long_name       ='Perturbation of specific density of dry-air'       , &
        units           ='kg-1 m3'                                           , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%ad_ptb                                         )
    end if
    if (nonhydrostatic) then
      call add_var(this%fields                                               , &
        name            ='u_lev_lon'                                         , &
        long_name       ='Zonal wind component'                              , &
        units           ='m s-1'                                             , &
        loc             ='lev_lon'                                           , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%u_lev_lon                                      )
      call add_var(this%fields                                               , &
        name            ='v_lev_lat'                                         , &
        long_name       ='Meridional wind component'                         , &
        units           ='m s-1'                                             , &
        loc             ='lev_lat'                                           , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%v_lev_lat                                      )
      call add_var(this%fields                                               , &
        name            ='mfx_lev_lon'                                       , &
        long_name       ='Zonal mass flux'                                   , &
        units           ='Pa m s-1'                                          , &
        loc             ='lev_lon'                                           , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%mfx_lev_lon                                    )
      call add_var(this%fields                                               , &
        name            ='mfy_lev_lat'                                       , &
        long_name       ='Meridional mass flux'                              , &
        units           ='Pa m s-1'                                          , &
        loc             ='lev_lat'                                           , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        ptr             =this%mfy_lev_lat                                    )
      call add_var(this%fields                                               , &
        name            ='adv_w_lev'                                         , &
        long_name       ='Advection tendency of w_lev'                       , &
        units           ='m s-2'                                             , &
        loc             ='lev'                                               , &
        mesh            =filter_mesh                                         , &
        halo            =filter_halo                                         , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%adv_w_lev                                      )
      call add_var(this%fields                                               , &
        name            ='adv_gz_lev'                                        , &
        long_name       ='Advection tendency of gz_lev'                      , &
        units           ='m2 s-2'                                            , &
        loc             ='lev'                                               , &
        mesh            =filter_mesh                                         , &
        halo            =filter_halo                                         , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%adv_gz_lev                                     )
    end if

  end subroutine aux_array_init

  subroutine aux_array_init_phys(this, filter_mesh, filter_halo, mesh, halo)

    class(aux_array_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: filter_mesh
    type(latlon_halo_type), intent(in) :: filter_halo(:)
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)

    if (physics_suite /= 'N/A' .and. physics_suite /= '') then
      call add_var(this%fields                                               , &
        name            ='dudt_phys'                                         , &
        long_name       ='Physics tendency of u'                             , &
        units           ='m s-2'                                             , &
        loc             ='lon'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%dudt_phys                                      )
      call add_var(this%fields                                               , &
        name            ='dvdt_phys'                                         , &
        long_name       ='Physics tendency of v'                             , &
        units           ='m s-2'                                             , &
        loc             ='lat'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%dvdt_phys                                      )
      call add_var(this%fields                                               , &
        name            ='dptdt_phys'                                        , &
        long_name       ='Physics tendency of pt'                            , &
        units           ='K s-1'                                             , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%dptdt_phys                                      )
      call add_var(this%fields                                               , &
        name            ='dqdt_phys'                                         , &
        long_name       ='Physics tendency of q'                             , &
        units           ='kg kg-1 s-1'                                       , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        dim4_name       ='tracers'                                           , &
        dim4_size       =ntracers                                            , &
        var4_names      =tracer_names                                        , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        ptr             =this%dqdt_phys                                      )
    end if

  end subroutine aux_array_init_phys

  subroutine aux_array_clear(this)

    class(aux_array_type), intent(inout) :: this

    class(*), pointer :: field
    integer i

    do i = 1, this%fields%size
      field => this%fields%value_at(i)
      select type (field)
      type is (latlon_field2d_type)
        call field%clear()
      type is (latlon_field3d_type)
        call field%clear()
      type is (latlon_field4d_type)
        call field%clear()
      end select
      deallocate(field)
    end do
    call this%fields%clear()

  end subroutine aux_array_clear

  subroutine aux_array_final(this)

    type(aux_array_type), intent(inout) :: this

    call this%clear()

  end subroutine aux_array_final

end module dynamics_types_mod
