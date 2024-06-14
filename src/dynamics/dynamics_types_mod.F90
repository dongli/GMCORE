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
    type(latlon_field3d_type) u
    type(latlon_field3d_type) v
    type(latlon_field3d_type) u_lon
    type(latlon_field3d_type) v_lat
    type(latlon_field3d_type) we_lev
    type(latlon_field3d_type) gz
    type(latlon_field3d_type) gz_lev
    type(latlon_field3d_type) dmg
    type(latlon_field3d_type) dmg_lev
    type(latlon_field3d_type) pt
    type(latlon_field3d_type) t
    type(latlon_field3d_type) tv
    type(latlon_field3d_type) mg
    type(latlon_field3d_type) mg_lev
    type(latlon_field2d_type) mgs
    type(latlon_field3d_type) ph
    type(latlon_field3d_type) ph_lev
    type(latlon_field2d_type) phs
    type(latlon_field3d_type) rhod
    ! Nonhydrostatic variable
    type(latlon_field3d_type) we
    type(latlon_field3d_type) w
    type(latlon_field3d_type) w_lev
    type(latlon_field3d_type) p
    type(latlon_field3d_type) p_lev
    type(latlon_field2d_type) ps
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
    type(latlon_field3d_type) du
    type(latlon_field3d_type) dv
    type(latlon_field3d_type) dgz
    type(latlon_field3d_type) dpt
    type(latlon_field2d_type) dmgs
#ifdef OUTPUT_H1_DTEND
    type(latlon_field3d_type) dudt_coriolis
    type(latlon_field3d_type) dvdt_coriolis
    type(latlon_field3d_type) dudt_wedudeta
    type(latlon_field3d_type) dvdt_wedvdeta
    type(latlon_field3d_type) dudt_dkedx
    type(latlon_field3d_type) dvdt_dkedy
    type(latlon_field3d_type) dudt_pgf
    type(latlon_field3d_type) dvdt_pgf
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
    type(latlon_field2d_type) landmask
    ! Topography
    type(latlon_field2d_type) gzs
    type(latlon_field2d_type) zs_std
    type(latlon_field2d_type) dzsdx
    type(latlon_field2d_type) dzsdy
    ! Reference surface pressure
    type(latlon_field2d_type) ref_ps
    type(latlon_field2d_type) ref_ps_smth
    type(latlon_field2d_type) ref_ps_perb
  contains
    procedure :: init_stage1 => static_init_stage1
    procedure :: init_stage2 => static_init_stage2
    procedure :: clear       => static_clear
    final static_final
  end type static_type

  type aux_array_type
    type(array_type) fields
    ! Smagorinsky damping variables
    type(latlon_field3d_type) smag_t      ! tension strain
    type(latlon_field3d_type) smag_s      ! shear strain on vertex
    type(latlon_field3d_type) kmh         ! nonlinear diffusion coef
    type(latlon_field3d_type) kmh_lon     ! nonlinear diffusion coef on zonal edge
    type(latlon_field3d_type) kmh_lat     ! nonlinear diffusion coef on meridional edge
    ! Other variables
    type(latlon_field3d_type) v_lon       ! Meridional wind speed at lon edge (m s-1)
    type(latlon_field3d_type) u_lat       ! Zonal wind speed at lat edge (m s-1)
    type(latlon_field3d_type) ke          ! Kinetic energy
    type(latlon_field3d_type) pv_lon      ! Potential vorticity on zonal edge
    type(latlon_field3d_type) pv_lat      ! Potential vorticity on merdional edge
    type(latlon_field3d_type) dmg_lon     ! Mass on zonal edge
    type(latlon_field3d_type) dmg_lat     ! Mass on merdional edge
    type(latlon_field3d_type) dmg_vtx     ! Mass on vertex
    type(latlon_field3d_type) pkh_lev     ! Exner pressure on half levels
    type(latlon_field3d_type) we_lev_lon  ! Vertical coordinate speed multiplied by 𝛛π/𝛛η on zonal edge
    type(latlon_field3d_type) we_lev_lat  ! Vertical coordinate speed multiplied by 𝛛π/𝛛η on merdional edge
    type(latlon_field3d_type) ptfx        ! Potential temperature on the zonal edge
    type(latlon_field3d_type) ptfy        ! Potential temperature on the merdional edge
    type(latlon_field3d_type) ptfz        ! Potential temperature on the vertical edge
    type(latlon_field3d_type) mfx_lon     ! Normal mass flux on zonal edge
    type(latlon_field3d_type) mfy_lat     ! Normal mass flux on merdional edge
    type(latlon_field3d_type) mfx_lat     ! Tangient mass flux on zonal edge
    type(latlon_field3d_type) mfy_lon     ! Tangient mass flux on merdional edge
    type(latlon_field3d_type) vor         ! Vorticity (s-1)
    type(latlon_field3d_type) pv          ! Potential vorticity
    type(latlon_field3d_type) div         ! Divergence (s-1)
    type(latlon_field3d_type) div2        ! Laplacian of divergence (s-1)
    type(latlon_field3d_type) dmf         ! Mass flux divergence on full level (Pa s-1)
    type(latlon_field3d_type) dmf_lev     ! Mass flux divergence on half level (Pa s-1)
    type(latlon_field3d_type) omg         ! Vertical pressure velocity (Pa s-1)
    ! Tendencies from physics
    type(latlon_field3d_type) dudt_phys
    type(latlon_field3d_type) dvdt_phys
    type(latlon_field3d_type) dptdt_phys
    type(latlon_field4d_type) dqdt_phys
    ! Perturbed quantities for calculating HPGF
    type(latlon_field3d_type) p_ptb
    type(latlon_field3d_type) gz_ptb
    type(latlon_field3d_type) dp_ptb
    type(latlon_field3d_type) ad_ptb
    ! Nonhydrostatic variables
    type(latlon_field3d_type) u_lev_lon
    type(latlon_field3d_type) v_lev_lat
    type(latlon_field3d_type) mfx_lev_lon
    type(latlon_field3d_type) mfy_lev_lat
    type(latlon_field3d_type) adv_w_lev
    type(latlon_field3d_type) adv_gz_lev
  contains
    procedure :: init      => aux_array_init
    procedure :: init_phys => aux_array_init_phys
    procedure :: clear     => aux_array_clear
    final aux_array_final
  end type aux_array_type

contains

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
      call append_field(this%fields                                          , &
        name            ='u'                                                 , &
        long_name       ='Zonal wind component'                              , &
        units           ='m s-1'                                             , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.false.                                             , &
        field           =this%u                                              )
    end if
    if (.not. advection) then
      call append_field(this%fields                                          , &
        name            ='v'                                                 , &
        long_name       ='Meridional wind component'                         , &
        units           ='m s-1'                                             , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.false.                                             , &
        field           =this%v                                              )
    end if
    call append_field(this%fields                                            , &
      name              ='u_lon'                                             , &
      long_name         ='Zonal wind component'                              , &
      units             ='m s-1'                                             , &
      loc               ='lon'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.true.                                              , &
      field             =this%u_lon                                          )
    call append_field(this%fields                                            , &
      name              ='v_lat'                                             , &
      long_name         ='Meridional wind component'                         , &
      units             ='m s-1'                                             , &
      loc               ='lat'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.true.                                              , &
      field             =this%v_lat                                          )
    if (baroclinic .or. advection) then
      call append_field(this%fields                                          , &
        name            ='we_lev'                                            , &
        long_name       ='Vertical mass flux'                                , &
        units           ='Pa s-1'                                            , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%we_lev                                         )
    end if
    call append_field(this%fields                                            , &
      name              ='gz'                                                , &
      long_name         ='Geopotential'                                      , &
      units             ='m2 s-2'                                            , &
      loc               ='cell'                                              , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h0'                                                , &
      restart           =.not. baroclinic                                    , & ! FIXME: Revise this.
      field             =this%gz                                             )
    if (nonhydrostatic) then
      call append_field(this%fields                                          , &
        name            ='gz_lev'                                            , &
        long_name       ='Geopotential'                                      , &
        units           ='m2 s-2'                                            , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.true.                                              , &
        field           =this%gz_lev                                         , &
        halo_cross_pole =.true.                                              )
    else
      call append_field(this%fields                                          , &
        name            ='gz_lev'                                            , &
        long_name       ='Geopotential'                                      , &
        units           ='m2 s-2'                                            , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%gz_lev                                         )
    end if
    call append_field(this%fields                                            , &
      name              ='dmg'                                               , &
      long_name         ='Dry-air weight'                                    , &
      units             ='Pa'                                                , &
      loc               ='cell'                                              , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      field             =this%dmg                                            )
    if (baroclinic .or. advection) then
      call append_field(this%fields                                          , &
        name            ='dmg_lev'                                           , &
        long_name       ='Dry-air weight'                                    , &
        units           ='Pa'                                                , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%dmg_lev                                        )
    end if
    if (baroclinic .or. advection) then
      call append_field(this%fields                                          , &
        name            ='pt'                                                , &
        long_name       ='Modified potential temperature'                    , &
        units           ='K'                                                 , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.true.                                              , &
        field           =this%pt                                             , &
        halo_cross_pole =.true.                                              )
    end if
    if (baroclinic .or. advection) then
      call append_field(this%fields                                          , &
        name            ='t'                                                 , &
        long_name       ='Temperature'                                       , &
        units           ='K'                                                 , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.false.                                             , &
        field           =this%t                                              )
    end if
    if (baroclinic .or. advection) then
      call append_field(this%fields                                          , &
        name            ='tv'                                                , &
        long_name       ='Virtual temperature'                               , &
        units           ='K'                                                 , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%tv                                             )
    end if
    call append_field(this%fields                                            , &
      name              ='mg'                                                , &
      long_name         ='Dry-air weight'                                    , &
      units             ='Pa'                                                , &
      loc               ='cell'                                              , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      field             =this%mg                                             )
    if (baroclinic .or. advection) then
      call append_field(this%fields                                          , &
        name            ='mg_lev'                                            , &
        long_name       ='Dry-air weight'                                    , &
        units           ='Pa'                                                , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%mg_lev                                         )
    end if
    if (baroclinic .or. advection) then
      call append_field(this%fields                                          , &
        name            ='mgs'                                               , &
        long_name       ='Dry-air weight on surface'                         , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.true.                                              , &
        field           =this%mgs                                            )
    end if
    if (baroclinic .or. advection) then
      call append_field(this%fields                                          , &
        name            ='ph'                                                , &
        long_name       ='Hydrostatic pressure'                              , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%ph                                             )
    end if
    if (baroclinic .or. advection) then
      call append_field(this%fields                                          , &
        name            ='ph_lev'                                            , &
        long_name       ='Hydrostatic pressure'                              , &
        units           ='Pa'                                                , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%ph_lev                                         )
    end if
    if (baroclinic) then
      call append_field(this%fields                                          , &
        name            ='phs'                                               , &
        long_name       ='Hydrostatic pressure on surface'                   , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.false.                                             , &
        field           =this%phs                                            )
      call this%phs%link(this%ph_lev, mesh%half_nlev)
    end if
    if (baroclinic) then
      call append_field(this%fields                                          , &
        name            ='rhod'                                              , &
        long_name       ='Dry-air density'                                   , &
        units           ='kg m-3'                                            , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.false.                                             , &
        field           =this%rhod                                           )
    end if
    if (nonhydrostatic) then
      call append_field(this%fields                                          , &
        name            ='we'                                                , &
        long_name       ='Vertical mass flux'                                , &
        units           ='Pa s-1'                                            , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%we                                             )
    end if
    if (nonhydrostatic) then
      call append_field(this%fields                                          , &
        name            ='w'                                                 , &
        long_name       ='Vertical wind speed'                               , &
        units           ='m s-1'                                             , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%w                                              )
    end if
    if (nonhydrostatic) then
      call append_field(this%fields                                          , &
        name            ='w_lev'                                             , &
        long_name       ='Vertical wind speed'                               , &
        units           ='m s-1'                                             , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.true.                                              , &
        field           =this%w_lev                                          , &
        halo_cross_pole =.true.                                              )
    end if
    if (hydrostatic) then
      call append_field(this%fields                                          , &
        name            ='p'                                                 , &
        long_name       ='Pressure'                                          , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%p                                              , &
        link            =this%ph                                             )
    else if (nonhydrostatic) then
      call append_field(this%fields                                          , &
        name            ='p'                                                 , &
        long_name       ='Pressure'                                          , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%p                                              )
    end if
    if (hydrostatic) then
      call append_field(this%fields                                          , &
        name            ='p_lev'                                             , &
        long_name       ='Pressure'                                          , &
        units           ='Pa'                                                , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%p_lev                                          , &
        link            =this%ph_lev                                         )
    else if (nonhydrostatic) then
      call append_field(this%fields                                          , &
        name            ='p_lev'                                             , &
        long_name       ='Pressure'                                          , &
        units           ='Pa'                                                , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%p_lev                                          )
    end if
    if (hydrostatic) then
      call append_field(this%fields                                          , &
        name            ='ps'                                                , &
        long_name       ='Surface pressure'                                  , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%ps                                             , &
        link            =this%phs                                            )
    else if (nonhydrostatic) then
      call append_field(this%fields                                          , &
        name            ='ps'                                                , &
        long_name       ='Surface pressure'                                  , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%ps                                             )
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

    call append_field(this%fields                                            , &
      name              ='dudt'                                              , &
      long_name         ='Dynamic tendency of u'                             , &
      units             ='m s-2'                                             , &
      loc               ='lon'                                               , &
      mesh              =filter_mesh                                         , &
      halo              =filter_halo                                         , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      field             =this%du                                             )
    call append_field(this%fields                                            , &
      name              ='dvdt'                                              , &
      long_name         ='Dynamic tendency of v'                             , &
      units             ='m s-2'                                             , &
      loc               ='lat'                                               , &
      mesh              =filter_mesh                                         , &
      halo              =filter_halo                                         , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      field             =this%dv)
    if (baroclinic) then
      call append_field(this%fields                                          , &
        name            ='dptdt'                                             , &
        long_name       ='Dynamic tendency of pt'                            , &
        units           ='K s-1'                                             , &
        loc             ='cell'                                              , &
        mesh            =filter_mesh                                         , &
        halo            =filter_halo                                         , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%dpt                                            )
      call append_field(this%fields                                          , &
        name            ='dmgsdt'                                            , &
        long_name       ='Dynamic tendency of mgs'                           , &
        units           ='Pa s-1'                                            , &
        loc             ='cell'                                              , &
        mesh            =filter_mesh                                         , &
        halo            =filter_halo                                         , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%dmgs                                           )
    end if

    if (nonhydrostatic .or. .not. baroclinic) then
      call append_field(this%fields                                          , &
        name            ='dgzdt'                                             , &
        long_name       ='Dynamic tendency of gz'                            , &
        units           ='m2 s-2'                                            , &
        loc             ='cell'                                              , &
        mesh            =filter_mesh                                         , &
        halo            =filter_halo                                         , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%dgz                                            )
    end if

#ifdef OUTPUT_H1_DTEND
    call append_field(this%fields                                            , &
      name              ='dudt_coriolis'                                     , &
      long_name         ='Dynamic tendency of u due to Coriolis force'       , &
      units             ='m s-2'                                             , &
      loc               ='lon'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      field             =this%dudt_coriolis                                  )
    call append_field(this%fields                                            , &
      name              ='dvdt_coriolis'                                     , &
      long_name         ='Dynamic tendency of v due to Coriolis force'       , &
      units             ='m s-2'                                             , &
      loc               ='lat'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      field             =this%dvdt_coriolis                                  )
    call append_field(this%fields                                            , &
      name              ='dudt_dkedx'                                        , &
      long_name         ='Dynamic tendency of u due to kinetic gradient'     , &
      units             ='m s-2'                                             , &
      loc               ='lon'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      field             =this%dudt_dkedx                                     )
    call append_field(this%fields                                            , &
      name              ='dvdt_dkedy'                                        , &
      long_name         ='Dynamic tendency of v due to kinetic gradient'     , &
      units             ='m s-2'                                             , &
      loc               ='lat'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      field             =this%dvdt_dkedy                                     )
    call append_field(this%fields                                            , &
      name              ='dudt_pgf'                                          , &
      long_name         ='Dynamic tendency of u due to PGF'                  , &
      units             ='m s-2'                                             , &
      loc               ='lon'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      field             =this%dudt_pgf                                       )
    call append_field(this%fields                                            , &
      name              ='dvdt_pgf'                                          , &
      long_name         ='Dynamic tendency of v due to PGF'                  , &
      units             ='m s-2'                                             , &
      loc               ='lat'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      field             =this%dvdt_pgf                                       )
    if (baroclinic) then
      call append_field(this%fields                                          , &
        name            ='dudt_wedudeta'                                     , &
        long_name       ='Dynamic tendency of u due to vertical advection'   , &
        units           ='m s-2'                                             , &
        loc             ='lon'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%dudt_wedudeta                                  )
      call append_field(this%fields                                          , &
        name            ='dvdt_wedvdeta'                                     , &
        long_name       ='Dynamic tendency of v due to vertical advection'   , &
        units           ='m s-2'                                             , &
        loc             ='lat'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%dvdt_wedvdeta                                  )
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

    call append_field(this%fields                                            , &
      name              ='gzs'                                               , &
      long_name         ='Surface geopotential'                              , &
      units             ='m2 s-2'                                            , &
      loc               ='cell'                                              , &
      mesh              =filter_mesh                                         , &
      halo              =filter_halo                                         , &
      output            ='h0'                                                , &
      restart           =.true.                                              , &
      field             =this%gzs                                            )
    if (.not. advection) then
      call append_field(this%fields                                          , &
        name            ='landmask'                                          , &
        long_name       ='Land mask'                                         , &
        units           ='1'                                                 , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.true.                                              , &
        field           =this%landmask                                       )
      call append_field(this%fields                                          , &
        name            ='zs_std'                                            , &
        long_name       ='Subgrid variance of zs'                            , &
        units           ='m2 s-2'                                            , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.true.                                              , &
        field           =this%zs_std                                         )
      call append_field(this%fields                                          , &
        name            ='dzsdx'                                             , &
        long_name       ='Zonal gradient of zs'                              , &
        units           ='1'                                                 , &
        loc             ='lon'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.true.                                              , &
        field           =this%dzsdx)
      call append_field(this%fields                                          , &
        name            ='dzsdy'                                             , &
        long_name       ='Meridional gradient of zs'                         , &
        units           ='1'                                                 , &
        loc             ='lat'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.true.                                              , &
        field           =this%dzsdy                                          )
    end if
    call append_field(this%fields                                            , &
      name              ='ref_ps'                                            , &
      long_name         ='Reference surface pressure'                        , &
      units             ='Pa'                                                , &
      loc               ='cell'                                              , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.true.                                              , &
      field             =this%ref_ps                                         )
    call append_field(this%fields                                            , &
      name              ='ref_ps_smth'                                       , &
      long_name         ='Smoothed reference surface pressure'               , &
      units             ='Pa'                                                , &
      loc               ='cell'                                              , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.true.                                              , &
      field             =this%ref_ps_smth                                    )
    call append_field(this%fields                                            , &
      name              ='ref_ps_perb'                                       , &
      long_name         ='Perturbation of reference surface pressure'        , &
      units             ='Pa'                                                , &
      loc               ='cell'                                              , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.true.                                              , &
      field             =this%ref_ps_perb                                    )

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
      call append_field(this%fields                                          , &
        name            ='smag_t'                                            , &
        long_name       ='Tension of horizontal wind for Smagorinsky damping', &
        units           ='s-2'                                               , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%smag_t                                         )
      call append_field(this%fields                                          , &
        name            ='smag_s'                                            , &
        long_name       ='Shear of horizontal wind for Smagorinsky damping'  , &
        units           ='s-2'                                               , &
        loc             ='vtx'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%smag_s                                         )
      call append_field(this%fields                                          , &
        name            ='kmh'                                               , &
        long_name       ='Horizontal eddy viscosity for Smagorinsky damping' , &
        units           ='s-1'                                               , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%kmh                                            )
      call append_field(this%fields                                          , &
        name            ='kmh_lon'                                           , &
        long_name       ='Horizontal eddy viscosity for Smagorinsky damping' , &
        units           ='s-1'                                               , &
        loc             ='lon'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%kmh_lon                                        )
      call append_field(this%fields                                          , &
        name            ='kmh_lat'                                           , &
        long_name       ='Horizontal eddy viscosity for Smagorinsky damping' , &
        units           ='s-1'                                               , &
        loc             ='lat'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%kmh_lat                                        )
    end if
    if (.not. advection) then
      call append_field(this%fields                                          , &
        name            ='v_lon'                                             , &
        long_name       ='Zonal wind component'                              , &
        units           ='m s-1'                                             , &
        loc             ='lon'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%v_lon                                          )
      call append_field(this%fields                                          , &
        name            ='u_lat'                                             , &
        long_name       ='Meridional wind component'                         , &
        units           ='m s-1'                                             , &
        loc             ='lat'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%u_lat                                          )
      call append_field(this%fields                                          , &
        name            ='ke'                                                , &
        long_name       ='Kinetic energy'                                    , &
        units           ='m2 s-2'                                            , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%ke                                             )
      call append_field(this%fields                                          , &
        name            ='pv_lon'                                            , &
        long_name       ='Potential vorticity'                               , &
        units           ='Pa-1 s-1'                                          , &
        loc             ='lon'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%pv_lon                                         )
      call append_field(this%fields                                          , &
        name            ='pv_lat'                                            , &
        long_name       ='Potential vorticity'                               , &
        units           ='Pa-1 s-1'                                          , &
        loc             ='lat'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%pv_lat                                         )
    end if
    call append_field(this%fields                                            , &
      name              ='dmg_lon'                                           , &
      long_name         ='Dry-air weight'                                    , &
      units             ='Pa'                                                , &
      loc               ='lon'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      field             =this%dmg_lon                                        )
    call append_field(this%fields                                            , &
      name              ='dmg_lat'                                           , &
      long_name         ='Dry-air weight'                                    , &
      units             ='Pa'                                                , &
      loc               ='lat'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      field             =this%dmg_lat                                        )
    call append_field(this%fields                                            , &
      name              ='dmg_vtx'                                           , &
      long_name         ='Dry-air weight'                                    , &
      units             ='Pa'                                                , &
      loc               ='vtx'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      field             =this%dmg_vtx                                        )
    if (baroclinic) then
      call append_field(this%fields                                          , &
        name            ='pkh_lev'                                           , &
        long_name       ='Hydrostatic pressure under Kappa exponent'         , &
        units           ='Pa'                                                , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%pkh_lev                                        )
      call append_field(this%fields                                          , &
        name            ='we_lev_lon'                                        , &
        long_name       ='Vertical mass flux'                                , &
        units           ='Pa s-1'                                            , &
        loc             ='lev_lon'                                           , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%we_lev_lon                                     )
      call append_field(this%fields                                          , &
        name            ='we_lev_lat'                                        , &
        long_name       ='Vertical mass flux'                                , &
        units           ='Pa s-1'                                            , &
        loc             ='lev_lat'                                           , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%we_lev_lat                                     )
      call append_field(this%fields                                          , &
        name            ='ptfx'                                              , &
        long_name       ='Zonal flux of pt'                                  , &
        units           ='K m s-1'                                           , &
        loc             ='lon'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%ptfx                                           )
      call append_field(this%fields                                          , &
        name            ='ptfy'                                              , &
        long_name       ='Meridional flux of pt'                             , &
        units           ='K m s-1'                                           , &
        loc             ='lat'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%ptfy                                           )
      call append_field(this%fields                                          , &
        name            ='ptfz'                                              , &
        long_name       ='Vertical flux of pt'                               , &
        units           ='K m s-1'                                           , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%ptfz                                           )
    end if
    call append_field(this%fields                                            , &
      name              ='mfx_lon'                                           , &
      long_name         ='Zonal mass flux'                                   , &
      units             ='Pa m s-1'                                          , &
      loc               ='lon'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      field             =this%mfx_lon                                        )
    call append_field(this%fields                                            , &
      name              ='mfy_lat'                                           , &
      long_name         ='Meridional mass flux'                              , &
      units             ='Pa m s-1'                                          , &
      loc               ='lat'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      field             =this%mfy_lat)
    call append_field(this%fields                                            , &
      name              ='mfx_lat'                                           , &
      long_name         ='Zonal mass flux'                                   , &
      units             ='Pa m s-1'                                          , &
      loc               ='lat'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      field             =this%mfx_lat                                        )
    call append_field(this%fields                                            , &
      name              ='mfy_lon'                                           , &
      long_name         ='Meridional mass flux'                              , &
      units             ='Pa m s-1'                                          , &
      loc               ='lon'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      field             =this%mfy_lon                                        )
    if (.not. advection) then
      call append_field(this%fields                                          , &
        name            ='vor'                                               , &
        long_name       ='Relative vorticity'                                , &
        units           ='s-1'                                               , &
        loc             ='vtx'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.false.                                             , &
        field           =this%vor                                            )
    end if
    call append_field(this%fields                                            , &
      name              ='pv'                                                , &
      long_name         ='Potential vorticity'                               , &
      units             ='Pa-1 s-1'                                          , &
      loc               ='vtx'                                               , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =''                                                  , &
      restart           =.false.                                             , &
      field             =this%pv                                             , &
      halo_cross_pole   =.true.                                              )
    if (.not. advection) then
      call append_field(this%fields                                          , &
        name            ='div'                                               , &
        long_name       ='Divergence'                                        , &
        units           ='s-1'                                               , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.false.                                             , &
        field           =this%div                                            )
    end if
    if (use_div_damp .and. div_damp_order == 4) then
      call append_field(this%fields                                          , &
        name            ='div2'                                              , &
        long_name       ='Gradient of divergence'                            , &
        units           ='m-1 s-1'                                           , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h0'                                                , &
        restart         =.false.                                             , &
        field           =this%div2                                           )
    end if
    call append_field(this%fields                                            , &
      name              ='dmf'                                               , &
      long_name         ='Mass flux divergence'                              , &
      units             ='Pa s-1'                                            , &
      loc               ='cell'                                              , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      field             =this%dmf                                            )
    if (nonhydrostatic) then
      call append_field(this%fields                                          , &
        name            ='dmf_lev'                                           , &
        long_name       ='Mass flux divergence'                              , &
        units           ='Pa s-1'                                            , &
        loc             ='lev'                                               , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%dmf_lev                                        )
    end if
    if (baroclinic) then
      call append_field(this%fields                                          , &
        name            ='omg'                                               , &
        long_name       ='Omega'                                             , &
        units           ='Pa s-1'                                            , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%omg                                            )
    end if
    if (pgf_scheme == 'ptb') then
      call append_field(this%fields                                          , &
        name            ='p_ptb'                                             , &
        long_name       ='Perturbation of pressure'                          , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%p_ptb                                          )
      call append_field(this%fields                                          , &
        name            ='gz_ptb'                                            , &
        long_name       ='Perturbation of geopotential'                      , &
        units           ='m2 s-2'                                            , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%gz_ptb                                         )
      call append_field(this%fields                                          , &
        name            ='dp_ptb'                                            , &
        long_name       ='Perturbation of dry-air weight'                    , &
        units           ='Pa'                                                , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%dp_ptb                                         )
      call append_field(this%fields                                          , &
        name            ='ad_ptb'                                            , &
        long_name       ='Perturbation of specific density of dry-air'       , &
        units           ='kg-1 m3'                                           , &
        loc             ='cell'                                              , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%ad_ptb                                         )
    end if
    if (nonhydrostatic) then
      call append_field(this%fields                                          , &
        name            ='u_lev_lon'                                         , &
        long_name       ='Zonal wind component'                              , &
        units           ='m s-1'                                             , &
        loc             ='lev_lon'                                           , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%u_lev_lon                                      )
      call append_field(this%fields                                          , &
        name            ='v_lev_lat'                                         , &
        long_name       ='Meridional wind component'                         , &
        units           ='m s-1'                                             , &
        loc             ='lev_lat'                                           , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%v_lev_lat                                      )
      call append_field(this%fields                                          , &
        name            ='mfx_lev_lon'                                       , &
        long_name       ='Zonal mass flux'                                   , &
        units           ='Pa m s-1'                                          , &
        loc             ='lev_lon'                                           , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%mfx_lev_lon                                    )
      call append_field(this%fields                                          , &
        name            ='mfy_lev_lat'                                       , &
        long_name       ='Meridional mass flux'                              , &
        units           ='Pa m s-1'                                          , &
        loc             ='lev_lat'                                           , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =''                                                  , &
        restart         =.false.                                             , &
        field           =this%mfy_lev_lat                                    )
      call append_field(this%fields                                          , &
        name            ='adv_w_lev'                                         , &
        long_name       ='Advection tendency of w_lev'                       , &
        units           ='m s-2'                                             , &
        loc             ='lev'                                               , &
        mesh            =filter_mesh                                         , &
        halo            =filter_halo                                         , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%adv_w_lev                                      )
      call append_field(this%fields                                          , &
        name            ='adv_gz_lev'                                        , &
        long_name       ='Advection tendency of gz_lev'                      , &
        units           ='m2 s-2'                                            , &
        loc             ='lev'                                               , &
        mesh            =filter_mesh                                         , &
        halo            =filter_halo                                         , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%adv_gz_lev                                     )
    end if

  end subroutine aux_array_init

  subroutine aux_array_init_phys(this, filter_mesh, filter_halo, mesh, halo)

    class(aux_array_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in), target :: filter_mesh
    type(latlon_halo_type), intent(in), target :: filter_halo(:)
    type(latlon_mesh_type), intent(in), target :: mesh
    type(latlon_halo_type), intent(in), target :: halo(:)

    type(latlon_mesh_type), pointer :: mesh_ptr
    type(latlon_halo_type), pointer :: halo_ptr(:)

    if (filter_ptend) then
      mesh_ptr => filter_mesh
      halo_ptr => filter_halo
    else
      mesh_ptr => mesh
      halo_ptr => halo
    end if

    if (physics_suite /= 'N/A' .and. physics_suite /= '') then
      call append_field(this%fields                                          , &
        name            ='dudt_phys'                                         , &
        long_name       ='Physics tendency of u'                             , &
        units           ='m s-2'                                             , &
        loc             ='lon'                                               , &
        mesh            =mesh_ptr                                            , &
        halo            =halo_ptr                                            , &
        output          ='h1'                                                , &
        restart         =.true.                                              , &
        field           =this%dudt_phys                                      )
      call append_field(this%fields                                          , &
        name            ='dvdt_phys'                                         , &
        long_name       ='Physics tendency of v'                             , &
        units           ='m s-2'                                             , &
        loc             ='lat'                                               , &
        mesh            =mesh_ptr                                            , &
        halo            =halo_ptr                                            , &
        output          ='h1'                                                , &
        restart         =.true.                                              , &
        field           =this%dvdt_phys                                      )
      call append_field(this%fields                                          , &
        name            ='dptdt_phys'                                        , &
        long_name       ='Physics tendency of pt'                            , &
        units           ='K s-1'                                             , &
        loc             ='cell'                                              , &
        mesh            =mesh_ptr                                            , &
        halo            =halo_ptr                                            , &
        output          ='h1'                                                , &
        restart         =.true.                                              , &
        field           =this%dptdt_phys                                     )
      call append_field(this%fields                                          , &
        name            ='dqdt_phys'                                         , &
        long_name       ='Physics tendency of q'                             , &
        units           ='kg kg-1 s-1'                                       , &
        loc             ='cell'                                              , &
        mesh            =mesh_ptr                                            , &
        dim4_name       ='tracers'                                           , &
        dim4_size       =ntracers                                            , &
        var4_names      =tracer_names                                        , &
        halo            =halo_ptr                                            , &
        output          ='h1'                                                , &
        restart         =.true.                                              , &
        field           =this%dqdt_phys                                      )
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
    end do
    call this%fields%clear()

  end subroutine aux_array_clear

  subroutine aux_array_final(this)

    type(aux_array_type), intent(inout) :: this

    call this%clear()

  end subroutine aux_array_final

end module dynamics_types_mod
