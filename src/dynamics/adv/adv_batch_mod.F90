! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================
! Description:
!
!   The tracer advection can be seperated into different batches. Each batch
!   can have different time step size. The wind and mass flux are accumulated
!   along model integration, and averaged to middle time level of advection time
!   step cycle.
!
!   The batch type allocates necessary arrays, and provides wind accumulation
!   subroutines.
!
! Note:
!
!   - It needs to verify the wind and mass flux accumulation manners:
!     Averaging wind and mass flux on n + 1/2 time level, or n time level.
!
! Authors:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
! ==============================================================================

module adv_batch_mod

  use flogger
  use container
  use const_mod
  use namelist_mod
  use math_mod
  use time_mod
  use latlon_mesh_mod
  use latlon_field_types_mod
  use latlon_operators_mod
  use latlon_parallel_mod
  use filter_mod
  use vert_coord_mod

  implicit none

  private

  public adv_batch_type
  public adv_fill_vhalo

  ! Different tracers can be combined into one batch, and advected in different
  ! time steps.
  type adv_batch_type
    character(30) :: scheme_h = 'N/A'
    character(30) :: scheme_v = 'N/A'
    character(10) :: loc  = 'cell'
    character(30) :: name = ''
    logical  :: initialized = .false.
    logical  :: dynamic     = .false.
    logical  :: passive     = .true.
    integer  :: ntracers    = 1
    integer  :: nstep       = 0     ! Number of dynamic steps for one adv step
    integer  :: step        = 0     ! Step counter
    real(r8) :: dt                  ! Advection time step size in seconds
    integer , allocatable :: idx(:) ! Global index of tracers in this batch
    type(latlon_mesh_type), pointer :: mesh => null()
    type(array_type) fields
    type(latlon_field3d_type) m     ! Dry-air weight at n time level
    type(latlon_field3d_type) mfx   ! Mass flux in x direction
    type(latlon_field3d_type) mfy   ! Mass flux in y direction
    type(latlon_field3d_type) mfz   ! Mass flux in z direction
    type(latlon_field3d_type) u
    type(latlon_field3d_type) v
    type(latlon_field3d_type) w     ! Only for advection tests
    type(latlon_field3d_type) qmfx
    type(latlon_field3d_type) qmfy
    type(latlon_field3d_type) qmfz
    ! FFSL variables
    type(latlon_field3d_type) u_frac
    type(latlon_field3d_type) w_frac
    type(latlon_field3d_type) mfx_frac
    type(latlon_field3d_type) mfz_frac
    type(latlon_field3d_type) cflx
    type(latlon_field3d_type) cfly
    type(latlon_field3d_type) cflz
    type(latlon_field3d_type) divx
    type(latlon_field3d_type) divy
    type(latlon_field3d_type) qmfx0
    type(latlon_field3d_type) qmfy0
    type(latlon_field3d_type) qx
    type(latlon_field3d_type) qy
    type(adv_batch_type), pointer :: bg => null() ! Background batch
  contains
    procedure :: init       => adv_batch_init
    procedure :: clear      => adv_batch_clear
    procedure :: copy_m_old => adv_batch_copy_m_old
    procedure :: set_wind   => adv_batch_set_wind
    procedure :: accum_wind => adv_batch_accum_wind
    procedure :: calc_cflxy_mass   => adv_batch_calc_cflxy_mass
    procedure :: calc_cflxy_tracer => adv_batch_calc_cflxy_tracer
    procedure :: calc_cflz_mass    => adv_batch_calc_cflz_mass
    procedure :: calc_cflz_tracer  => adv_batch_calc_cflz_tracer
    procedure, private :: prepare_ffsl_h => adv_batch_prepare_ffsl_h
    procedure, private :: prepare_ffsl_v => adv_batch_prepare_ffsl_v
    final :: adv_batch_final
  end type adv_batch_type

contains

  subroutine adv_batch_init(this, filter, filter_mesh, filter_halo, mesh, halo, scheme, batch_loc, batch_name, dt, dynamic, passive, idx, bg)

    class(adv_batch_type), intent(inout) :: this
    type(filter_type), intent(in) :: filter
    type(latlon_mesh_type), intent(in), target :: filter_mesh
    type(latlon_halo_type), intent(in) :: filter_halo(:)
    type(latlon_mesh_type), intent(in), target :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)
    character(*), intent(in) :: scheme
    character(*), intent(in) :: batch_loc
    character(*), intent(in) :: batch_name
    real(r8), intent(in) :: dt
    logical, intent(in) :: dynamic
    logical, intent(in) :: passive
    integer, intent(in), optional :: idx(:)
    class(adv_batch_type), intent(in), target, optional :: bg

    call this%clear()

    this%mesh => filter_mesh
    this%fields = array(30)

    this%loc      = batch_loc
    this%name     = batch_name
    this%dt       = dt
    this%dynamic  = dynamic
    this%passive  = passive
    this%nstep    = dt / dt_dyn
    this%step     = 0

    if (present(bg)) this%bg => bg

    ! Discriminate horizontal and vertical schemes.
    if (count_string(scheme, ':') == 1) then
      this%scheme_h = split_string(scheme, ':', 1)
      this%scheme_v = split_string(scheme, ':', 2)
    else
      this%scheme_h = scheme
      this%scheme_v = scheme
    end if

    select case (batch_loc)
    case ('cell')
      if (this%passive) then
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_m'                             , &
          long_name       ='Dry-air weight'                                    , &
          units           ='Pa'                                                , &
          loc             ='cell'                                              , &
          mesh            =filter_mesh                                         , &
          halo            =filter_halo                                         , &
          output          =merge('h0', '  ', advection)                        , &
          restart         =.false.                                             , &
          field           =this%m                                              )
      end if
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_u'                             , &
        long_name         ='U wind component'                                  , &
        units             ='m s-1'                                             , &
        loc               ='lon'                                               , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        output            =merge('h0', '  ', advection)                        , &
        restart           =.false.                                             , &
        field             =this%u                                              )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_v'                             , &
        long_name         ='V wind component'                                  , &
        units             ='m s-1'                                             , &
        loc               ='lat'                                               , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        output            =merge('h0', '  ', advection)                        , &
        restart           =.false.                                             , &
        field             =this%v                                              )
      if (advection .and. nlev > 1) then
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_w'                             , &
          long_name       ='Vertical coordinate velocity'                      , &
          units           ='s-1'                                               , &
          loc             ='lev'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          output          =merge('h0', '  ', advection)                        , &
          restart         =.false.                                             , &
          field           =this%w                                              )
      end if
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_mfz'                           , &
        long_name         ='Vertical mass flux'                                , &
        units             ='Pa s-1'                                            , &
        loc               ='lev'                                               , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        output            =merge('h0', '  ', advection)                        , &
        restart           =.false.                                             , &
        field             =this%mfz                                            )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_mfx'                           , &
        long_name         ='Mass flux in x direction'                          , &
        units             ='Pa m s-1'                                          , &
        loc               ='lon'                                               , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        output            =merge('h0', '  ', advection)                        , &
        restart           =.false.                                             , &
        field             =this%mfx                                            )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_mfy'                           , &
        long_name         ='Mass flux in y direction'                          , &
        units             ='Pa m s-1'                                          , &
        loc               ='lat'                                               , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        output            =merge('h0', '  ', advection)                        , &
        restart           =.false.                                             , &
        field             =this%mfy                                            )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_qmfx'                          , &
        long_name         ='Tracer mass flux in x direction'                   , &
        units             ='Pa kg kg-1 m s-1'                                  , &
        loc               ='lon'                                               , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        output            =merge('h0', '  ', advection)                        , &
        restart           =.false.                                             , &
        field             =this%qmfx                                           )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_qmfy'                          , &
        long_name         ='Tracer mass flux in y direction'                   , &
        units             ='Pa kg kg-1 m s-1'                                  , &
        loc               ='lat'                                               , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        output            =merge('h0', '  ', advection)                        , &
        restart           =.false.                                             , &
        field             =this%qmfy                                           )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_qmfz'                          , &
        long_name         ='Tracer mass flux in z direction'                   , &
        units             ='Pa kg kg-1 m s-1'                                  , &
        loc               ='lev'                                               , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        output            =merge('h0', '  ', advection)                        , &
        restart           =.false.                                             , &
        field             =this%qmfz                                           )
      select case (this%scheme_h)
      case ('ffsl')
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_mfx_frac'                      , &
          long_name       ='Fractional mass flux in x direction'               , &
          units           ='Pa m-1 s-1'                                        , &
          loc             ='lon'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          output          ='h0'                                                , &
          restart         =.false.                                             , &
          field           =this%mfx_frac                                       )
        if (.not. this%passive) then
          call append_field(this%fields                                        , &
            name          =trim(this%name) // '_u_frac'                        , &
            long_name     ='Fractional U wind component'                       , &
            units         ='m s-1'                                             , &
            loc           ='lon'                                               , &
            mesh          =mesh                                                , &
            halo          =halo                                                , &
            output        ='h0'                                                , &
            restart       =.false.                                             , &
            field         =this%u_frac                                         )
          if (advection .and. nlev > 1) then
            call append_field(this%fields                                      , &
              name        =trim(this%name) // '_w_frac'                        , &
              long_name   ='Fractional vertical coordinate velocity'           , &
              units       ='s-1'                                               , &
              loc         ='lev'                                               , &
              mesh        =mesh                                                , &
              halo        =halo                                                , &
              output      ='h0'                                                , &
              restart     =.false.                                             , &
              field       =this%w_frac                                         )
          end if
        end if
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_cflx'                          , &
          long_name       ='CFL number in x direction'                         , &
          units           =''                                                  , &
          loc             ='lon'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          output          =merge('h0', '  ', advection)                        , &
          restart         =.false.                                             , &
          field           =this%cflx                                           )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_cfly'                          , &
          long_name       ='CFL number in y direction'                         , &
          units           =''                                                  , &
          loc             ='lat'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          output          =merge('h0', '  ', advection)                        , &
          restart         =.false.                                             , &
          field           =this%cfly                                           )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_divx'                          , &
          long_name       ='Mass flux divergence in x direction'               , &
          units           ='Pa s-1'                                            , &
          loc             ='lon'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%divx                                           )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_divy'                          , &
          long_name       ='Mass flux divergence in y direction'               , &
          units           ='Pa s-1'                                            , &
          loc             ='lat'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%divy                                           )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_qmfx0'                         , &
          long_name       ='Inner tracer mass flux in x direction'             , &
          units           ='Pa kg kg-1 m s-1'                                  , &
          loc             ='lon'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%qmfx0                                          )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_qmfy0'                         , &
          long_name       ='Inner tracer mass flux in y direction'             , &
          units           ='Pa kg kg-1 m s-1'                                  , &
          loc             ='lat'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%qmfy0                                          )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_qx'                            , &
          long_name       ='Tracer dry mixing ratio after advection in x direction', &
          units           ='kg kg-1'                                           , &
          loc             ='cell'                                              , &
          mesh            =filter_mesh                                         , &
          halo            =filter_halo                                         , &
          restart         =.false.                                             , &
          halo_cross_pole =.true.                                              , &
          field           =this%qx                                             )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_qy'                            , &
          long_name       ='Tracer dry mixing ratio after advection in y direction', &
          units           ='kg kg-1'                                           , &
          loc             ='cell'                                              , &
          mesh            =filter_mesh                                         , &
          halo            =filter_halo                                         , &
          restart         =.false.                                             , &
          field           =this%qy                                             )
      end select
      select case (this%scheme_v)
      case ('ffsl')
        if (this%passive) then
          call append_field(this%fields                                        , &
            name          =trim(this%name) // '_mfz_frac'                      , &
            long_name     ='Fractional vertical mass flux'                     , &
            units         ='Pa m-2 s-1'                                        , &
            loc           ='lev'                                               , &
            mesh          =mesh                                                , &
            halo          =halo                                                , &
            output        ='h0'                                                , &
            restart       =.false.                                             , &
            field         =this%mfz_frac                                        )
        end if
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_cflz'                          , &
          long_name       ='CFL number in z direction'                         , &
          units           =''                                                  , &
          loc             ='lev'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          output          =merge('h0', '  ', advection)                        , &
          restart         =.false.                                             , &
          field           =this%cflz                                           )
      end select
    case ('lev')
      ! Only for nonhydrostatic dynamic calculation.
      call append_field(this%fields                                            , &
        name            =trim(this%name) // '_m'                               , &
        long_name       ='Dry-air weight'                                      , &
        units           ='Pa'                                                  , &
        loc             ='lev'                                                 , &
        mesh            =filter_mesh                                           , &
        halo            =filter_halo                                           , &
        restart         =.false.                                               , &
        field           =this%m                                                )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_u'                             , &
        long_name         ='U wind component'                                  , &
        units             ='m s-1'                                             , &
        loc               ='lev_lon'                                           , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%u                                              )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_v'                             , &
        long_name         ='V wind component'                                  , &
        units             ='m s-1'                                             , &
        loc               ='lev_lat'                                           , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%v                                              )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_mfz'                           , &
        long_name         ='Vertical mass flux'                                , &
        units             ='Pa s-1'                                            , &
        loc               ='cell'                                              , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%mfz                                            )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_mfx'                           , &
        long_name         ='Mass flux in x direction'                          , &
        units             ='Pa m s-1'                                          , &
        loc               ='lev_lon'                                           , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%mfx                                            )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_mfy'                           , &
        long_name         ='Mass flux in y direction'                          , &
        units             ='Pa m s-1'                                          , &
        loc               ='lev_lat'                                           , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%mfy                                            )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_qmfx'                          , &
        long_name         ='Tracer mass flux in x direction'                   , &
        units             ='Pa kg kg-1 m s-1'                                  , &
        loc               ='lev_lon'                                           , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%qmfx                                           )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_qmfy'                          , &
        long_name         ='Tracer mass flux in y direction'                   , &
        units             ='Pa kg kg-1 m s-1'                                  , &
        loc               ='lev_lat'                                           , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%qmfy                                           )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_qmfz'                          , &
        long_name         ='Tracer mass flux in z direction'                   , &
        units             ='Pa kg kg-1 m s-1'                                  , &
        loc               ='cell'                                              , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%qmfz                                           )
      select case (this%scheme_h)
      case ('ffsl')
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_mfx_frac'                      , &
          long_name       ='Fractional mass flux in x direction'               , &
          units           ='Pa m-1 s-1'                                        , &
          loc             ='lev_lon'                                           , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%mfx_frac                                       )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_cflx'                          , &
          long_name       ='CFL number in x direction'                         , &
          units           =''                                                  , &
          loc             ='lev_lon'                                           , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%cflx                                           )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_cfly'                          , &
          long_name       ='CFL number in y direction'                         , &
          units           =''                                                  , &
          loc             ='lev_lat'                                           , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%cfly                                           )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_divx'                          , &
          long_name       ='Mass flux divergence in x direction'               , &
          units           ='Pa s-1'                                            , &
          loc             ='lev_lon'                                           , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%divx                                           )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_divy'                          , &
          long_name       ='Mass flux divergence in y direction'               , &
          units           ='Pa s-1'                                            , &
          loc             ='lev_lat'                                           , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%divy                                           )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_qmfx0'                         , &
          long_name       ='Inner tracer mass flux in x direction'             , &
          units           ='Pa kg kg-1 m s-1'                                  , &
          loc             ='lev_lon'                                           , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%qmfx0                                          )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_qmfy0'                         , &
          long_name       ='Inner tracer mass flux in y direction'             , &
          units           ='Pa kg kg-1 m s-1'                                  , &
          loc             ='lev_lat'                                           , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%qmfy0                                          )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_qx'                            , &
          long_name       ='Tracer dry mixing ratio after advection in x direction', &
          units           ='kg kg-1'                                           , &
          loc             ='lev'                                               , &
          mesh            =filter_mesh                                         , &
          halo            =filter_halo                                         , &
          restart         =.false.                                             , &
          halo_cross_pole =.true.                                              , &
          field           =this%qx                                             )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_qy'                            , &
          long_name       ='Tracer dry mixing ratio after advection in y direction', &
          units           ='kg kg-1'                                           , &
          loc             ='lev'                                               , &
          mesh            =filter_mesh                                         , &
          halo            =filter_halo                                         , &
          restart         =.false.                                             , &
          field           =this%qy                                             )
      end select
      select case (this%scheme_v)
      case ('ffsl')
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_mfz_frac'                      , &
          long_name       ='Fractional vertical mass flux'                     , &
          units           ='Pa m-2 s-1'                                        , &
          loc             ='cell'                                              , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%mfz_frac                                       )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_cflz'                          , &
          long_name       ='CFL number in z direction'                         , &
          units           =''                                                  , &
          loc             ='cell'                                              , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%cflz                                           )
      end select
    case default
      call log_error('Invalid grid location ' // trim(batch_loc) // '!', __FILE__, __LINE__)
    end select

    if (present(idx)) then
      this%ntracers = size(idx)
      allocate(this%idx(this%ntracers))
      this%idx = idx
    end if

    call time_add_alert(batch_name, seconds=dt)

    this%initialized = .true.

  end subroutine adv_batch_init

  subroutine adv_batch_clear(this)

    class(adv_batch_type), intent(inout) :: this

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

    if (allocated(this%idx)) deallocate(this%idx)

    this%loc         = 'cell'
    this%name        = ''
    this%dt          = 0
    this%initialized = .false.
    this%dynamic     = .false.
    this%ntracers    = 0
    this%nstep       = 0
    this%step        = 0

  end subroutine adv_batch_clear

  subroutine adv_batch_copy_m_old(this, m)

    class(adv_batch_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: m

    call this%m%copy(m)
    call adv_fill_vhalo(this%m, no_negvals=.true.)
    call fill_halo(this%m)

  end subroutine adv_batch_copy_m_old

  subroutine adv_batch_set_wind(this, u, v, w, mfx, mfy, mfz, m, dt)

    class(adv_batch_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: u
    type(latlon_field3d_type), intent(in) :: v
    type(latlon_field3d_type), intent(in), optional :: w
    type(latlon_field3d_type), intent(in), optional :: mfx
    type(latlon_field3d_type), intent(in), optional :: mfy
    type(latlon_field3d_type), intent(in), optional :: mfz
    type(latlon_field3d_type), intent(in), optional :: m
    real(r8), intent(in), optional :: dt

    call this%u%link(u)
    call this%v%link(v)
    if (present(w  )) call this%w  %link(w  )
    if (present(mfx)) call this%mfx%link(mfx)
    if (present(mfy)) call this%mfy%link(mfy)
    if (present(mfz)) call this%mfz%link(mfz)
    if (present(m  )) call this%copy_m_old(m)

    if (this%scheme_h == 'ffsl') call this%prepare_ffsl_h(dt)
    if (this%scheme_v == 'ffsl') call this%prepare_ffsl_v(dt)

  end subroutine adv_batch_set_wind

  subroutine adv_batch_accum_wind(this, u, v, mfx, mfy, mfz, dt)

    ! FIXME: We do not need to accumulate u and v.

    class(adv_batch_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: u
    type(latlon_field3d_type), intent(in) :: v
    type(latlon_field3d_type), intent(in) :: mfx
    type(latlon_field3d_type), intent(in) :: mfy
    type(latlon_field3d_type), intent(in), optional :: mfz
    real(r8), intent(in), optional :: dt

    integer i, j, k

    if (this%step == 0) then
      ! Reset step.
      this%u  %d = 0
      this%v  %d = 0
      this%mfx%d = 0
      this%mfy%d = 0
      if (present(mfz)) this%mfz%d = 0
      this%step = 1
    end if
    if (this%step == this%nstep) then
      ! This is the end step.
      this%u   %d = (this%u  %d + u  %d) / this%nstep
      this%v   %d = (this%v  %d + v  %d) / this%nstep
      this%mfx %d = (this%mfx%d + mfx%d) / this%nstep
      this%mfy %d = (this%mfy%d + mfy%d) / this%nstep
      if (present(mfz)) this%mfz%d = (this%mfz%d + mfz%d) / this%nstep
    else
      ! Accumulating.
      this%u  %d = this%u  %d + u  %d
      this%v  %d = this%v  %d + v  %d
      this%mfx%d = this%mfx%d + mfx%d
      this%mfy%d = this%mfy%d + mfy%d
      if (present(mfz)) this%mfz%d = this%mfz%d + mfz%d
    end if
    this%step = this%step + 1
    if (this%step > this%nstep) then
      this%step = 0
      if (this%scheme_h == 'ffsl') call this%prepare_ffsl_h(dt)
      if (this%scheme_v == 'ffsl') call this%prepare_ffsl_v(dt)
    end if

  end subroutine adv_batch_accum_wind

  subroutine adv_batch_calc_cflxy_mass(this, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: dt

    integer i, j, k

    associate (mesh   => this%mesh  , &
               u      => this%u     , & ! in
               v      => this%v     , & ! in
               cflx   => this%cflx  , & ! out
               cfly   => this%cfly  , & ! out
               u_frac => this%u_frac)   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids - 1, mesh%half_ide
          cflx%d(i,j,k) = u%d(i,j,k) * dt / mesh%de_lon(j)
          u_frac%d(i,j,k) = u%d(i,j,k) - int(cflx%d(i,j,k)) * mesh%de_lon(j) / dt
        end do
      end do
      do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          cfly%d(i,j,k) = v%d(i,j,k) * dt / mesh%de_lat(j)
        end do
      end do
    end do
    end associate

  end subroutine adv_batch_calc_cflxy_mass

  subroutine adv_batch_calc_cflxy_tracer(this, mx, my, mfx, mfy, cflx, cfly, mfx_frac, dt)

    class(adv_batch_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: mx
    type(latlon_field3d_type), intent(in) :: my
    type(latlon_field3d_type), intent(in) :: mfx
    type(latlon_field3d_type), intent(in) :: mfy
    type(latlon_field3d_type), intent(inout) :: cflx
    type(latlon_field3d_type), intent(inout) :: cfly
    type(latlon_field3d_type), intent(inout) :: mfx_frac
    real(r8), intent(in) :: dt

    real(r8) dm
    integer ks, ke, i, j, k, l

    associate (mesh => this%mesh)
    select case (this%loc)
    case ('cell', 'lev')
      ks = merge(mesh%full_kds, mesh%half_kds, this%loc == 'cell')
      ke = merge(mesh%full_kde, mesh%half_kde, this%loc == 'cell')
      do k = ks, ke
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids - 1, mesh%half_ide
            dm = mfx%d(i,j,k) * mesh%le_lon(j) * dt
            mfx_frac%d(i,j,k) = mfx%d(i,j,k)
            if (dm >= 0) then
              do l = i, mesh%full_ims, -1
                if (dm < mx%d(l,j,k) * mesh%area_cell(j)) exit
                dm = dm - mx%d(l,j,k) * mesh%area_cell(j)
                mfx_frac%d(i,j,k) = mfx_frac%d(i,j,k) - mx%d(l,j,k) * mesh%de_lon(j) / dt
              end do
              cflx%d(i,j,k) = i - l + dm / mx%d(l,j,k) / mesh%area_cell(j)
            else
              do l = i + 1, mesh%full_ime
                if (dm > -mx%d(l,j,k) * mesh%area_cell(j)) exit
                dm = dm + mx%d(l,j,k) * mesh%area_cell(j)
                mfx_frac%d(i,j,k) = mfx_frac%d(i,j,k) + mx%d(l,j,k) * mesh%de_lon(j) / dt
              end do
              cflx%d(i,j,k) = i + 1 - l + dm / mx%d(l,j,k) / mesh%area_cell(j)
            end if
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            dm = mfy%d(i,j,k) * mesh%le_lat(j) * dt
            if (dm >= 0) then
              do l = j, mesh%full_jms, -1
                if (dm < my%d(i,l,k) * mesh%area_cell(l)) exit
                dm = dm - my%d(i,l,k) * mesh%area_cell(l)
              end do
              cfly%d(i,j,k) = j - l + dm / my%d(i,l,k) / mesh%area_cell(l)
            else
              do l = j + 1, mesh%full_jme
                if (dm > -my%d(i,l,k) * mesh%area_cell(l)) exit
                dm = dm + my%d(i,l,k) * mesh%area_cell(l)
              end do
              cfly%d(i,j,k) = j + 1 - l + dm / my%d(i,l,k) / mesh%area_cell(l)
            end if
            ! Clip CFL number that are out of range. This should be very rare and in the polar region.
            cfly%d(i,j,k) = min(max(cfly%d(i,j,k), -1.0_r8), 1.0_r8)
          end do
        end do
      end do
    end select
    end associate

  end subroutine adv_batch_calc_cflxy_tracer

  subroutine adv_batch_calc_cflz_mass(this, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: dt

    real(r8) dz
    integer i, j, k, l

    associate (mesh   => this%mesh   , &
               w      => this%w      , & ! in
               cflz   => this%cflz   , & ! out
               w_frac => this%w_frac )   ! out
    do k = mesh%half_kds + 1, mesh%half_kde - 1
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          dz = w%d(i,j,k) * dt
          w_frac%d(i,j,k) = w%d(i,j,k)
          if (dz >= 0) then
            do l = k - 1, mesh%full_kms, -1
              if (dz < mesh%full_dlev(l)) exit
              dz = dz - mesh%full_dlev(l)
              w_frac%d(i,j,k) = w_frac%d(i,j,k) - mesh%full_dlev(l) / dt
            end do
            cflz%d(i,j,k) = k - 1 - l + dz / mesh%full_dlev(l)
          else
            do l = k, mesh%full_kme
              if (dz > -mesh%full_dlev(l)) exit
              dz = dz + mesh%full_dlev(l)
              w_frac%d(i,j,k) = w_frac%d(i,j,k) + mesh%full_dlev(l) / dt
            end do
            cflz%d(i,j,k) = k - l + dz / mesh%full_dlev(l)
          end if
        end do
      end do
    end do
    end associate

  end subroutine adv_batch_calc_cflz_mass

  subroutine adv_batch_calc_cflz_tracer(this, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: dt

    real(r8) dm
    integer i, j, k, l

    associate (mesh     => this%mesh    , &
               m        => this%m       , & ! in
               mfz      => this%mfz     , & ! in
               cflz     => this%cflz    , & ! out
               mfz_frac => this%mfz_frac)   ! out
    select case (this%loc)
    case ('cell')
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            dm = mfz%d(i,j,k) * mesh%half_dlev(k) * dt
            mfz_frac%d(i,j,k) = mfz%d(i,j,k)
            if (dm >= 0) then
              do l = k - 1, mesh%full_kms, -1
                if (dm < m%d(i,j,l) * mesh%full_dlev(l)) exit
                dm = dm - m%d(i,j,l) * mesh%full_dlev(l)
                mfz_frac%d(i,j,k) = mfz_frac%d(i,j,k) - m%d(i,j,l) / dt
              end do
              cflz%d(i,j,k) = k - 1 - l + dm / m%d(i,j,l) / mesh%full_dlev(l)
            else
              do l = k, mesh%full_kme
                if (dm > -m%d(i,j,l) * mesh%full_dlev(l)) exit
                dm = dm + m%d(i,j,l) * mesh%full_dlev(l)
                mfz_frac%d(i,j,k) = mfz_frac%d(i,j,k) + m%d(i,j,l) / dt
              end do
              cflz%d(i,j,k) = k - l + dm / m%d(i,j,l) / mesh%full_dlev(l)
            end if
          end do
        end do
      end do
    case ('lev')
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            dm = mfz%d(i,j,k) * mesh%full_dlev(k) * dt
            mfz_frac%d(i,j,k) = mfz%d(i,j,k)
            if (dm >= 0) then
              do l = k, mesh%half_kms, -1
                if (dm < m%d(i,j,l) * mesh%half_dlev(l)) exit
                dm = dm - m%d(i,j,l) * mesh%half_dlev(l)
                mfz_frac%d(i,j,k) = mfz_frac%d(i,j,k) - m%d(i,j,l) / dt
              end do
              cflz%d(i,j,k) = k - l + dm / m%d(i,j,l) / mesh%half_dlev(l)
            else
              do l = k, mesh%half_kme
                if (dm > -m%d(i,j,l) * mesh%half_dlev(l)) exit
                dm = dm + m%d(i,j,l) * mesh%half_dlev(l)
                mfz_frac%d(i,j,k) = mfz_frac%d(i,j,k) + m%d(i,j,l) / dt
              end do
              cflz%d(i,j,k) = k - l + dm / m%d(i,j,l) / mesh%half_dlev(l)
            end if
          end do
        end do
      end do
    end select
    end associate

  end subroutine adv_batch_calc_cflz_tracer

  subroutine adv_batch_prepare_ffsl_h(this, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in), optional :: dt

    real(r8) dt_opt

    dt_opt = this%dt; if (present(dt)) dt_opt = dt

    associate (m        => this%m       , & ! in
               mfx      => this%mfx     , & ! in
               mfy      => this%mfy     , & ! in
               u        => this%u       , & ! in
               v        => this%v       , & ! in
               mfx_frac => this%mfx_frac, & ! out
               cflx     => this%cflx    , & ! out
               cfly     => this%cfly    , & ! out
               divx     => this%divx    , & ! out
               divy     => this%divy    )   ! out
    ! Calculate horizontal CFL number and divergence along each axis.
    if (this%passive) then
      call this%calc_cflxy_tracer(m, m, mfx, mfy, cflx, cfly, mfx_frac, dt_opt)
    else
      call this%calc_cflxy_mass(dt_opt)
    end if
    call divx_operator(u, divx)
    call divy_operator(v, divy)
    end associate

  end subroutine adv_batch_prepare_ffsl_h

  subroutine adv_batch_prepare_ffsl_v(this, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in), optional :: dt

    real(r8) dt_opt, dm
    integer i, j, k, l

    dt_opt = this%dt; if (present(dt)) dt_opt = dt

    if (this%passive) then
      call this%calc_cflz_tracer(dt_opt)
    else
      call this%calc_cflz_mass(dt_opt)
    end if

  end subroutine adv_batch_prepare_ffsl_v

  subroutine adv_batch_final(this)

    type(adv_batch_type), intent(inout) :: this

    call this%clear()

  end subroutine adv_batch_final

  subroutine adv_fill_vhalo(f, no_negvals)

    type(latlon_field3d_type), intent(inout) :: f
    logical, intent(in) :: no_negvals

    integer kds, kde, kms, kme, i, j, k

    select case (f%loc)
    case ('cell')
      kds = f%mesh%full_kds
      kde = f%mesh%full_kde
      kms = f%mesh%full_kms
      kme = f%mesh%full_kme
    case ('lev')
      kds = f%mesh%half_kds
      kde = f%mesh%half_kde
      kms = f%mesh%half_kms
      kme = f%mesh%half_kme
    case default
      stop 'Unhandled branch in adv_fill_vhalo!'
    end select

    ! Set upper and lower boundary conditions.
    do k = kds - 1, kms, -1
      do j = f%mesh%full_jds, f%mesh%full_jde
        do i = f%mesh%full_ids, f%mesh%full_ide
          ! f%d(i,j,k) = f%d(i,j,kds)
          ! f%d(i,j,k) = 2 * f%d(i,j,k+1) - f%d(i,j,k+2)
          f%d(i,j,k) = 3 * f%d(i,j,k+1) - 3 * f%d(i,j,k+2) + f%d(i,j,k+3)
          ! f%d(i,j,k) = 4 * f%d(i,j,k+1) - 6 * f%d(i,j,k+2) + 4 * f%d(i,j,k+3) - f%d(i,j,k+4)
        end do
      end do
    end do
    do k = kde + 1, kme
      do j = f%mesh%full_jds, f%mesh%full_jde
        do i = f%mesh%full_ids, f%mesh%full_ide
          ! f%d(i,j,k) = f%d(i,j,kde)
          ! f%d(i,j,k) = 2 * f%d(i,j,k-1) - f%d(i,j,k-2)
          f%d(i,j,k) = 3 * f%d(i,j,k-1) - 3 * f%d(i,j,k-2) + f%d(i,j,k-3)
          ! f%d(i,j,k) = 4 * f%d(i,j,k-1) - 6 * f%d(i,j,k-2) + 4 * f%d(i,j,k-3) - f%d(i,j,k-4)
        end do
      end do
    end do
    if (no_negvals) then
      do j = f%mesh%full_jds, f%mesh%full_jde
        do i = f%mesh%full_ids, f%mesh%full_ide
          if (any(f%d(i,j,kms:kds-1) < 0)) f%d(i,j,kds-1:kms) = f%d(i,j,kds)
          if (any(f%d(i,j,kde+1:kme) < 0)) f%d(i,j,kde+1:kme) = f%d(i,j,kde)
        end do
      end do
    end if

  end subroutine adv_fill_vhalo

end module adv_batch_mod
