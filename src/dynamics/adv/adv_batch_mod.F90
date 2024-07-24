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

  ! Different tracers can be combined into one batch, and advected in different
  ! frequencies.
  type adv_batch_type
    character(30) :: scheme_h = 'N/A'
    character(30) :: scheme_v = 'N/A'
    character(10) :: loc  = 'cell'
    character(30) :: name = ''
    logical  :: dynamic   = .false.
    logical  :: use_ieva  = .false.
    integer  :: ntracers  = 1
    integer  :: nstep     = 0       ! Number of dynamic steps for one adv step
    integer  :: step      = 0       ! Step counter
    real(r8) :: dt                  ! Advection time step size in seconds
    integer , allocatable :: idx(:) ! Global index of tracers in this batch
    type(latlon_mesh_type), pointer :: mesh => null()
    type(array_type) fields
    type(latlon_field3d_type) old_m
    type(latlon_field3d_type) mfx
    type(latlon_field3d_type) mfy
    type(latlon_field3d_type) mz
    type(latlon_field3d_type) u
    type(latlon_field3d_type) v
    type(latlon_field3d_type) we     ! Explicit part of vertical mass flux
    type(latlon_field3d_type) mfx0
    type(latlon_field3d_type) mfy0
    type(latlon_field3d_type) mx0
    type(latlon_field3d_type) my0
    type(latlon_field3d_type) mz0
    type(latlon_field3d_type) dmf
    type(latlon_field2d_type) dmgs
    type(latlon_field3d_type) qmfx
    type(latlon_field3d_type) qmfy
    type(latlon_field3d_type) qmfz
    ! IEVA variables
    type(latlon_field3d_type) we_imp ! Implicit part of vertical mass flux
    type(latlon_field3d_type) cfl_h  ! Horizontal CFL number for IEVA
    real(r8), allocatable, dimension(:) :: ieva_cflx_scale
    real(r8), allocatable, dimension(:) :: ieva_a
    real(r8), allocatable, dimension(:) :: ieva_b
    real(r8), allocatable, dimension(:) :: ieva_c
    real(r8), allocatable, dimension(:) :: ieva_r
    real(r8), allocatable, dimension(:) :: ieva_qm
    ! FFSL variables
    type(latlon_field3d_type) cflx
    type(latlon_field3d_type) cfly
    type(latlon_field3d_type) cflz
    type(latlon_field3d_type) divx
    type(latlon_field3d_type) divy
    type(latlon_field3d_type) qx
    type(latlon_field3d_type) qy
  contains
    procedure :: init       => adv_batch_init
    procedure :: clear      => adv_batch_clear
    procedure :: copy_old_m => adv_batch_copy_old_m
    procedure :: set_wind   => adv_batch_set_wind
    procedure :: accum_wind => adv_batch_accum_wind
    procedure, private :: prepare_ffsl_h => adv_batch_prepare_ffsl_h
    procedure, private :: prepare_ffsl_v => adv_batch_prepare_ffsl_v
    procedure, private :: prepare_ieva => adv_batch_prepare_ieva
    final :: adv_batch_final
  end type adv_batch_type

contains

  subroutine adv_batch_init(this, filter, filter_mesh, filter_halo, mesh, halo, scheme, batch_loc, batch_name, dt, dynamic, ieva, idx)

    class(adv_batch_type), intent(inout) :: this
    type(filter_type), intent(in) :: filter
    type(latlon_mesh_type), intent(in) :: filter_mesh
    type(latlon_halo_type), intent(in) :: filter_halo(:)
    type(latlon_mesh_type), intent(in), target :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)
    character(*), intent(in) :: scheme
    character(*), intent(in) :: batch_loc
    character(*), intent(in) :: batch_name
    real(r8), intent(in) :: dt
    logical, intent(in) :: dynamic
    logical, intent(in) :: ieva
    integer, intent(in), optional :: idx(:)

    real(r8) min_width
    integer j

    call this%clear()

    this%mesh => mesh
    this%fields = array(30)

    this%loc      = batch_loc
    this%name     = batch_name
    this%dt       = dt
    this%dynamic  = dynamic
    this%use_ieva = ieva
    this%nstep    = dt / dt_dyn
    this%step     = 0

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
      if (.not. this%dynamic) then
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_old_m'                         , &
          long_name       ='Saved dry-air weight'                              , &
          units           ='Pa'                                                , &
          loc             ='cell'                                              , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.true.                                              , &
          field           =this%old_m                                          )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_mfx0'                          , &
          long_name       ='Mass flux in x direction'                          , &
          units           ='Pa m s-1'                                          , &
          loc             ='lon'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.true.                                              , &
          field           =this%mfx0                                           )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_mfy0'                          , &
          long_name       ='Mass flux in y direction'                          , &
          units           ='Pa m s-1'                                          , &
          loc             ='lat'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.true.                                              , &
          field           =this%mfy0                                           )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_mx0'                           , &
          long_name       ='Dry-air weight'                                    , &
          units           ='Pa'                                                , &
          loc             ='lon'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.true.                                              , &
          field           =this%mx0                                            )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_my0'                           , &
          long_name       ='Dry-air weight'                                    , &
          units           ='Pa'                                                , &
          loc             ='lat'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.true.                                              , &
          field           =this%my0                                            )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_mz0'                           , &
          long_name       ='Dry-air weight'                                    , &
          units           ='Pa'                                                , &
          loc             ='lev'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.true.                                              , &
          field           =this%mz0                                            )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_dmf'                           , &
          long_name       ='Mass flux divergence'                              , &
          units           ='Pa m s-1'                                          , &
          loc             ='cell'                                              , &
          mesh            =filter_mesh                                         , &
          halo            =filter_halo                                         , &
          restart         =.false.                                             , &
          field           =this%dmf                                            )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_dmgs'                          , &
          long_name       ='Surface dry-air weight tendency'                   , &
          units           ='Pa s-1'                                            , &
          loc             ='cell'                                              , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%dmgs                                           )
      end if
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_u'                             , &
        long_name         ='U wind component'                                  , &
        units             ='m s-1'                                             , &
        loc               ='lon'                                               , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%u                                              )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_v'                             , &
        long_name         ='V wind component'                                  , &
        units             ='m s-1'                                             , &
        loc               ='lat'                                               , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%v                                              )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_we'                            , &
        long_name         ='Vertical mass flux'                                , &
        units             ='Pa s-1'                                            , &
        loc               ='lev'                                               , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%we                                             )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_mfx'                           , &
        long_name         ='Mass flux in x direction'                          , &
        units             ='Pa m s-1'                                          , &
        loc               ='lon'                                               , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%mfx                                            )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_mfy'                           , &
        long_name         ='Mass flux in y direction'                          , &
        units             ='Pa m s-1'                                          , &
        loc               ='lat'                                               , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%mfy                                            )
      call append_field(this%fields                                            , &
          name            =trim(this%name) // '_mz'                            , &
          long_name       ='Dry-air weight on half level'                      , &
          units           ='Pa'                                                , &
          loc             ='lev'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%mz                                             )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_qmfx'                          , &
        long_name         ='Tracer mass flux in x direction'                   , &
        units             ='Pa kg kg-1 m s-1'                                  , &
        loc               ='lon'                                               , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%qmfx                                           )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_qmfy'                          , &
        long_name         ='Tracer mass flux in y direction'                   , &
        units             ='Pa kg kg-1 m s-1'                                  , &
        loc               ='lat'                                               , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%qmfy                                           )
      call append_field(this%fields                                            , &
        name              =trim(this%name) // '_qmfz'                          , &
        long_name         ='Tracer mass flux in z direction'                   , &
        units             ='Pa kg kg-1 m s-1'                                  , &
        loc               ='lev'                                               , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%qmfz                                           )
      if (this%use_ieva) then
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_cfl_h'                         , &
          long_name       ='Horizontal CFL number'                             , &
          units           =''                                                  , &
          loc             ='lev'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          output          ='h0'                                                , &
          restart         =.false.                                             , &
          field           =this%cfl_h                                          )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_we_imp'                        , &
          long_name       ='Implicit part of vertical mass flux'               , &
          units           ='Pa s-1'                                            , &
          loc             ='lev'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          output          ='h0'                                                , &
          field           =this%we_imp                                         )
        allocate(this%ieva_cflx_scale(mesh%full_jds:mesh%full_jde))
        allocate(this%ieva_a         (mesh%full_kds:mesh%full_kde))
        allocate(this%ieva_b         (mesh%full_kds:mesh%full_kde))
        allocate(this%ieva_c         (mesh%full_kds:mesh%full_kde))
        allocate(this%ieva_r         (mesh%full_kds:mesh%full_kde))
        allocate(this%ieva_qm        (mesh%full_kds:mesh%full_kde))
      end if
      select case (this%scheme_h)
      case ('ffsl')
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_cflx'                          , &
          long_name       ='CFL number in x direction'                         , &
          units           =''                                                  , &
          loc             ='lon'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%cflx                                           )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_cfly'                          , &
          long_name       ='CFL number in y direction'                         , &
          units           =''                                                  , &
          loc             ='lat'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
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
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_cflz'                          , &
          long_name       ='CFL number in z direction'                         , &
          units           =''                                                  , &
          loc             ='lev'                                               , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%cflz                                           )
      end select
    case ('lev')
      ! Only for nonhydrostatic dynamic calculation.
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
        name              =trim(this%name) // '_we'                            , &
        long_name         ='Vertical mass flux'                                , &
        units             ='Pa s-1'                                            , &
        loc               ='cell'                                              , &
        mesh              =mesh                                                , &
        halo              =halo                                                , &
        restart           =.false.                                             , &
        field             =this%we                                             )
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
          name            =trim(this%name) // '_mz'                            , &
          long_name       ='Dry-air weight on half level'                      , &
          units           ='Pa'                                                , &
          loc             ='cell'                                              , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          field           =this%mz                                             )
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
      if (this%use_ieva) then
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_cfl_h'                         , &
          long_name       ='Horizontal CFL number'                             , &
          units           =''                                                  , &
          loc             ='cell'                                              , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          output          ='h0'                                                , &
          restart         =.false.                                             , &
          field           =this%cfl_h                                          )
        call append_field(this%fields                                          , &
          name            =trim(this%name) // '_we_imp'                        , &
          long_name       ='Implicit part of vertical mass flux'               , &
          units           ='Pa s-1'                                            , &
          loc             ='cell'                                              , &
          mesh            =mesh                                                , &
          halo            =halo                                                , &
          restart         =.false.                                             , &
          output          ='h0'                                                , &
          field           =this%we_imp                                         )
        allocate(this%ieva_cflx_scale(mesh%full_jds:mesh%full_jde))
        allocate(this%ieva_a         (mesh%half_kds:mesh%half_kde))
        allocate(this%ieva_b         (mesh%half_kds:mesh%half_kde))
        allocate(this%ieva_c         (mesh%half_kds:mesh%half_kde))
        allocate(this%ieva_r         (mesh%half_kds:mesh%half_kde))
        allocate(this%ieva_qm        (mesh%half_kds:mesh%half_kde))
      end if
      select case (this%scheme_h)
      case ('ffsl')
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

    if (this%use_ieva) then
      min_width = minval(filter%width_lon(mesh%full_jds_no_pole:mesh%full_jde_no_pole))
      call global_min(proc%comm_model, min_width)
      this%ieva_cflx_scale = 1
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        this%ieva_cflx_scale(j) = min_width / filter%width_lon(j)
      end do
    end if

    call time_add_alert(batch_name, seconds=dt)

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

    if (allocated(this%idx            )) deallocate(this%idx            )
    if (allocated(this%ieva_cflx_scale)) deallocate(this%ieva_cflx_scale)
    if (allocated(this%ieva_a         )) deallocate(this%ieva_a         )
    if (allocated(this%ieva_b         )) deallocate(this%ieva_b         )
    if (allocated(this%ieva_c         )) deallocate(this%ieva_c         )
    if (allocated(this%ieva_r         )) deallocate(this%ieva_r         )
    if (allocated(this%ieva_qm        )) deallocate(this%ieva_qm        )

    this%loc         = 'cell'
    this%name        = ''
    this%dt          = 0
    this%dynamic     = .false.
    this%ntracers    = 0
    this%nstep       = 0
    this%step        = 0

  end subroutine adv_batch_clear

  subroutine adv_batch_copy_old_m(this, m)

    class(adv_batch_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: m

    this%old_m%d = m%d

  end subroutine adv_batch_copy_old_m

  subroutine adv_batch_set_wind(this, u_lon, v_lat, we_lev, mfx_lon, mfy_lat, dmg_lev, dt)

    class(adv_batch_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: u_lon
    type(latlon_field3d_type), intent(in) :: v_lat
    type(latlon_field3d_type), intent(in) :: we_lev
    type(latlon_field3d_type), intent(in) :: mfx_lon
    type(latlon_field3d_type), intent(in) :: mfy_lat
    type(latlon_field3d_type), intent(in) :: dmg_lev
    real(r8), intent(in) :: dt

    call this%u  %link(u_lon  )
    call this%v  %link(v_lat  )
    call this%we %link(we_lev )
    call this%mfx%link(mfx_lon)
    call this%mfy%link(mfy_lat)
    call this%mz %link(dmg_lev)

    if (this%use_ieva) call this%prepare_ieva()
    if (this%scheme_h == 'ffsl') call this%prepare_ffsl_h(dt)
    if (this%scheme_v == 'ffsl') call this%prepare_ffsl_v(dt)

  end subroutine adv_batch_set_wind

  subroutine adv_batch_accum_wind(this, dmg_lon, dmg_lat, dmg_lev, mfx_lon, mfy_lat)

    class(adv_batch_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: dmg_lon
    type(latlon_field3d_type), intent(in) :: dmg_lat
    type(latlon_field3d_type), intent(in) :: dmg_lev
    type(latlon_field3d_type), intent(in) :: mfx_lon
    type(latlon_field3d_type), intent(in) :: mfy_lat

    integer i, j, k

    if (this%step == -1) then
      ! Reset step.
      this%mfx%d = this%mfx0%d
      this%mfy%d = this%mfy0%d
      this%u  %d = this%mx0 %d
      this%v  %d = this%my0 %d
      this%mz %d = this%mz0 %d
      this%step = 1
    end if
    if (this%step == 0) then
      ! This is the first step.
      this%mfx%d = mfx_lon%d
      this%mfy%d = mfy_lat%d
      this%u  %d = dmg_lon%d
      this%v  %d = dmg_lat%d
      this%mz %d = dmg_lev%d
    else if (this%step == this%nstep) then
      ! This is the end step.
      this%mfx %d = (this%mfx%d + mfx_lon%d) / (this%nstep + 1)
      this%mfy %d = (this%mfy%d + mfy_lat%d) / (this%nstep + 1)
      this%u   %d = (this%u  %d + dmg_lon%d) / (this%nstep + 1)
      this%v   %d = (this%v  %d + dmg_lat%d) / (this%nstep + 1)
      this%mz  %d = (this%mz %d + dmg_lev%d) / (this%nstep + 1)
      this%mfx0%d = mfx_lon%d
      this%mfy0%d = mfy_lat%d
      this%mx0 %d = dmg_lon%d
      this%my0 %d = dmg_lat%d
      this%mz0 %d = dmg_lev%d
    else
      ! Accumulating.
      this%mfx %d = this%mfx%d + mfx_lon%d
      this%mfy %d = this%mfy%d + mfy_lat%d
      this%u   %d = this%u  %d + dmg_lon%d
      this%v   %d = this%v  %d + dmg_lat%d
      this%mz  %d = this%mz %d + dmg_lev%d
    end if
    this%step = merge(0, this%step + 1, this%dynamic)
    if (this%dynamic .or. this%step > this%nstep) then
      if (.not. this%dynamic) this%step = -1
      associate (mesh => this%u%mesh, &
                 dt   => this%dt    , &
                 mfx  => this%mfx   , &
                 mfy  => this%mfy   , &
                 mx   => this%u     , &
                 my   => this%v     , &
                 mz   => this%mz    , &
                 u    => this%u     , &
                 v    => this%v     , &
                 we   => this%we    , &
                 dmf  => this%dmf   , &
                 dmgs => this%dmgs  )
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids - 1, mesh%half_ide
            u%d(i,j,k) = mfx%d(i,j,k) / mx%d(i,j,k)
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            v%d(i,j,k) = mfy%d(i,j,k) / my%d(i,j,k)
          end do
        end do
      end do
      ! Diagnose horizontal mass flux divergence.
      call div_operator(mfx, mfy, dmf)
      ! Diagnose surface hydrostatic pressure tendency.
      dmgs%d = 0
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            dmgs%d(i,j) = dmgs%d(i,j) - dmf%d(i,j,k)
          end do
        end do
      end do
      ! Diagnose vertical mass flux.
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            we%d(i,j,k) = -vert_coord_calc_dmgdt_lev(k, dmgs%d(i,j)) - sum(dmf%d(i,j,1:k-1))
          end do
        end do
      end do
      if (this%use_ieva) call this%prepare_ieva()
      if (this%scheme_h == 'ffsl') call this%prepare_ffsl_h(dt)
      if (this%scheme_v == 'ffsl') call this%prepare_ffsl_v(dt)
      end associate
    end if

  end subroutine adv_batch_accum_wind

  subroutine adv_batch_prepare_ffsl_h(this, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: dt

    real(r8) work(this%u%mesh%full_ids:this%u%mesh%full_ide,this%u%mesh%half_nlev)
    real(r8) pole(this%u%mesh%half_nlev)
    integer ks, ke, i, j, k

    associate (mesh => this%u%mesh, &
               mfx  => this%mfx   , &
               mfy  => this%mfy   , &
               u    => this%u     , &
               v    => this%v     , &
               cflx => this%cflx  , &
               cfly => this%cfly  , &
               divx => this%divx  , &
               divy => this%divy  )
    ! Calculate horizontal CFL number and divergence along each axis.
    select case (this%loc)
    case ('cell', 'lev')
      ks = merge(mesh%full_kds, mesh%half_kds, this%loc == 'cell')
      ke = merge(mesh%full_kde, mesh%half_kde, this%loc == 'cell')
      do k = ks, ke
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids - 1, mesh%half_ide
            cflx%d(i,j,k) = u%d(i,j,k) * dt / mesh%de_lon(j)
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            cfly%d(i,j,k) = v%d(i,j,k) * dt / mesh%de_lat(j)
          end do
        end do
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide
            divx%d(i,j,k) = (u%d(i,j,k) - u%d(i-1,j,k)) * mesh%le_lon(j) / mesh%area_cell(j)
            divy%d(i,j,k) = (v%d(i,j,k) * mesh%le_lat(j) - v%d(i,j-1,k) * mesh%le_lat(j-1)) / mesh%area_cell(j)
          end do
        end do
      end do
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = v%d(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work(:,ks:ke), pole(ks:ke))
        pole(ks:ke) = pole(ks:ke) * mesh%le_lat(j) / global_mesh%area_pole_cap
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            divy%d(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = -v%d(i,j-1,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work(:,ks:ke), pole(ks:ke))
        pole(ks:ke) = pole(ks:ke) * mesh%le_lat(j-1) / global_mesh%area_pole_cap
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            divy%d(i,j,k) = pole(k)
          end do
        end do
      end if
    case ('vtx')
    end select
    end associate

  end subroutine adv_batch_prepare_ffsl_h

  subroutine adv_batch_prepare_ffsl_v(this, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: dt

    real(r8) deta
    integer i, j, k, l

    associate (mesh => this%u%mesh, &
               mz   => this%mz    , &
               we   => this%we    , &
               cflz => this%cflz  )
    select case (this%loc)
    case ('cell')
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            deta = we%d(i,j,k) / mz%d(i,j,k) * mesh%half_dlev(k) * dt
            if (deta >= 0) then
              do l = k, mesh%full_kde
                if (deta < mesh%full_dlev(l)) exit
                deta = deta - mesh%full_dlev(l)
              end do
              cflz%d(i,j,k) = l - k + deta / mesh%half_dlev(l)
            else
              do l = k - 1, mesh%full_kds, -1
                if (deta > -mesh%full_dlev(l)) exit
                deta = deta + mesh%full_dlev(l)
              end do
              cflz%d(i,j,k) = l - (k - 1) + deta / mesh%half_dlev(l)
            end if
          end do
        end do
      end do
    case ('lev')
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide
            deta = we%d(i,j,k) / mz%d(i,j,k) * mesh%full_dlev(k) * dt
            if (deta >= 0) then
              do l = k, mesh%half_kde
                if (deta < mesh%half_dlev(l)) exit
                deta = deta - mesh%half_dlev(l)
              end do
              cflz%d(i,j,k) = l - k + deta / mesh%full_dlev(l)
            else
              do l = k, mesh%half_kds, -1
                if (deta > -mesh%half_dlev(l)) exit
                deta = deta + mesh%half_dlev(l)
              end do
              cflz%d(i,j,k) = l - k + deta / mesh%full_dlev(l)
            end if
          end do
        end do
      end do
    end select
    end associate

  end subroutine adv_batch_prepare_ffsl_v

  subroutine adv_batch_prepare_ieva(this)

    class(adv_batch_type), intent(inout) :: this

    real(r8) work(this%u%mesh%full_ids:this%u%mesh%full_ide,this%u%mesh%half_nlev)
    real(r8) pole(this%u%mesh%half_nlev)
    real(r8) cfl_min, cfl_max, cfl_v, b
    integer i, j, k, ks, ke, ku

    associate (mesh   => this%u%mesh, &
               dt     => this%dt    , &
               u      => this%u     , & ! in
               v      => this%v     , & ! in
               mz     => this%mz    , & ! in
               we     => this%we    , & ! in
               cfl_h  => this%cfl_h , & ! working array
               we_imp => this%we_imp)   ! out
    ks = merge(this%u%mesh%half_kds, this%u%mesh%full_kds, we%loc == 'lev') + 1
    ke = merge(this%u%mesh%half_kde, this%u%mesh%full_kde, we%loc == 'lev') - 1
    ! Calculate horizontal CFL number.
    select case (we%loc)
    case ('lev')
      do k = ks, ke
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide
            ku = merge(k - 1, k, we%d(i,j,k) >= 0)
            cfl_h%d(i,j,k) = dt * (this%ieva_cflx_scale(j) *                             &
              (max(u%d(i,j,ku), 0.0_r8) - min(u%d(i-1,j,ku), 0.0_r8)) / mesh%de_lon(j) + &
              (max(v%d(i,j,ku), 0.0_r8) - min(v%d(i,j-1,ku), 0.0_r8)) / mesh%le_lon(j)   &
            )
          end do
        end do
      end do
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            ku = merge(k - 1, k, we%d(i,j,k) >= 0)
            work(i,k) = max(v%d(i,j,ku), 0.0_r8)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * dt * mesh%le_lat(j) / mesh%area_pole_cap
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            cfl_h%d(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            ku = merge(k - 1, k, we%d(i,j,k) >= 0)
            work(i,k) = min(v%d(i,j-1,ku), 0.0_r8)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = -pole * dt * mesh%le_lat(j-1) / mesh%area_pole_cap
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            cfl_h%d(i,j,k) = pole(k)
          end do
        end do
      end if
    case ('cell')
      do k = ks, ke
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide
            ku = merge(k, k + 1, we%d(i,j,k) >= 0)
            cfl_h%d(i,j,k) = dt * (this%ieva_cflx_scale(j) *                             &
              (max(u%d(i,j,ku), 0.0_r8) - min(u%d(i-1,j,ku), 0.0_r8)) / mesh%de_lon(j) + &
              (max(v%d(i,j,ku), 0.0_r8) - min(v%d(i,j-1,ku), 0.0_r8)) / mesh%le_lon(j)   &
            )
          end do
        end do
      end do
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            ku = merge(k, k + 1, we%d(i,j,k) >= 0)
            work(i,k) = max(v%d(i,j,ku), 0.0_r8)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j) / mesh%area_pole_cap
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            cfl_h%d(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            ku = merge(k, k + 1, we%d(i,j,k) >= 0)
            work(i,k) = min(v%d(i,j-1,ku), 0.0_r8)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = -pole * mesh%le_lat(j-1) / mesh%area_pole_cap
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            cfl_h%d(i,j,k) = pole(k)
          end do
        end do
      end if
    end select
    ! Partition vertical mass flux into explicit and implicit parts.
    do k = ks, ke
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          cfl_max = ieva_cfl_max - ieva_eps * cfl_h%d(i,j,k)
          cfl_min = ieva_cfl_min * cfl_max / ieva_cfl_max
          cfl_v = dt * abs(we%d(i,j,k) / mz%d(i,j,k))
          if (cfl_v <= cfl_min) then
            b = 1
          else if (cfl_min < cfl_v .and. cfl_v <= 2 * cfl_max - cfl_min) then
            b = 1.0_r8 / (1 + (cfl_v - cfl_min)**2 / (4 * cfl_max * (cfl_max - cfl_min)))
          else
            b = cfl_max / cfl_v
          end if
          this%we_imp%d(i,j,k) = (1 - b) * this%we%d(i,j,k)
          this%we    %d(i,j,k) = b * this%we%d(i,j,k)
          if (b < 0 .or. b > 1) then
            print *, '---', i, j, k
            print *, 'cfl_h        =', cfl_h%d(i,j,k)
            print *, 'cfl_v        =', cfl_v
            print *, 'u            =', u%d(i,j,k)
            print *, 'v            =', v%d(i,j,k)
            print *, 'cfl_max      =', cfl_max
            print *, 'cfl_min      =', cfl_min
            print *, 'b            =', b
            print *, 'we_imp       =', this%we_imp%d(i,j,k)
            print *, 'we           =', this%we%d(i,j,k)
            stop 999
          end if
          this%cfl_h%d(i,j,k) = b
        end do
      end do
    end do
    end associate

  end subroutine adv_batch_prepare_ieva

  subroutine adv_batch_final(this)

    type(adv_batch_type), intent(inout) :: this

    call this%clear()

  end subroutine adv_batch_final

end module adv_batch_mod
