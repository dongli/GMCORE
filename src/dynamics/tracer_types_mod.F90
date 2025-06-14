! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module tracer_types_mod

  use const_mod
  use namelist_mod
  use latlon_field_types_mod

  implicit none

  integer nbatches
  integer ntracers
  integer ntracers_water

  integer idx_qv
  integer idx_qc, idx_nc
  integer idx_qi, idx_ni
  integer idx_qr, idx_nr
  integer idx_qs, idx_ns
  integer idx_qg
  integer idx_qh
  integer idx_qo3
  integer idx_qso2

  character(strlen_name     ), allocatable :: batch_names(:)
  character(strlen_name     ), allocatable :: tracer_batches(:)
  character(strlen_name     ), allocatable :: tracer_names(:)
  character(strlen_long_name), allocatable :: tracer_long_names(:)
  character(strlen_units    ), allocatable :: tracer_units(:)
  real(r8), allocatable :: batch_dts(:)
  ! Follow other model's definition:
  ! 0 - Generic tracer
  ! 1 - Prognostic chemical tracer
  ! 2 - Diagnostic chemical tracer
  integer, allocatable :: tracer_types(:)
  logical, allocatable :: is_water_tracer(:)

  type tracers_type
    logical :: initialized = .false.
    type(latlon_field4d_type) q
    ! Some diagnostics:
    type(latlon_field3d_type) qm      ! Total moisture or water substances
    type(latlon_field3d_type) qm_lev  ! Total moisture or water substances on half level
  contains
    procedure :: init_stage1 => tracers_init_stage1
    procedure :: init_stage2 => tracers_init_stage2
    procedure :: clear => tracers_clear
    final :: tracers_final
  end type tracers_type

  type(tracers_type), allocatable, target :: tracers(:) ! (blocks)

contains

  subroutine tracers_init_stage1(this, filter_mesh, filter_halo, mesh, halo)

    class(tracers_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: filter_mesh
    type(latlon_halo_type), intent(in) :: filter_halo(:)
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)

    character(field_name_len     ) name
    character(field_long_name_len) long_name
    character(field_units_len    ) units

    call this%clear()

    name      = 'qm'
    long_name = 'Total moist tracer dry mixing ratioo'
    units     = 'kg kg-1'
    call this%qm%init(name, long_name, units, 'cell', mesh, halo, restart=.true.)

    if (nonhydrostatic) then
      name      = 'qm_lev'
      long_name = 'Total moist tracer dry mixing ratioo on half level'
      units     = 'kg kg-1'
      call this%qm_lev%init(name, long_name, units, 'lev', mesh, halo, restart=.true.)
    end if

    this%initialized = .true.

  end subroutine tracers_init_stage1

  subroutine tracers_init_stage2(this, filter_mesh, filter_halo, mesh, halo)

    class(tracers_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: filter_mesh
    type(latlon_halo_type), intent(in) :: filter_halo(:)
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)

    if (ntracers > 0) then
      call this%q%init(                               &
        name              ='q'                      , &
        long_name         ='Tracer dry mixing ratio', &
        units             ='kg kg-1'                , &
        loc               ='cell'                   , &
        mesh              =filter_mesh              , &
        halo              =filter_halo              , &
        halo_cross_pole   =.true.                   , &
        output            ='h0'                     , &
        restart           =.true.                   , &
        dim4_name         ='tracers'                , &
        dim4_size         =ntracers                 , &
        var4_names        =tracer_names(:ntracers)  )
    end if

  end subroutine tracers_init_stage2

  subroutine tracers_clear(this)

    class(tracers_type), intent(inout) :: this

    call this%q     %clear()
    call this%qm    %clear()
    call this%qm_lev%clear()

    this%initialized = .false.

  end subroutine tracers_clear

  subroutine tracers_final(this)

    type(tracers_type), intent(inout) :: this

    call this%clear()

  end subroutine tracers_final

end module tracer_types_mod
