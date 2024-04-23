! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module regrid_mod

  use container
  use namelist_mod
  use block_mod
  use latlon_mesh_mod
  use latlon_field_types_mod
  use latlon_interp_mod
  use debug_mod

  implicit none

  integer, parameter :: max_regrid_fields = 100
  integer, parameter :: level_plev = 1
  integer, parameter :: level_zlev = 2

  type regrid_type
    integer :: level = 0
    integer :: nfields = 0
    type(latlon_mesh_type) mesh
    type(latlon_field3d_type), allocatable :: fields(:)
  contains
    procedure :: init => regrid_type_init
    procedure :: clear => regrid_type_clear
    final regrid_type_final
  end type regrid_type

  integer regrid_nlev
  real(r8), allocatable :: regrid_plev(:)

  type(latlon_mesh_type) regrid_global_mesh
  type(regrid_type), allocatable :: regrids(:)
  logical :: regrid_initialized = .false.

contains

  subroutine regrid_init()

    integer i, j

    call regrid_final()

    if (output_nlev == 0) return

    regrid_nlev = output_nlev
    allocate(regrid_plev(regrid_nlev))
    regrid_plev = output_plev_hPa * 100

    call regrid_global_mesh%init_global(nlon, nlat, regrid_nlev, lon_hw=0, lat_hw=0, lev_hw=0)
    do i = regrid_global_mesh%full_kds, regrid_global_mesh%full_kde
      regrid_global_mesh%full_lev(i) = regrid_plev(i)
    end do

    allocate(regrids(size(blocks)))

    do i = 1, size(blocks)
      call regrids(i)%init(blocks(i), regrid_global_mesh, level_plev)
      j = 0
      j = j + 1; call regrids(i)%fields(j)%init('u', 'U wind component', 'm/s', 'cell', regrids(i)%mesh)
      j = j + 1; call regrids(i)%fields(j)%init('v', 'V wind component', 'm/s', 'cell', regrids(i)%mesh)
      j = j + 1; call regrids(i)%fields(j)%init('t', 'Temperature'     , 'K'  , 'cell', regrids(i)%mesh)
      regrids(i)%nfields = j
    end do

    regrid_initialized = .true.

  end subroutine regrid_init

  subroutine regrid_final()

    call regrid_global_mesh%clear()

    if (allocated(regrid_plev)) deallocate(regrid_plev)
    if (allocated(regrids)) deallocate(regrids)

  end subroutine regrid_final

  subroutine regrid_run(itime)

    integer, intent(in) :: itime

    integer i, j

    do i = 1, size(regrids)
      associate (dstate => blocks(i)%dstate(itime))
      j = 0
      j = j + 1; call latlon_interp_plev(dstate%ph, dstate%u, regrid_plev, regrids(i)%fields(j))
      j = j + 1; call latlon_interp_plev(dstate%ph, dstate%v, regrid_plev, regrids(i)%fields(j))
      j = j + 1; call latlon_interp_plev(dstate%ph, dstate%t, regrid_plev, regrids(i)%fields(j))
      end associate
    end do

  end subroutine regrid_run

  subroutine regrid_type_init(this, block, global_mesh, level)

    class(regrid_type), intent(inout) :: this
    type(block_type), intent(in) :: block
    type(latlon_mesh_type), intent(in) :: global_mesh
    integer, intent(in) :: level

    call this%clear()

    this%level = level
    this%nfields = 0

    call this%mesh%init_from_parent(global_mesh, block%id, &
      block%mesh%full_ids, block%mesh%full_ide, &
      block%mesh%full_jds, block%mesh%full_jde)

    allocate(this%fields(max_regrid_fields))

  end subroutine regrid_type_init

  subroutine regrid_type_clear(this)

    class(regrid_type), intent(inout) :: this

    call this%mesh%clear()
    if (allocated(this%fields)) deallocate(this%fields)

  end subroutine regrid_type_clear

  subroutine regrid_type_final(this)

    type(regrid_type), intent(inout) :: this

    call this%clear()

  end subroutine regrid_type_final

end module regrid_mod
