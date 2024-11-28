module testbed_objects_mod

  use testbed_types_mod

  implicit none

  type testbed_objects_type
    type(physics_mesh_type), pointer :: mesh
    type(testbed_state_type) state
    type(testbed_tend_type ) tend
  end type testbed_objects_type

  type(testbed_objects_type), allocatable :: objects(:)

contains

  subroutine testbed_objects_init(mesh)

    type(physics_mesh_type), intent(in), target :: mesh(:)

    integer nblk, iblk

    call testbed_objects_final()

    nblk = size(mesh)
    allocate(objects(nblk))
    do iblk = 1, nblk
      objects(iblk)%mesh => mesh(iblk)
      call objects(iblk)%state%init(objects(iblk)%mesh)
      call objects(iblk)%tend %init(objects(iblk)%mesh)
    end do

  end subroutine testbed_objects_init

  subroutine testbed_objects_final()

    integer iblk

    if (allocated(objects)) then
      do iblk = 1, size(objects)
        call objects(iblk)%state%clear()
        call objects(iblk)%tend %clear()
      end do
    end if

  end subroutine testbed_objects_final

end module testbed_objects_mod