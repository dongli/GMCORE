! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module cam_physics_objects_mod

  use physics_types_mod
  use cam_physics_types_mod

  implicit none

  type cam_physics_objects_type
    type(physics_mesh_type), pointer :: mesh
    type(cam_state_type) state
    type(cam_tend_type ) tend
  end type cam_physics_objects_type

  type(cam_physics_objects_type), allocatable, target :: objects(:)

contains

  subroutine cam_physics_objects_init(mesh)

    type(physics_mesh_type), intent(in), target :: mesh(:)

    integer nblk, iblk

    call cam_physics_objects_final()

    nblk = size(mesh)
    allocate(objects(nblk))
    do iblk = 1, nblk
      objects(iblk)%mesh => mesh(iblk)
      call objects(iblk)%state%init(mesh(iblk))
      call objects(iblk)%tend %init(mesh(iblk))
    end do

  end subroutine cam_physics_objects_init

  subroutine cam_physics_objects_final()

    integer iblk

    if (allocated(objects)) then
      do iblk = 1, size(objects)
        call objects(iblk)%state%clear()
        call objects(iblk)%tend %clear()
      end do
    end if

  end subroutine cam_physics_objects_final

end module cam_physics_objects_mod
