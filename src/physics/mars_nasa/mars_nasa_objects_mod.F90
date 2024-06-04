! ==============================================================================
! This file is part of GoMars since 2023.
!
! GoMars is a Martian general circulation model developed in Institute of
! Atmospheric Physics (IAP), Chinese Academy of Sciences (CAS).
!
! GMCORE is a dynamical core for atmospheric model used in GoMars.
!
! GoMars and GMCORE are distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module mars_nasa_objects_mod

  use physics_mesh_mod
  use mars_nasa_const_mod
  use mars_nasa_physics_types_mod

  implicit none

  type mars_nasa_objects_type
    type(physics_mesh_type), pointer :: mesh
    type(mars_nasa_state_type) state
    type(mars_nasa_tend_type) tend
  end type mars_nasa_objects_type

  type(mars_nasa_objects_type), allocatable, target :: objects(:)

contains

  subroutine mars_nasa_objects_init(mesh)

    type(physics_mesh_type), intent(in), target :: mesh(:)

    integer nblk, iblk

    call mars_nasa_objects_final()

    nblk = size(mesh)
    allocate(objects(nblk))
    do iblk = 1, nblk
      objects(iblk)%mesh => mesh(iblk)
      call objects(iblk)%state%init(objects(iblk)%mesh)
      call objects(iblk)%tend %init(objects(iblk)%mesh)
    end do

  end subroutine mars_nasa_objects_init

  subroutine mars_nasa_objects_final()

    integer iblk

    if (allocated(objects)) then
      do iblk = 1, size(objects)
        call objects(iblk)%state%clear()
        call objects(iblk)%tend %clear()
      end do
    end if

  end subroutine mars_nasa_objects_final

end module mars_nasa_objects_mod
