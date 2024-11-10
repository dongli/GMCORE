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

module gomars_v1_objects_mod

  use physics_mesh_mod
  use gomars_v1_types_mod

  implicit none

  private

  public physics_mesh_type
  public gomars_v1_objects_type
  public objects
  public gomars_v1_objects_init
  public gomars_v1_objects_final

  type gomars_v1_objects_type
    type(physics_mesh_type), pointer :: mesh
    type(gomars_v1_state_type) state
    type(gomars_v1_tend_type) tend
  end type gomars_v1_objects_type

  type(gomars_v1_objects_type), allocatable :: objects(:)

contains

  subroutine gomars_v1_objects_init(mesh)

    type(physics_mesh_type), intent(in), target :: mesh(:)

    integer nblk, iblk

    call gomars_v1_objects_final()

    nblk = size(mesh)
    allocate(objects(nblk))
    do iblk = 1, nblk
      objects(iblk)%mesh => mesh(iblk)
      call objects(iblk)%state%init(objects(iblk)%mesh)
      call objects(iblk)%tend %init(objects(iblk)%mesh)
    end do

  end subroutine gomars_v1_objects_init

  subroutine gomars_v1_objects_final()

    integer iblk

    if (allocated(objects)) then
      do iblk = 1, size(objects)
        call objects(iblk)%state%clear()
        call objects(iblk)%tend %clear()
      end do
      deallocate(objects)
    end if

  end subroutine gomars_v1_objects_final

end module gomars_v1_objects_mod