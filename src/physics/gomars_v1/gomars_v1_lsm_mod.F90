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

module gomars_v1_lsm_mod

  use gomars_v1_const_mod

  implicit none

  private

  public gomars_v1_lsm_init
  public gomars_v1_lsm_final
  public sthick, sdepth

  real(r8), allocatable :: sthick(:)
  real(r8), allocatable :: sdepth(:)

contains

  subroutine gomars_v1_lsm_init()

    call gomars_v1_lsm_final()

    allocate(sthick(nsoil))
    allocate(sdepth(nsoil))

  end subroutine gomars_v1_lsm_init

  subroutine gomars_v1_lsm_final()

    if (allocated(sthick)) deallocate(sthick)
    if (allocated(sdepth)) deallocate(sdepth)

  end subroutine gomars_v1_lsm_final

end module gomars_v1_lsm_mod