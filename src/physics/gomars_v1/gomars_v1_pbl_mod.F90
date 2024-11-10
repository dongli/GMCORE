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

module gomars_v1_pbl_mod

  use gomars_v1_const_mod

  implicit none

  private

  public gomars_v1_pbl_init
  public gomars_v1_pbl_final

  real(r8), parameter :: rmu    = 1.0_r8
  real(r8), parameter :: z0     = 0.01_r8
  real(r8), parameter :: epsl0  = 0.1_r8
  real(r8), parameter :: vk     = 0.4_r8
  real(r8), parameter :: alphl0 = 0.1_r8
  real(r8), parameter :: ric    = 0.195_r8
  real(r8), parameter :: factl  = 0.25_r8
  real(r8), parameter :: factm  = 1.2_r8

  real(r8), parameter :: rmu1mu = (1 - rmu) / rmu
  real(r8) dtmu

contains

  subroutine gomars_v1_pbl_init()

    dtmu = dt * rmu

  end subroutine gomars_v1_pbl_init

  subroutine gomars_v1_pbl_final()

  end subroutine gomars_v1_pbl_final

end module gomars_v1_pbl_mod