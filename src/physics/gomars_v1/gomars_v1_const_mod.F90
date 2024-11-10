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

module gomars_v1_const_mod

  use mars_orbit_const_mod
  use const_mod

  implicit none

  ! Avogadro constant (mol-1)
  real(r8), parameter :: nav        = na
  ! Perfect gas constant (J K-1 mol-1)
  real(r8), parameter :: rgp        = ru
  ! Boltzmann constant (J K-1)
  real(r8), parameter :: kbz        = kb
  ! Molecular weight of water (kg mol-1)
  real(r8), parameter :: mh2o       = m_h2o
  ! Weight of a water molecule (kg)
  real(r8), parameter :: m0         = mh2o / nav
  ! Dust particle density (kg m-3)
  real(r8), parameter :: dpden_dt   = 2.5e3_r8
  ! Water ice particle density (kg m-3)
  real(r8), parameter :: dpden_ice  = 917.0_r8
  ! Volume of a water molecule (m3)
  real(r8), parameter :: vo1        = mh2o / dpden_ice
  ! Activation energy for desorption of water on a dust-like surface (J mol-1)
  real(r8), parameter :: desorp     = 0.288e-19_r8
  ! Esitimated activation energy for surface diffusion of water molecules (J mol-1)
  real(r8), parameter :: surfdif   = desorp / 10.0_r8
  ! Jump frequency of a water molecule (s-1)
  real(r8), parameter :: nus       = 1.0e13_r8
  ! Contact parameter of water ice on dust (m=cos(theta)) Franck's number
  real(r8), parameter :: mteta     = 0.95_r8
  ! Standard deviation of the dust distribution (reff = 0.5)
  real(r8), parameter :: dev_dt    = 0.63676_r8
  ! Standard deviation of the water ice distribution (reff = 0.1)
  real(r8), parameter :: dev_ice   = 0.3087_r8

  ! Pressure of tropopause (i.e. model top pressure) (Pa)
  real(r8) ptrop
  ! Reference surface pressure (Pa)
  real(r8) psl
  ! Physics time step (s)
  real(r8) dt
  ! Amount of surface ice required to change albedo (kg m-2)
  real(r8) icethresh_kgm2

  ! Number of soil layers
  integer , parameter :: nl         = 40

end module gomars_v1_const_mod