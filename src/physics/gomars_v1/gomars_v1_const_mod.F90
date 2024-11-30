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
  use tracer_mod, only: ntracers

  implicit none

  ! Number of soil layers
  integer , parameter :: nl        = 40
  ! Number of aerosol tracers
  integer , parameter :: naer      = 5
  ! Tracer size bin for radiation
  integer , parameter :: nbin_rt   = 20
  integer , parameter :: nratio    = 15
  ! Number of dust particle sizes
  integer , parameter :: ndp       = 2
  integer , parameter :: ndp_dt    = ndp
  ! Number of spectral intervals in the infrared band
  integer , parameter :: nspecti   = 5
  ! Number of spectral intervals in the visible band
  integer , parameter :: nspectv   = 7
  ! Number of Gauss quadrature points
  integer , parameter :: ngauss    = 17
  integer , parameter :: npref     = 11
  integer , parameter :: ntref     = 7
  integer , parameter :: l_taumax  = 35
  integer , parameter :: l_pint    = 51
  integer , parameter :: l_refh2o  = 10
  integer , parameter :: nrefi     = 4
  integer , parameter :: nrefv     = 6
  ! Number of vertical layers
  integer :: nlev = 0

  ! Avogadro constant (mol-1)
  real(r8), parameter :: nav        = na
  ! Perfect gas constant (J K-1 mol-1)
  real(r8), parameter :: rgp        = ru
  ! Boltzmann constant (J K-1)
  real(r8), parameter :: kbz        = kb
  ! Molecular weight of water (kg mol-1)
  real(r8), parameter :: mh2o       = m_h2o
  ! Weight of a water molecule (kg)
  real(r8), parameter :: m0         = m_h2o / nav
  ! Ratio between air mass and water mass
  real(r8), parameter :: mwratio    = m_co2 / m_h2o
  ! Dust particle density (kg m-3)
  real(r8), parameter :: dpden_dt   = 2.5e3_r8
  !
  real(r8), parameter :: scale_dt   = 3.0_r8 * kb / (amu * m_co2)
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
  ! Albedo of surface ice at northern hemisphere
  real(r8), parameter :: alicen    = 0.600_r8
  ! Albedo of surface ice at southern hemisphere
  real(r8), parameter :: alices    = 0.500_r8
  ! Emissivity of bare ground outsidt 15 micron band width
  real(r8), parameter :: egognd    = 1.0_r8
  ! Emissivity of bare ground inside 15 micron band width
  real(r8), parameter :: eg15gnd   = 1.0_r8
  !
  real(r8), parameter :: eg15co2s  = 1.0_r8
  !
  real(r8), parameter :: eg15co2n  = 0.8_r8
  !
  real(r8), parameter :: egoco2s   = 1.0_r8
  !
  real(r8), parameter :: egoco2n   = 0.8_r8
  ! CO2 latent heat
  real(r8), parameter :: xlhtc     = 5.902e5_r8
  ! A radiation code conversion factor
  real(r8), parameter :: cmk       = 3.51e22_r8
  !
  real(r8), parameter :: scaveff   = 0.6_r8
  ! 1/3
  real(r8), parameter :: athird    = 1.0_r8 / 3.0_r8

  ! Pressure of tropopause (i.e. model top pressure) (Pa)
  real(r8) ptrop
  ! Reference surface pressure (Pa)
  real(r8) psl
  ! Physics time step (s)
  real(r8) dt
  ! Amount of surface ice required to change albedo (kg m-2)
  real(r8) icethresh_kgm2

  real(r8), parameter :: sqrdy = sqrt(4 * pi / mars_sol_seconds)

end module gomars_v1_const_mod