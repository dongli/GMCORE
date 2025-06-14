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
! Introduction:
!
!   The solar fluxes for Mars at given time are calculated from the 1985 Wehrli
!   Standard Extraterrestrial Solar Irradiance Spectrum.
!
!   A black-body at 5780 K produces an integrated flux over these intervals of
!   1359 W m-2 at 1 AU compared to 1355 W m-2 for the obserged sun.
! ==============================================================================

module gomars_v2_solar_mod

  use gomars_v2_const_mod
  use gomars_v2_spectra_mod
  use gomars_v2_objects_mod
  use gomars_v2_orbit_mod

  implicit none

  private

  public gomars_v2_solar_init
  public gomars_v2_solar_final
  public update_solar_flux
  public fsol_spec_mars
  public fsol_mars

  ! Solar flux within each visible spectral interval at 1AU (W m-2)
  real(r8), allocatable, dimension(:) :: fsol_spec_1au
  real(r8), allocatable, dimension(:) :: fsol_spec_mars
  real(r8) fsol_1au
  real(r8) fsol_mars

contains

  subroutine gomars_v2_solar_init()

    if (spec_vs%n /= 7) then
      stop 'gomars_v2_solar_mod only matches with visible spectra with 7 bands!'
    end if

    allocate(fsol_spec_1au (spec_vs%n))
    allocate(fsol_spec_mars(spec_vs%n))

    ! Sum equals 1356 W m-2 (values from Wehrli, 1985)
    fsol_spec_1au = [12.7_r8, 24.2_r8, 54.6_r8, 145.9_r8, 354.9_r8, 657.5_r8, 106.3_r8]
    fsol_1au = sum(fsol_spec_1au)

  end subroutine gomars_v2_solar_init

  subroutine gomars_v2_solar_final()

    if (allocated(fsol_spec_1au )) deallocate(fsol_spec_1au )
    if (allocated(fsol_spec_mars)) deallocate(fsol_spec_mars)

  end subroutine gomars_v2_solar_final

  subroutine update_solar_flux(ls)

    real(r8), intent(in) :: ls ! Solar longitude (rad)

    real(r8) d2

    d2 = solar_dist(ls)**2
    fsol_spec_mars = fsol_spec_1au / d2
    fsol_mars = fsol_1au / d2

  end subroutine update_solar_flux

end module gomars_v2_solar_mod
