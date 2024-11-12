
subroutine setspv(wnov, dwnv, wavev, solarf, tauray)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  ! PURPOSE:
  !    Set up the spectral intervals in the visible (solar).  Based
  ! on Chris McKay's SETSPV code.
  !
  ! INPUT PARAMETERS
  ! L_NSPECTV  - 
  !
  ! OUTPUT PARAMETERS
  ! WNOV       - Array of wavenumbers at the spectral interval
  !              center for the visible.  Array is NSPECTV
  !              elements long.
  ! DWNV       - Array of "delta wavenumber", i.e., the width,
  !              in wavenumbers (cm^-1) of each visible spectral
  !              interval.  NSPECTV elements long.
  ! WAVEV      - Array (NSPECTV elements long) of the wavelenght
  !              (in microns) at the center of each visible spectral
  !              interval.
  ! SOLARF     - Array (NSPECTV elements) of solar flux (W/M^2) in
  !              each spectral interval.  Values are for 1 AU, and
  !              are scaled to the Mars distance elsewhere.
  ! TAURAY     - Array (NSPECTV elements) of the wavelength independent
  !              part of Rayleigh Scattering.  The pressure dependent 
  !              part is computed elsewhere (OPTCV).
 
  use gomars_v1_const_mod

  implicit none

  real(r8), intent(out) :: wnov  (l_nspectv)
  real(r8), intent(out) :: dwnv  (l_nspectv)
  real(r8), intent(out) :: wavev (l_nspectv)
  real(r8), intent(out) :: solarf(l_nspectv)
  real(r8), intent(out) :: tauray(l_nspectv)

  real(r8) wl
  integer n, m

  real(8), parameter :: pray0 = 9.423d+4 ! Rayleigh scattering reference pressure in hPa.
 
  ! BWNV - Bin wavenumber of the edges of the visible spectral bins
  ! units are inverse centimeters.  Dimension needs to be changed
  ! if the number of visible bins changes.

  ! Bin wavenumber - wavenumber [cm^(-1)] at the edges of the visible
  ! spectral bins.  Go from smaller to larger wavenumbers, the same as
  ! in the IR.

  ! 2222.22D0    ->   4.50 microns
  ! 3087.37D0    ->   3.24 microns
  ! 4030.63D0    ->   2.48 microns
  ! 5370.57D0    ->   1.86 microns
  ! 7651.11D0    ->   1.31 microns
  ! 12500.00D0   ->   0.80 microns
  ! 25000.00D0   ->   0.40 microns
  ! 41666.67D0   ->   0.24 microns

  real(8) bwnv(l_nspectv+1)

  ! Solar flux within each spectral interval, at 1AU (W/M^2)
  ! Sum equals 1356 W/m^2 (values from Wehrli, 1985)

  real(8) solar(l_nspectv)

  ! Set up mean wavenumbers and wavenumber deltas.  Units of 
  ! wavenumbers is cm^(-1); units of wavelengths is microns.

  bwnv = [ &
    2222.22d0, 3087.37d0, 4030.63d0,  &
    5370.57d0, 7651.11d0, 12500.00d0, &
    25000.00d0, 41666.67d0            &
  ]

  solar = [12.7, 24.2, 54.6, 145.9, 354.9, 657.5, 106.3]

  do m = 1, l_nspectv
    wnov(m)  = 0.5 * (bwnv(m+1) + bwnv(m))
    dwnv(m)  = bwnv(m+1) - bwnv(m)
    wavev(m) = 1.0e+4 / wnov(m)
  end do

  do n = 1, l_nspectv
    solarf(n) = solar(n)
  end do

  ! Set up the wavelength independent part of Rayleigh Scattering.
  ! The pressure dependent part will be computed elsewhere (OPTCV).
  ! WAVEV is in microns.  There is no Rayleigh scattering in the IR.

  do n = 1, l_nspectv
    wl = wavev(n)
    tauray(n) = (8.7_r8 / g) * (1.527_r8 * (1.0 + 0.013 / wl**2)/ wl**4) / pray0
  end do

end subroutine setspv