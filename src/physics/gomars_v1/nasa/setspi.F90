subroutine setspi(wnoi, dwni, wavei, planckir)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  ! Set up the spectral intervals in the infrared.  Based on
  ! Chris McKay's SETSPI code.

  use gomars_v1_const_mod

  implicit none

  real(r8), intent(out) :: wnoi    (nspecti)
  real(r8), intent(out) :: dwni    (nspecti)
  real(r8), intent(out) :: wavei   (nspecti)
  real(r8), intent(out) :: planckir(nspecti,8501)

  real(r8) a, b, ans, y, bpa, bma, t
  real(r8) wn1, wn2
  integer n, nw, nt, m

  ! C1 and C2 values from Goody and Yung (2nd edition)  MKS units
  ! These values lead to a "sigma" (sigma*T^4) of 5.67032E-8 W m^-2 K^-4
  real(r8), parameter :: c1 = 3.741832d-16 ! W m^-2
  real(r8), parameter :: c2 = 1.438786d-2  ! m K
  
  real(r8) :: x(12) = [-0.981560634246719D0,  -0.904117256370475D0, &
                       -0.769902674194305D0,  -0.587317954286617D0, &
                       -0.367831498998180D0,  -0.125233408511469D0, &
                        0.125233408511469D0,   0.367831498998180D0, &
                        0.587317954286617D0,   0.769902674194305D0, &
                        0.904117256370475D0,   0.981560634246719D0]

  real(r8) :: w(12) = [ 0.047175336386512D0,   0.106939325995318D0, &
                        0.160078328543346D0,   0.203167426723066D0, &
                        0.233492536538355D0,   0.249147045813403D0, &
                        0.249147045813403D0,   0.233492536538355D0, &
                        0.203167426723066D0,   0.160078328543346D0, &
                        0.106939325995318D0,   0.047175336386512D0]

  ! BWNI - Bin wavenumber of the edges of the IR spectral bins
  ! units are inverse centimeters.  Dimension needs to be changed
  ! if the number of IR bins changes.

  ! Bin wavenumber - wavenumber [cm^(-1)] at the edges of the IR
  ! spectral bins.

  real(r8) bwni(nspecti+1)

  bwni = [ &
      10.000d0, & ! 1000.0 microns
     166.667d0, & !   60.0 microns
     416.667d0, & !   24.0 microns
     833.333d0, & !   12.0 microns
    1250.000d0, & !    8.0 microns
    2222.222d0  & !    4.5 microns
  ]

  ! Set up mean wavenumbers and wavenumber deltas.  Units of 
  ! wavenumbers is cm^(-1); units of wavelengths is microns.
  do m = 1, nspecti
    wnoi (m) = 0.5 * (bwni(m+1) + bwni(m))
    dwni (m) = bwni(m+1) - bwni(m)
    wavei(m) = 1.0e+4 / wnoi(m)
  end do

  !  For each IR wavelength interval, compute the integral of B(T), the
  !  Planck function, divided by the wavelength interval, in cm-1.  The
  !  integration is in MKS units, the final answer is the same as the
  !  original planck.f; W m^-2 wavenumber^-1, where wavenumber is in CM^-1.
  do nw = 1, nspecti
    a = 1.0d-2 / bwni(nw+1)
    b = 1.0d-2 / bwni(nw)
    bpa = (b + a) / 2.0
    bma = (b - a) / 2.0
    do nt = 500, 9000
      t   = dble(nt) / 1.0d+1
      ans = 0.0d0
      do m = 1, 12
        y   = bma * x(m) + bpa
        ans = ans + w(m) * c1 / (y**5 * (exp(c2 / (y * t)) - 1.0d0))
      end do
      planckir(nw,nt-499) = ans * bma / (pi * dwni(nw))
    end do
  end do

end subroutine setspi