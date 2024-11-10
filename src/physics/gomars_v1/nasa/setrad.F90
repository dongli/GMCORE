
subroutine setrad(l_ntref, l_npref, l_pint, l_refh2o, l_nspecti, l_nspectv, l_ngauss, &
                  tgasref, pfgasref, co2v, co2i, qextv, qscatv, wv, gv, &
                  qexti, qscati, wi, gi, fzeroi, fzerov)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  ! PURPOSE:
  !    Set up values used by the radiation code, such as the CO2 gas
  ! absorption coefficients.  True constants are defined, and the 
  ! time-independent quantities used by the radiation code are 
  ! calculated. 
  !
  ! INPUT PARAMETERS
  ! DTAU(L,M)      - Dust optical depth of layer L, and for aerosol 
  !                  species M.
  ! ptrop          - Pressure of the tropopause (mb)
  !
  ! OUTPUT PARAMETERS
  !
  ! AEROSOL RADIATIVE OPTICAL CONSTANTS
  ! Values are at the wavelenght interval center
  !
  ! MIE SCATTERING - Size distribution weighted
  ! Qextv    - Extinction efficiency - in the visible.
  ! Qscatv   - Scattering efficiency - in the visible.
  ! WV       - Single scattering albedo - in the visible.
  ! GV       - Asymmetry parameter - in the visible.
  !
  ! Qexti    - Extinction efficiency - in the infrared.
  ! Qscati   - Scattering efficiency - in the infrared.
  ! WI       - Single scattering albedo - in the infrared.
  ! GI       - Asymmetry parameter - in the infrared.

  use gomars_v1_const_mod

  implicit none

  integer, intent(in) :: l_ntref
  integer, intent(in) :: l_npref
  integer, intent(in) :: l_pint
  integer, intent(in) :: l_refh2o
  integer, intent(in) :: l_nspecti
  integer, intent(in) :: l_nspectv
  integer, intent(in) :: l_ngauss
  real(r8), intent(out) :: co2i(l_ntref,l_pint,l_refh2o,l_nspecti,l_ngauss)
  real(r8), intent(out) :: co2v(l_ntref,l_pint,l_refh2o,l_nspectv,l_ngauss)
  real(r8), intent(out) :: tgasref(l_ntref)
  real(r8), intent(out) :: pfgasref(l_pint)
  real(r8), intent(out) :: qextv(l_nspectv)
  real(r8), intent(out) :: qscatv(l_nspectv)
  real(r8), intent(out) :: wv(l_nspectv)
  real(r8), intent(out) :: gv(l_nspectv)
  real(r8), intent(out) :: qexti(l_nspecti)
  real(r8), intent(out) :: qscati(l_nspecti)
  real(r8), intent(out) :: wi(l_nspecti)
  real(r8), intent(out) :: gi(l_nspecti)
  real(r8), intent(out) :: fzeroi(l_nspecti)
  real(r8), intent(out) :: fzerov(l_nspectv)

  integer n, ns, nt, np, nw, ng
  real(r8) pgasref(l_npref)

  ! Visible dust properties:  M. Wolff Planck-weighted values (T=6000K)
  ! Log-normal size distribution:  Reff = 1.5 microns, Veff = 0.5

  ! Qext - M. Wolff values
  ! Visulal wave lengths (The order is increasing wave number)
  real(r8) qev1(l_nspectv)
  real(r8) qsv1(l_nspectv)
  real(r8) gv1 (l_nspectv)

  ! M. Wolff Planck-weighted values (T=215K)
  ! Infrared wave lengths.  (The order is increasing wave number.)
  real(8) qei1(l_nspecti)
  real(8) qsi1(l_nspecti)
  real(8) gi1 (l_nspecti)

  qev1 = [   &
    1.834d0, &
    2.296d0, &
    2.672d0, &
    2.829d0, &
    2.698d0, &
    2.452d0, &
    2.261d0  &
  ]

  qsv1 = [   &
    1.695d0, &
    2.031d0, &
    2.583d0, &
    2.744d0, &
    2.626d0, &
    2.225d0, &
    1.525d0  &
  ]

  gv1 = [    &
    0.551d0, &
    0.640d0, &
    0.661d0, &
    0.678d0, &
    0.690d0, &
    0.743d0, &
    0.868d0  &
  ]

  qei1 = [   &
    0.008d0, &
    0.262d0, &
    0.491d0, &
    1.017d0, &
    0.444d0  &
  ]

  qsi1 = [   &
    0.001d0, &
    0.037d0, &
    0.122d0, &
    0.351d0, &
    0.336d0  &
  ]

  gi1 = [    &
    0.004D0, &
    0.030D0, &
    0.095D0, &
    0.214D0, &
    0.316D0  &
  ]

  ! Set the reference pressure and temperature arrays.  These are
  ! the pressures and temperatures at which we have k-coefficients.

  pgasref( 1) = 1.0e-6
  pgasref( 2) = 1.0e-5
  pgasref( 3) = 1.0e-4
  pgasref( 4) = 1.0e-3
  pgasref( 5) = 1.0e-2
  pgasref( 6) = 1.0e-1
  pgasref( 7) = 1.0
  pgasref( 8) = 1.0e+1
  pgasref( 9) = 1.0e+2
  pgasref(10) = 1.0e+3
  pgasref(11) = 1.0e+4

  tgasref(1)  =  50.0
  tgasref(2)  = 100.0
  tgasref(3)  = 150.0
  tgasref(4)  = 200.0
  tgasref(5)  = 250.0
  tgasref(6)  = 300.0
  tgasref(7)  = 350.0

  ! Fill the (VISIBLE) arrays Qextv, Qscatv, WV, GV
  do n = 1, l_nspectv
    qextv (n) = qev1(n)
    qscatv(n) = qsv1(n)
    if (qscatv(n) >= qextv(n)) then
      qscatv(n) = 0.99999 * qextv(n)
    end if
    wv(n) = qscatv(n) / qextv(n)
    gv(n) = gv1(n)
  end do

  ! Fill the (INFRARED) arrays Qexti, Qscati, WI, GI
  do n = 1, l_nspecti
    qexti (n) = qei1(n)
    qscati(n) = qsi1(n)
    if (qscati(n) >= qexti(n)) then
      qscati(n) = 0.99999 * qexti(n)
    end if
    wi(n) = qscati(n) / qexti(n)
    gi(n) = gi1(n)
  end do

  ! Interpolate CO2 k coefficients to the finer pressure grid.
  call laginterp(l_ntref, l_npref, l_pint, l_nspecti, l_nspectv, l_ngauss, &
                 pgasref, pfgasref, co2i, co2v, fzeroi, fzerov)

end subroutine setrad
