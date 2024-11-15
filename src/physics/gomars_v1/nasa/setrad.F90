
subroutine setrad(tgasref, pfgasref, co2v, co2i, qextv, qscatv, wv, gv, &
                  qexti, qscati, wi, gi, fzeroi, fzerov)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  ! Set up values used by the radiation code, such as the CO2 gas
  ! absorption coefficients.  True constants are defined, and the 
  ! time-independent quantities used by the radiation code are 
  ! calculated. 

  use gomars_v1_const_mod

  implicit none

  real(r8), intent(out) :: co2i(ntref,l_pint,l_refh2o,nspecti,ngauss)
  real(r8), intent(out) :: co2v(ntref,l_pint,l_refh2o,nspectv,ngauss)
  real(r8), intent(out) :: tgasref(ntref)
  real(r8), intent(out) :: pfgasref(l_pint)
  real(r8), intent(out) :: qextv(nspectv)
  real(r8), intent(out) :: qscatv(nspectv)
  real(r8), intent(out) :: wv(nspectv)
  real(r8), intent(out) :: gv(nspectv)
  real(r8), intent(out) :: qexti(nspecti)
  real(r8), intent(out) :: qscati(nspecti)
  real(r8), intent(out) :: wi(nspecti)
  real(r8), intent(out) :: gi(nspecti)
  real(r8), intent(out) :: fzeroi(nspecti)
  real(r8), intent(out) :: fzerov(nspectv)

  integer n, ns, nt, np, nw, ng
  real(r8) pgasref(npref)

  ! Visible dust properties:  M. Wolff Planck-weighted values (T=6000K)
  ! Log-normal size distribution:  Reff = 1.5 microns, Veff = 0.5

  ! Qext - M. Wolff values
  ! Visulal wave lengths (The order is increasing wave number)
  real(r8) qev1(nspectv)
  real(r8) qsv1(nspectv)
  real(r8) gv1 (nspectv)

  ! M. Wolff Planck-weighted values (T=215K)
  ! Infrared wave lengths.  (The order is increasing wave number.)
  real(8) qei1(nspecti)
  real(8) qsi1(nspecti)
  real(8) gi1 (nspecti)

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
  do n = 1, nspectv
    qextv (n) = qev1(n)
    qscatv(n) = qsv1(n)
    if (qscatv(n) >= qextv(n)) then
      qscatv(n) = 0.99999 * qextv(n)
    end if
    wv(n) = qscatv(n) / qextv(n)
    gv(n) = gv1(n)
  end do

  ! Fill the (INFRARED) arrays Qexti, Qscati, WI, GI
  do n = 1, nspecti
    qexti (n) = qei1(n)
    qscati(n) = qsi1(n)
    if (qscati(n) >= qexti(n)) then
      qscati(n) = 0.99999 * qexti(n)
    end if
    wi(n) = qscati(n) / qexti(n)
    gi(n) = gi1(n)
  end do

  ! Interpolate CO2 k coefficients to the finer pressure grid.
  call laginterp(ntref, npref, l_pint, nspecti, nspectv, ngauss, &
                 pgasref, pfgasref, co2i, co2v, fzeroi, fzerov)

end subroutine setrad
