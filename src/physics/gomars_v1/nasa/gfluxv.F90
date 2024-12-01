subroutine gfluxv(dtdel, tdel, taucumin, wdel, cdel, cosz, sol, als, &
                  btop, bsfc, fmidp, fmidm, diffv, fluxup, fluxdn, detau)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_rad_mod

  implicit none

  real(r8), intent(in ) :: dtdel    (nlayrad)
  real(r8), intent(in ) :: tdel     (nlevrad)
  real(r8), intent(in ) :: taucumin (2*nlev+3)
  real(r8), intent(in ) :: wdel     (nlayrad)
  real(r8), intent(in ) :: cdel     (nlayrad)
  real(r8), intent(in ) :: cosz
  real(r8), intent(in ) :: sol
  real(r8), intent(in ) :: als
  real(r8), intent(in ) :: btop
  real(r8), intent(in ) :: bsfc
  real(r8), intent(out) :: fmidp    (nlayrad)
  real(r8), intent(out) :: fmidm    (nlayrad)
  real(r8), intent(out) :: diffv
  real(r8), intent(out) :: fluxup
  real(r8), intent(out) :: fluxdn
  real(r8), intent(out) :: detau

  integer , parameter :: nlp = 101 ! Must be larger than 2*nlev+3
  real(r8), parameter :: sqrt3 = sqrt(3.0_r8)
  integer k, l
  real(r8) factor
  real(r8) tau   (nlevrad)
  real(r8) w0    (nlayrad)
  real(r8) cosbar(nlayrad)
  real(r8) dtau  (nlayrad)
  real(r8) taucum(2*nlev+3)
  real(r8) alpha (nlp)
  real(r8) lambda(nlp)
  real(r8) gamma (nlp)
  real(r8) g1    (nlp)
  real(r8) g2    (nlp)
  real(r8) g3    (nlp)
  real(r8) g4
  real(r8) denom
  real(r8) am
  real(r8) ap
  real(r8) cpm1  (nlp)
  real(r8) cmm1  (nlp)
  real(r8) cp    (nlp)
  real(r8) cm    (nlp)
  real(r8) ep
  real(r8) e1    (nlp)
  real(r8) e2    (nlp)
  real(r8) e3    (nlp)
  real(r8) e4    (nlp)
  real(r8) x1    (nlp)
  real(r8) x2    (nlp)
  real(r8) taumid
  real(r8) cpmid
  real(r8) cmmid

  ! Delta-Eddington scaling
  factor    = 1 - wdel(1) * cdel(1)**2
  tau   (1) = tdel(1) * factor
  taucum(1) = 0
  taucum(2) = taucumin(2) * factor
  taucum(3) = taucum(2) + (taucumin(3) - taucumin(2)) * factor

  do l = 1, nlayrad - 1
    factor      = 1 - wdel(l) * cdel(l)**2
    w0    (l  ) = wdel(l) * (1 - cdel(l)**2) / factor
    cosbar(l  ) = cdel(l) / (1 + cdel(l))
    dtau  (l  ) = dtdel(l) * factor
    tau   (l+1) = tau(l) + dtau(l)
    k           = 2 * l + 2
    taucum(k  ) = tau(l+1)
    taucum(k+1) = taucum(k) + (taucumin(k+1) - taucumin(k)) * factor
  end do

  l           = nlayrad
  factor      = 1 - wdel(l) * cdel(l)**2
  w0    (l  ) = wdel(l) * (1 - cdel(l)**2) / factor
  cosbar(l  ) = cdel(l) / (1 + cdel(l))
  dtau  (l  ) = dtdel(l) * factor
  tau   (l+1) = tau(l) + dtau(l)
  k           = 2 * l + 1
  taucum(k  ) = tau(l+1)
  detau       = taucum(k)

  do l = 1, nlayrad
    alpha (l) = sqrt((1 - w0(l)) / (1 - w0(l) * cosbar(l)))
    g1    (l) = (sqrt3 * 0.5_r8) * (2 - w0(l) * (1 + cosbar(l)))
    g2    (l) = (sqrt3 * 0.5_r8 * w0(l)) * (1 - cosbar(l))
    g3    (l) = 0.5_r8 * (1 - sqrt3 * cosbar(l) * cosz)
    lambda(l) = sqrt(g1(l)**2 - g2(l)**2)
    gamma (l) = (g1(l) - lambda(l)) / g2(l)
  end do

  do l = 1, nlayrad
    g4    = 1 - g3(l)
    denom = lambda(l)**2 - 1.0_r8 / cosz**2
    ! There is a potential problem if w0=0 and ubarv=cosz, then denom will
    ! vanish. This only happens physically when the scattering goes to zero.
    ! Prevent this with an if statement.
    if (denom == 0) denom = 1.0e-10_r8

    am = sol * w0(l) * (g4    * (g1(l) + 1.0_r8 / cosz) + g2(l) * g3(l)) / denom
    ap = sol * w0(l) * (g3(l) * (g1(l) - 1.0_r8 / cosz) + g2(l) * g4   ) / denom

    factor  = exp(-min(tau(l) / cosz, maxexp))
    cpm1(l) = ap * factor
    cmm1(l) = am * factor

    factor = exp(-min(tau(l+1) / cosz, maxexp))
    cp(l)  = ap * factor
    cm(l)  = am * factor
  end do

  ! Calculate the exponential terms for the tridiaonal rotated layered method.
  do l = 1, nlayrad
    ep = exp(min(taumax, lambda(l) * dtau(l))) ! Clipped exponential
    e1(l) = ep + gamma(l) / ep
    e2(l) = ep - gamma(l) / ep
    e3(l) = gamma(l) * ep + 1.0_r8 / ep
    e4(l) = gamma(l) * ep - 1.0_r8 / ep
  end do

  call dsolver(nlayrad, gamma, cp, cm, cpm1, cmm1, e1, e2, e3, e4, btop, bsfc, als, x1, x2)

  do l = 1, nlayrad - 1
    ep = exp(min(taumax, lambda(l) * (taucum(2*l+1) - taucum(2*l))))
    g4 = 1 - g3(l)
    denom = lambda(l)**2 - 1.0_r8 / cosz**2
    if (denom == 0) denom = 1.0e-10_r8
    am = sol * w0(l) * (g4    * (g1(l) + 1.0_r8 / cosz) + g2(l) * g3(l)) / denom
    ap = sol * w0(l) * (g3(l) * (g1(l) - 1.0_r8 / cosz) + g2(l) * g4   ) / denom

    taumid = taucum(2*l+1)
    cpmid = ap * exp(-min(taumid / cosz, maxexp))
    cmmid = am * exp(-min(taumid / cosz, maxexp))

    fmidp(l) = x1(l) * ep + gamma(l) * x2(l) / ep + cpmid
    fmidm(l) = x1(l) * ep * gamma(l) + x2(l) / ep + cmmid

    ! Add the direct flux to the downwelling term.
    fmidm(l) = fmidm(l) + cosz * sol * exp(-min(taumid / cosz, maxexp))
  end do

  ! Flux at the top layer
  ep = 1
  g4 = 1 - g3(1)
  denom = lambda(1)**2 - 1.0_r8 / cosz**2
  if (denom == 0) denom = 1.0e-10_r8
  am = sol * w0(1) * (g4    * (g1(1) + 1.0_r8 / cosz) + g2(1) * g3(1)) / denom
  ap = sol * w0(1) * (g3(1) * (g1(1) - 1.0_r8 / cosz) + g2(1) * g4   ) / denom

  cpmid = ap
  cmmid = am

  fluxup = x1(1) * ep + gamma(1) * x2(1) / ep + cpmid
  fluxdn = x1(1) * ep * gamma(1) + x2(1) / ep + cmmid

  ! Add the direct flux to the downwelling term.
  fluxdn = fluxdn + cosz * sol * exp(-min(taucum(1) / cosz, maxexp))

  ! This is for the "special" bottom layer, where we take dtau instead of dtau/2.
  l = nlayrad
  ep = exp(min(taumax, lambda(l) * (taucum(2*nlev+3) - taucum(2*nlev+2))))
  g4 = 1 - g3(l)
  denom = lambda(l)**2 - 1.0_r8 / cosz**2
  if (denom == 0) denom = 1.0e-10_r8
  am = sol * w0(l) * (g4    * (g1(l) + 1.0_r8 / cosz) + g2(l) * g3(l)) / denom
  ap = sol * w0(l) * (g3(l) * (g1(l) - 1.0_r8 / cosz) + g2(l) * g4   ) / denom

  taumid = min(taucum(2*nlev+3), taumax)
  cpmid = ap * exp(-min(taumid / cosz, maxexp))
  cmmid = am * exp(-min(taumid / cosz, maxexp))

  fmidp(l) = x1(l) * ep + gamma(l) * x2(l) / ep + cpmid
  fmidm(l) = x1(l) * ep * gamma(l) + x2(l) / ep + cmmid

  diffv = fmidm(l)

  fmidm(l) = fmidm(l) + cosz * sol * exp(-min(taumid / cosz, maxexp))

end subroutine gfluxv