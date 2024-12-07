subroutine gfluxi(is, tlev, dtau, taucum, wbar, cosb, albi, btop, bsfc, ftopup, fmidp, fmidm)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_rad_mod

  implicit none

  integer , intent(in ) :: is
  real(r8), intent(in ) :: tlev  (2*nlev+3)
  real(r8), intent(in ) :: dtau  (nlayrad )
  real(r8), intent(in ) :: taucum(2*nlev+3)
  real(r8), intent(in ) :: wbar  (nlayrad )
  real(r8), intent(in ) :: cosb  (nlayrad )
  real(r8), intent(in ) :: albi
  real(r8), intent(in ) :: btop
  real(r8), intent(in ) :: bsfc
  real(r8), intent(out) :: ftopup
  real(r8), intent(out) :: fmidp (nlayrad)
  real(r8), intent(out) :: fmidm (nlayrad)

  integer, parameter :: nlp = 101 ! Must be larger than 2*nlev+3
  integer l, nt, nt2
  real(r8) term, ep, em, dtauk, cpmid, cmmid, fluxup, fluxdn
  real(r8) alpha (nlp)
  real(r8) lambda(nlp)
  real(r8) gamma (nlp)
  real(r8) b0    (nlp)
  real(r8) b1    (nlp)
  real(r8) cp    (nlp)
  real(r8) cm    (nlp)
  real(r8) cpm1  (nlp)
  real(r8) cmm1  (nlp)
  real(r8) e1    (nlp)
  real(r8) e2    (nlp)
  real(r8) e3    (nlp)
  real(r8) e4    (nlp)
  real(r8) x1    (nlp)
  real(r8) x2    (nlp)

  do l = 1, nlayrad - 1
    alpha (l) = sqrt((1 - wbar(l)) / (1 - wbar(l) * cosb(l)))
    lambda(l) = alpha(l) * (1 - wbar(l) * cosb(l)) / ubari
    nt  = tlev(2*l  ) * 10 - 499
    nt2 = tlev(2*l+2) * 10 - 499
    b0(l) = planckir(is,nt)
    b1(l) = (planckir(is,nt2) - planckir(is,nt)) / dtau(l)
  end do

  l = nlayrad
  alpha (l) = sqrt((1 - wbar(l)) / (1 - wbar(l) * cosb(l)))
  lambda(l) = alpha(l) * (1 - wbar(l) * cosb(l)) / ubari
  nt  = tlev(2*l  ) * 10 - 499
  nt2 = tlev(2*l+1) * 10 - 499
  b0(l) = planckir(is,nt2)
  b1(l) = (planckir(is,nt) - planckir(is,nt2)) / dtau(l)

  do l = 1, nlayrad
    gamma(l) = (1 - alpha(l)) / (1 + alpha(l))
    term = ubari / (1 - wbar(l) * cosb(l))
    cp(l) = b0(l) + b1(l) * dtau(l) + b1(l) * term
    cm(l) = b0(l) + b1(l) * dtau(l) - b1(l) * term
    cpm1(l) = b0(l) + b1(l) * term
    cmm1(l) = b0(l) - b1(l) * term
  end do

  do l = 1, nlayrad
    ep = exp(min(lambda(l) * dtau(l), maxexp))
    em = 1.0_r8 / ep
    e1(l) = ep + gamma(l) * em
    e2(l) = ep - gamma(l) * em
    e3(l) = gamma(l) * ep + em
    e4(l) = gamma(l) * ep - em
  end do

  call dsolver(nlayrad, gamma, cp, cm, cpm1, cmm1, e1, e2, e3, e4, btop, bsfc, albi, x1, x2)

  do l = 1, nlayrad
    dtauk = taucum(2*l+1) - taucum(2*l)
    ep = exp(min(lambda(l) * dtauk, maxexp))
    em = 1.0_r8 / ep
    term = ubari / (1 - wbar(l) * cosb(l))
    cpmid = b0(l) + b1(l) * dtauk + b1(l) * term
    cmmid = b0(l) + b1(l) * dtauk - b1(l) * term
    fmidp(l) = x1(l) * ep + gamma(l) * x2(l) * em + cpmid
    fmidm(l) = x1(l) * ep * gamma(l) + x2(l) * em + cmmid
    ! For flux, integrate over the hemisphere treating intensity constant.
    fmidp(l) = fmidp(l) * pi
    fmidm(l) = fmidm(l) * pi
  end do

  l = nlayrad
  ep = exp(min(lambda(l) * dtau(l), taumax))
  em = 1.0_r8 / ep
  term = ubari / (1 - wbar(l) * cosb(l))
  cpmid = b0(l) + b1(l) * dtau(l) + b1(l) * term
  cmmid = b0(l) + b1(l) * dtau(l) - b1(l) * term
  fmidp(l) = x1(l) * ep + gamma(l) * x2(l) * em + cpmid
  fmidm(l) = x1(l) * ep * gamma(l) + x2(l) * em + cmmid
  fmidp(l) = fmidp(l) * pi
  fmidm(l) = fmidm(l) * pi

  ep = 1
  em = 1
  term = ubari / (1 - wbar(1) * cosb(1))
  cpmid = b0(1) + b1(1) * term
  cmmid = b0(1) - b1(1) * term
  fluxup = x1(1) * ep + gamma(1) * x2(1) * em + cpmid
  fluxdn = x1(1) * ep * gamma(1) + x2(1) * em + cmmid

  ftopup = (fluxup - fluxdn) * pi

end subroutine gfluxi