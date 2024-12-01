subroutine dsolver(nl, gamma, cp, cm, cpm1, cmm1, e1, e2, e3, e4, btop, bsfc, als, x1, x2)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod, only: r8

  implicit none

  integer, parameter :: nmax = 201

  integer , intent(in ) :: nl
  real(r8), intent(in ) :: gamma(nl)
  real(r8), intent(in ) :: cp   (nl)
  real(r8), intent(in ) :: cm   (nl)
  real(r8), intent(in ) :: cpm1 (nl)
  real(r8), intent(in ) :: cmm1 (nl)
  real(r8), intent(in ) :: e1   (nl)
  real(r8), intent(in ) :: e2   (nl)
  real(r8), intent(in ) :: e3   (nl)
  real(r8), intent(in ) :: e4   (nl)
  real(r8), intent(in ) :: btop
  real(r8), intent(in ) :: bsfc
  real(r8), intent(in ) :: als
  real(r8), intent(out) :: x1   (nl)
  real(r8), intent(out) :: x2   (nl)

  integer i, l

  real(r8) a(nmax), b(nmax), c(nmax), d(nmax), x(nmax)

  a(1) = 0
  b(1) = gamma(1) + 1
  c(1) = gamma(1) - 1
  d(1) = btop - cmm1(1)

  l = 0
  do i = 2, 2 * nl - 2, 2
    l    = l + 1
    a(i) = (e1(l) + e3(l)) * (gamma(l+1) - 1)
    b(i) = (e2(l) + e4(l)) * (gamma(l+1) - 1)
    c(i) = 2 * (1 - gamma(l+1)**2)
    d(i) = (gamma(l+1) - 1) * (cpm1(l+1) - cp(l)) + (1 - gamma(l+1)) * (cm(l) - cmm1(l+1))
  end do

  l = 0
  do i = 3, 2 * nl - 1, 2
    l    = l + 1
    a(i) = 2 * (1 - gamma(l)**2)
    b(i) = (e1(l) - e3(l)) * (gamma(l+1) + 1)
    c(i) = (e1(l) + e3(l)) * (gamma(i+1) - 1)
    d(i) = e3(l) * (cpm1(l+1) - cp(l)) + e1(l) * (cm(l) - cmm1(l+1))
  end do

  a(2*nl) = e1(nl) - als * e3(nl)
  b(2*nl) = e2(nl) - als * e4(nl)
  c(2*nl) = 0
  d(2*nl) = bsfc - cp(nl) + als * cm(nl)

  call dtridgl(2*nl, a, b, c, d, x)

  do l = 1, nl
    i = 2 * l
    x1(l) = x(i-1) + x(i)
    x2(l) = x(i-1) - x(i)

    if (x2(l) /= 0) then
      if (abs(x2(l) / (x(i-1) + 1.0e-20_r8)) < 1.0e-30) x2(l) = 0
    end if
  end do

end subroutine dsolver