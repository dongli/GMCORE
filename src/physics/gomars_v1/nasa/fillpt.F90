subroutine fillpt( &
  pl             , &
  tl             , &
  tg             , &
  tstrat         , &
  plev_rad       , &
  tlev_rad       , &
  pmid_rad       , &
  tmid_rad       )

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  ! Put the pressure and temperature arrays onto the radiation grid.

  use gomars_v1_const_mod

  implicit none

  real(r8), intent(in ) :: pl      (2*nlev+3)
  real(r8), intent(in ) :: tl      (2*nlev+3)
  real(r8), intent(in ) :: tg
  real(r8), intent(in ) :: tstrat
  real(r8), intent(out) :: plev_rad(2*nlev+3)
  real(r8), intent(out) :: tlev_rad(2*nlev+3)
  real(r8), intent(out) :: pmid_rad(2*nlev+3)
  real(r8), intent(out) :: tmid_rad(2*nlev+3)

  integer n, k, l

  n = 2 * nlev + 3

  plev_rad(1) = ptrop * 0.5_r8
  plev_rad(2) = ptrop * 0.5_r8
  do l = 3, n
    plev_rad(l) = pl(l)
  end do

  tlev_rad(1:3) = tstrat
  do l = 4, n - 1, 2
    tlev_rad(l) = tl(l)
  end do
  do l = 5, n - 2, 2
    tlev_rad(l) = tlev_rad(l+1) + (tlev_rad(l-1) - tlev_rad(l+1)) * &
                               log(plev_rad(l  ) / plev_rad(l+1)) / &
                               log(plev_rad(l-1) / plev_rad(l+1))
  end do
  tlev_rad(n) = tg

  ! tmid_rad and pmid_rad are used by optci and optcv subroutines to get the
  ! index for CO2 K-coefficient interpolation.
  tmid_rad(1) = tlev_rad(2)
  tmid_rad(2) = tlev_rad(2)
  pmid_rad(1) = plev_rad(1)
  pmid_rad(2) = plev_rad(2)
  do k = 1, nlev
    l = 2 * k + 1
    tmid_rad(l  ) = tlev_rad(l)
    tmid_rad(l+1) = tlev_rad(l)
    pmid_rad(l  ) = plev_rad(l)
    pmid_rad(l+1) = plev_rad(l)
  end do
  tmid_rad(n) = tlev_rad(n)
  pmid_rad(n) = plev_rad(n)

end subroutine fillpt