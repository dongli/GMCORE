subroutine fillpt( &
  p              , &
  p_lev          , &
  t              , &
  t_lev          , &
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

  real(r8), intent(in ) :: p       (nlev  )
  real(r8), intent(in ) :: p_lev   (nlev+1)
  real(r8), intent(in ) :: t       (nlev  )
  real(r8), intent(in ) :: t_lev   (nlev+1)
  real(r8), intent(in ) :: tg
  real(r8), intent(in ) :: tstrat
  real(r8), intent(out) :: plev_rad(2*nlev+3)
  real(r8), intent(out) :: tlev_rad(2*nlev+3)
  real(r8), intent(out) :: pmid_rad(2*nlev+3)
  real(r8), intent(out) :: tmid_rad(2*nlev+3)

  integer n, k, l

  n = 2 * nlev + 3

  plev_rad(1:2) = ptrop * 0.5_r8 ! FIXME: Could we set plev_rad(1) to 0?
  do k = 1, nlev
    l = 2 * k + 2
    plev_rad(l) = p(k)
  end do
  do k = 1, nlev + 1
    l = 2 * k + 1
    plev_rad(l) = p_lev(k)
  end do

  tlev_rad(1:3) = tstrat
  do k = 1, nlev
    l = 2 * k + 2
    tlev_rad(l) = t(k)
  end do
  do k = 2, nlev
    l = 2 * k + 1
    tlev_rad(l) = t_lev(k)
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