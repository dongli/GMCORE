subroutine potemp2( &
  ps              , &
  tstrat          , &
  aadj            , &
  badj            , &
  plogadj         , &
  pl              , &
  om              , &
  tl              , &
  teta            )

  use gomars_v1_const_mod
  use vert_coord_mod

  implicit none

  real(r8), intent(in   ) :: ps
  real(r8), intent(in   ) :: tstrat
  real(r8), intent(in   ) :: aadj   (2*nlev+3)
  real(r8), intent(in   ) :: badj   (2*nlev+3)
  real(r8), intent(  out) :: plogadj(2*nlev+3)
  real(r8), intent(  out) :: pl     (2*nlev+3)
  real(r8), intent(  out) :: om     (2*nlev+3)
  real(r8), intent(inout) :: tl     (2*nlev+3)
  real(r8), intent(  out) :: teta   (2*nlev+3)

  integer n, k
  real(r8) dps, slope

  n = 2 * nlev + 3

  pl(2) = ptrop * 0.5_r8
  om(2) = (pl(2) / p0)**rd_o_cpd
  do k = 3, n, 2
    pl(k) = vert_coord_calc_mg_lev((k - 1) / 2, ps)
    om(k) = (pl(k) / p0)**rd_o_cpd
  end do
  do k = 4, n - 1, 2
    pl(k) = vert_coord_calc_mg((k - 2) / 2, ps)
    om(k) = (pl(k) / p0)**rd_o_cpd
  end do

  do k = 2, n - 1, 2
    teta(k) = tl(k) / om(k)
  end do

  ! Use a Taylor series expression for log(pl(k+1))-log(pl(k)).
  dps = pl(n) - psl
  do k = 3, n - 1
    plogadj(k) = aadj(k) + badj(k) * dps
  end do

  ! Interpolate to find potential temperature at layer boundaries.
  do k = 3, n - 4, 2
    slope = (teta(k+1) - teta(k-1)) / (plogadj(k) + plogadj(k-1))
    teta(k) = teta(k-1) + slope * plogadj(k-1)
  end do

  ! Extrapolate from above two midpoints.
  k = n - 2
  slope = (teta(k-1) - teta(k-3)) / (plogadj(k-2) + plogadj(k-3))
  teta(k) = teta(k-1) + slope * plogadj(k-1)
  k = n
  slope = (teta(k-1) - teta(k-2)) / plogadj(k-2)
  teta(k) = teta(k-1) + slope * plogadj(k-1)

  ! Change temperature at half levels.
  do k = 3, n, 2
    tl(k) = teta(k) * om(k)
  end do

end subroutine potemp2