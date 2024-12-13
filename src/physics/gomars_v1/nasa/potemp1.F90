subroutine potemp1(ps, tstrat, t, aadj, badj, plogadj, pl, om, tl, teta)

  use gomars_v1_const_mod
  use vert_coord_mod

  implicit none

  real(r8), intent(in ) :: ps
  real(r8), intent(in ) :: tstrat
  real(r8), intent(in ), dimension(  nlev  ) :: t
  real(r8), intent(in ), dimension(2*nlev+3) :: aadj
  real(r8), intent(in ), dimension(2*nlev+3) :: badj
  real(r8), intent(out), dimension(2*nlev+3) :: plogadj
  real(r8), intent(out), dimension(2*nlev+3) :: pl
  real(r8), intent(out), dimension(2*nlev+3) :: om
  real(r8), intent(out), dimension(2*nlev+3) :: tl
  real(r8), intent(out), dimension(2*nlev+3) :: teta

  integer n, k
  real(r8) dps, slope

  n = 2 * nlev + 3

  pl(2) = ptrop * 0.5_r8
  om(2) = (pl(2) / p0)**rd_o_cpd ! FIXME: Should we must use ps instead of p0?
  ! Set pressure on half levels.
  do k = 3, n, 2
    pl(k) = vert_coord_calc_mg_lev((k - 1) / 2, ps)
    om(k) = (pl(k) / p0)**rd_o_cpd
  end do
  ! Set pressure on full levels.
  do k = 4, n - 1, 2
    pl(k) = vert_coord_calc_mg((k - 2) / 2, ps)
    om(k) = (pl(k) / p0)**rd_o_cpd
  end do
  
  tl(2) = tstrat
  ! Set temperature on full levels.
  do k = 4, n - 1, 2
    tl(k) = t((k - 2) / 2)
  end do
  ! Set potential temperature on full levels.
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

  ! Calculate temperature from potential temperature at half levels.
  do k = 3, n, 2
    tl(k) = teta(k) * om(k)
  end do

end subroutine potemp1