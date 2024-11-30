subroutine tpindex(t, p, qh2o, coef, idx_t, idx_p, idx_h2o, wratio)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_rad_mod, only: pref => pfgasref, tref => tgasref, wrefh2o

  implicit none

  real(r8), intent(in ) :: t
  real(r8), intent(in ) :: p
  real(r8), intent(in ) :: qh2o
  real(r8), intent(out) :: coef(4)
  integer , intent(out) :: idx_t      ! Temperature-grid index
  integer , intent(out) :: idx_p      ! Pressure-grid index
  integer , intent(out) :: idx_h2o    ! Water abundance index
  real(r8), intent(out) :: wratio     ! Water abundance ratio

  integer i
  real(r8) ct, cp, plog

  ! Get the upper and lower Temperature-grid indicies that bound the
  ! requested temperature.  If the requested temperature is outside
  ! the T-grid, set up to extrapolate from the appropriate end.
  if (t <= tref(1)) then
    idx_t = 1
  else
    idx_t = 0
    do i = 1, ntref - 1
      if (t > tref(i) .and. t <= tref(i+1)) then
        idx_t = i
        exit
      end if
    end do
    if (idx_t == 0) then
      idx_t = ntref - 1
    end if
  end if
  ct = (t - tref(idx_t)) / (tref(idx_t+1) - tref(idx_t))

  ! Get the upper and lower Pressure-grid indicies that bound the
  ! requested pressure.  If the requested pressure is outside
  ! the P-grid, set up to extrapolate from the appropiate end.
  plog = log10(p)
  idx_p = 0
  do i = 2, npint - 1
    if (plog <= pref(i)) then
      idx_p = i - 1
      exit
    end if
  end do
  if (idx_p == 0) then
    idx_p = npint - 1
  end if
  cp = (plog - pref(idx_p)) / (pref(idx_p+1) - pref(idx_p))

  ! Fill the interpolation coefficients.
  coef(1) = (1 - cp) * (1 - ct)
  coef(2) = cp * (1 - ct)
  coef(3) = cp * ct
  coef(4) = (1 - cp) * ct

  ! Get the indices for water abundance. There are 10 sets of k-coefficients
  ! with differing amounts of water vs CO2.
  if (qh2o <= wrefh2o(1)) then
    idx_h2o = 1
    wratio = 0
  else if (qh2o >= wrefh2o(nrefh2o)) then
    idx_h2o = nrefh2o
    wratio = 0
  else
    do i = 2, nrefh2o
      if (qh2o >= wrefh2o(i-1) .and. qh2o < wrefh2o(i)) then
        idx_h2o = i
        wratio = (qh2o - wrefh2o(i-1)) / (wrefh2o(i) - wrefh2o(i-1))
        exit
      end if
    end do
  end if

end subroutine tpindex