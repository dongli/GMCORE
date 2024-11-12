subroutine solarza(lon, cos_lat, sin_lat, cos_decl, sin_decl, time_of_day, cosz)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  ! Compute the cosine of the solar zenith angle.  Take into account
  ! if the sun rises or sets during the time period under consideration.

  use gomars_v1_const_mod

  implicit none

  real(r8), intent(in) :: lon         ! Longitude (rad)
  real(r8), intent(in) :: cos_lat     ! Cosine of latitude
  real(r8), intent(in) :: sin_lat     ! Sine of latitude
  real(r8), intent(in) :: cos_decl    ! Cosine of solar declination angle
  real(r8), intent(in) :: sin_decl    ! Sine of solar declination angle
  real(r8), intent(in) :: time_of_day ! Time fraction of a day (1)
  real(r8), intent(out) :: cosz       ! Cosine of solar zenith angle

  real(r8) lon0, adt, c1, c2, c3, c4
  real(r8) t1, cosz1
  real(r8) t2, cosz2
  real(r8) t0, tr1, trr

  lon0 = lon - pi
  c1  = sin_decl * sin_lat
  c2  = cos_decl * cos_lat
  adt = pi2 * (dt / earth_day_seconds)
  t1  = pi2 * time_of_day
  t2  = t1 + adt
  c3  = t1 + lon0
  c4  = t2 + lon0
  cosz1 = c1 + c2 * cos(c3)
  cosz2 = c1 + c2 * cos(c4)

  if (cosz1 >= 0 .and. cosz2 >= 0) then
    ! Sun above the horizon for the entire time period.
    cosz = c1 + c2 * (sin(c4) - sin(c3)) / adt
  else if (cosz1 <= 0 .and. cosz2 <= 0) then
    ! Sun below the horizon for the entire time period.
    cosz = 0
  else
    ! Sun rises or sets during the time period.
    tr1 = acos(-c1 / c2)
    trr = tr1 - lon0
    if (trr > t2) then
      tr1 = pi2 - tr1
      trr = tr1 - lon0
      if (trr < t1) trr = trr + pi2
      if (trr > t2) trr = trr - pi2
    else if (trr < t1) then
      trr = trr + pi2
      if (.not. (trr >= t1 .and. trr <= t2)) then
        trr = trr + pi2
        if (.not. (trr >= t1 .and. trr <= t2)) then
          tr1 = pi2 - tr1
          trr = tr1 - lon0
          if (trr < t1) trr = trr + pi2
          if (trr > t2) trr = trr - pi2
        end if
      end if
    end if
    if (cosz1 < 0 .and. cosz2 > 0) then
      ! Sun rises after time t1.
      if (trr < t1 .or. trr > t2) trr = t2
      cosz = (c1 * (t2 - trr) + c2 * (sin(t2 + lon0) - sin(trr + lon0))) / adt
    else
      ! Sun sets before time t2.
      if (trr < t1 .or. trr > t2) trr = t1
      cosz = (c1 * (trr - t1) + c2 * (sin(trr + lon0) - sin(t1 + lon0))) / adt
    end if
  end if

end subroutine solarza