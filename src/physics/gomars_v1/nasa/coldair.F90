subroutine coldair(tstrat, dp, pl, tl, gt, co2ice, zin, dmadt, atmcond)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  ! Check if the air temperature has dropped below the CO2 frost point. If so,
  ! the amount of CO2 condensation is calculated and appropriate adjustments
  ! are made.

  use gomars_v1_const_mod

  implicit none

  real(r8), intent(inout) :: tstrat
  real(r8), intent(in) :: dp(nlev)
  real(r8), intent(in) :: pl(2*nlev+3)
  real(r8), intent(inout) :: tl(2*nlev+3)
  real(r8), intent(inout) :: gt
  real(r8), intent(inout) :: co2ice
  real(r8), intent(in) :: zin(nl)
  real(r8), intent(out) :: dmadt
  real(r8), intent(inout) :: atmcond(nlev)

  integer k, l
  real(r8) psat
  real(r8) tsat
  real(r8) condens
  real(r8) acondns
  real(r8) tinp
  real(r8) tgp

  condens = 0.0_r8

  ! Calculate stratospheric condensation.
  psat = pl(2)
  tsat = 3182.48_r8 / (23.3494_r8 - log(psat / 100.0_r8))
  if (tstrat < tsat) then
    acondns = cpd * (tsat - tstrat) * (ptrop / g) / xlhtc
    tstrat = tsat
    condens = condens + acondns / dt
  end if

  ! Calculate tropospheric condensation.
  do l = 1, nlev
    k = 2 * l + 2
    psat = pl(k)
    tsat = 3182.48_r8 / (23.3494_r8 - log(psat / 100.0_r8))
    if (tl(k) < tsat) then
      acondns = cpd * (tsat - tl(k)) * (dp(l) / g) / xlhtc
      tl(k) = tsat
      condens = condens + acondns / dt
      atmcond(l) = atmcond(l) + acondns
    end if
  end do

  if (condens > 0) then
    ! CO2 frost point at this surface pressure.
    psat = pl(2*nlev+3)
    tsat = 3182.48_r8 / (23.3494_r8 - log(psat / 100.0_r8))
    if (co2ice > 0) then
      ! Case 1: CO2 ice already on the ground.
      ! Add condensation to existing CO2 ice mass.
      co2ice = co2ice + dt * condens
      gt = tsat
    else
      ! Case 2: No CO2 ice on the ground; Ground is warmer.
      ! Ground temperature drops when CO2 ice sublimes on warmer surface.
      tinp = sqrdy / zin(1)
      tgp = gt - tinp * dt * condens * xlhtc
      ! Calculate how much CO2 ice sublimes on hitting the ground.
      if (tgp >= tsat) then
        ! Case 2A: All CO2 ice sublimes; No net condensation.
        condens = 0
        gt = tgp
        co2ice = 0
      else
        ! Case 2B: Ground cooled to CO2 frost point and some CO2 ice remains.
        ! Calculate how much CO2 ice remains.
        condens = condens * (tsat - tgp) / (gt - tgp)
        co2ice = condens * dt
        gt = tsat
      end if
    end if
  end if

  dmadt = condens

end subroutine coldair