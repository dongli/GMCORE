subroutine coldair( &
  tstrat          , &
  dp              , &
  pl              , &
  tl              , &
  tg              , &
  co2ice_sfc      , &
  q               , &
  tmg             , &
  tmfdns          , &
  zin             , &
  dmadt           , &
  atmcond         )

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  ! Check if the air temperature has dropped below the CO2 frost point. If so,
  ! the amount of CO2 condensation is calculated and appropriate adjustments
  ! are made.

  use formula_mod
  use gomars_v1_const_mod
  use gomars_v1_namelist_mod, only: co2scav
  use gomars_v1_tracers_mod

  implicit none

  real(r8), intent(inout) :: tstrat
  real(r8), intent(in   ) :: dp     (nlev)
  real(r8), intent(in   ) :: pl     (2*nlev+3)
  real(r8), intent(inout) :: tl     (2*nlev+3)
  real(r8), intent(inout) :: tg
  real(r8), intent(inout) :: co2ice_sfc
  real(r8), intent(inout) :: q      (nlev,ntracers)
  real(r8), intent(  out) :: tmg    (     ntracers)
  real(r8), intent(  out) :: tmfdns (     ntracers)
  real(r8), intent(in   ) :: zin    (nsoil)
  real(r8), intent(  out) :: dmadt
  real(r8), intent(inout) :: atmcond(nlev)

  integer k, l, m, n
  integer k_scavup, k_scavdn
  real(r8) tsat
  real(r8) condens
  real(r8) acondns
  real(r8) tinp
  real(r8) tgp

  n = 2 * nlev + 3

  condens = 0

  ! Calculate stratospheric condensation.
  tsat = dewpoint_temperature_mars(pl(2))
  if (tstrat < tsat) then
    acondns = cpd * (tsat - tstrat) * (ptrop / g) / xlhtc
    tstrat = tsat
    condens = condens + acondns / dt
  end if

  ! Calculate tropospheric condensation.
  do k = 1, nlev
    l = 2 * k + 2
    tsat = dewpoint_temperature_mars(pl(l))
    if (tl(l) < tsat) then
      acondns    = cpd * (tsat - tl(l)) * (dp(k) / g) / xlhtc
      tl(l)      = tsat
      condens    = condens + acondns / dt
      atmcond(k) = atmcond(k) + acondns
    end if
  end do

  if (condens > 0) then
    ! CO2 frost point at this surface pressure.
    tsat = dewpoint_temperature_mars(pl(n))
    if (co2ice_sfc > 0) then
      ! Case 1: CO2 ice already on the ground.
      ! Add condensation to existing CO2 ice mass.
      co2ice_sfc = co2ice_sfc + dt * condens
      tg = tsat
    else
      ! Case 2: No CO2 ice on the ground; Ground is warmer.
      ! Ground temperature drops when CO2 ice sublimes on warmer surface.
      tinp = sqrdy / zin(1)
      tgp = tg - tinp * dt * condens * xlhtc
      ! Calculate how much CO2 ice sublimes on hitting the ground.
      if (tgp >= tsat) then
        ! Case 2A: All CO2 ice sublimes; No net condensation.
        condens = 0
        tg = tgp
        co2ice_sfc = 0
      else
        ! Case 2B: Ground cooled to CO2 frost point and some CO2 ice remains.
        ! Calculate how much CO2 ice remains.
        condens = condens * (tsat - tgp) / (tg - tgp)
        tg = tsat
        co2ice_sfc = condens * dt
      end if
    end if
  end if

  dmadt = condens

  ! Aerosol scavenging by CO2 snow fall.
  if (co2scav) then
    k_scavup = 0
    k_scavdn = 0
    ! Determine the portion of atmosphere affected by CO2 condensation.
    do k = 1, nlev
      if (atmcond(k) > 0) then
        k_scavup = k
        exit
      end if
    end do
    do k = nlev, 1, -1
      if (atmcond(k) > 0) then
        k_scavdn = k
        exit
      end if
    end do
    if (k_scavup /= 0 .and. k_scavdn /= 0) then
      do k = k_scavup, k_scavdn
        if (k_scavdn == nlev) then
          ! If condensation occurs down to the surface, put all aerosols on the surface (in fact, only the cloud mass matters).
          tmg   (iMa_vap) = tmg   (iMa_vap) + scaveff * q(k,iMa_cld) * dp(k) / g
          tmg   (iMa_dt ) = tmg   (iMa_dt ) + scaveff * q(k,iMa_dt ) * dp(k) / g
          tmg   (iMa_cor) = tmg   (iMa_cor) + scaveff * q(k,iMa_cor) * dp(k) / g
          tmfdns(iMa_dt ) = tmfdns(iMa_dt ) + scaveff * q(k,iMa_dt ) * dp(k) / g / dt
          tmfdns(iMa_cor) = tmfdns(iMa_cor) + scaveff * q(k,iMa_cor) * dp(k) / g / dt
          tmfdns(iMa_cld) = tmfdns(iMa_cld) + scaveff * q(k,iMa_cld) * dp(k) / g / dt
        else
          ! If condensation occurs in a restricted portion, put aerosols in the highest layer unaffected by CO2 condensation.
          do m = 1, ntracers
            q(k_scavdn+1,m) = q(k_scavdn+1,m) + scaveff * q(k,m) * dp(k) / dp(k_scavdn+1)
          end do
        end if
        do m = 1, naer
          q(k,m) = q(k,m) * (1 - scaveff)
        end do
        atmcond(k) = 0
      end do
    end if
  end if

end subroutine coldair