subroutine coldair( &
  ps              , &
  tstrat          , &
  dp              , &
  p               , &
  t               , &
  tg              , &
  co2ice_sfc      , &
  q               , &
  tm_sfc          , &
  tmflx_sfc_dn    , &
  zin             , &
  dmsdt           )

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

  real(r8), intent(in   ) :: ps
  real(r8), intent(inout) :: tstrat
  real(r8), intent(in   ) :: dp          (nlev)
  real(r8), intent(in   ) :: p           (nlev)
  real(r8), intent(inout) :: t           (nlev)
  real(r8), intent(inout) :: tg
  real(r8), intent(inout) :: co2ice_sfc
  real(r8), intent(inout) :: q           (nlev,ntracers)
  real(r8), intent(  out) :: tm_sfc      (     ntracers)
  real(r8), intent(  out) :: tmflx_sfc_dn(     ntracers)
  real(r8), intent(in   ) :: zin         (nsoil)
  real(r8), intent(  out) :: dmsdt

  integer k, m
  integer k_scavup, k_scavdn
  real(r8) tsat
  real(r8) dm(nlev) ! Condensation mass on each layer.
  real(r8) tinp
  real(r8) tgp

  dmsdt = 0

  ! Calculate stratospheric condensation.
  tsat = dewpoint_temperature_mars(ptrop * 0.5_r8)
  if (tstrat < tsat) then
    dm(1)  = cpd * (tsat - tstrat) * (ptrop / g) / xlhtc ! FIXME: Is the units kg m-1?
    tstrat = tsat
    dmsdt  = dmsdt + dm(1) / dt
  end if

  ! Calculate tropospheric condensation.
  do k = 1, nlev
    tsat = dewpoint_temperature_mars(p(k))
    if (t(k) < tsat) then
      dm(k) = cpd * (tsat - t(k)) * (dp(k) / g) / xlhtc
      t(k)  = tsat
      dmsdt = dmsdt + dm(k) / dt
    else
      dm(k) = 0
    end if
  end do

  if (dmsdt > 0) then
    ! CO2 frost point at this surface pressure.
    tsat = dewpoint_temperature_mars(ps)
    if (co2ice_sfc > 0) then
      ! Case 1: CO2 ice already on the ground.
      ! Add condensation to existing CO2 ice mass.
      co2ice_sfc = co2ice_sfc + dmsdt * dt
      tg = tsat
    else
      ! Case 2: No CO2 ice on the ground; Ground is warmer.
      ! Ground temperature drops when CO2 ice sublimes on warmer surface.
      tinp = sqrdy / zin(1)
      tgp  = tg - tinp * dt * dmsdt * xlhtc
      ! Calculate how much CO2 ice sublimes on hitting the ground.
      if (tgp >= tsat) then
        ! Case 2A: All CO2 ice sublimes; No net condensation.
        dmsdt      = 0
        tg         = tgp
        co2ice_sfc = 0
      else
        ! Case 2B: Ground cooled to CO2 frost point and some CO2 ice remains.
        ! Calculate how much CO2 ice remains.
        dmsdt      = dmsdt * (tsat - tgp) / (tg - tgp)
        tg         = tsat
        co2ice_sfc = dmsdt * dt
      end if
    end if
  end if

  ! Aerosol scavenging by CO2 snow fall.
  if (co2scav) then
    k_scavup = 0
    k_scavdn = 0
    ! Determine the portion of atmosphere affected by CO2 condensation.
    do k = 1, nlev
      if (dm(k) > 0) then
        k_scavup = k
        exit
      end if
    end do
    do k = nlev, 1, -1
      if (dm(k) > 0) then
        k_scavdn = k
        exit
      end if
    end do
    if (k_scavup /= 0 .and. k_scavdn /= 0) then
      do k = k_scavup, k_scavdn
        if (k_scavdn == nlev) then
          ! If condensation occurs down to the surface, put all aerosols on the surface (in fact, only the cloud mass matters).
          tm_sfc      (iMa_vap) = tm_sfc      (iMa_vap) + scaveff * q(k,iMa_cld) * dp(k) / g
          tm_sfc      (iMa_dst) = tm_sfc      (iMa_dst) + scaveff * q(k,iMa_dst) * dp(k) / g
          tm_sfc      (iMa_cor) = tm_sfc      (iMa_cor) + scaveff * q(k,iMa_cor) * dp(k) / g
          tmflx_sfc_dn(iMa_dst) = tmflx_sfc_dn(iMa_dst) + scaveff * q(k,iMa_dst) * dp(k) / g / dt
          tmflx_sfc_dn(iMa_cor) = tmflx_sfc_dn(iMa_cor) + scaveff * q(k,iMa_cor) * dp(k) / g / dt
          tmflx_sfc_dn(iMa_cld) = tmflx_sfc_dn(iMa_cld) + scaveff * q(k,iMa_cld) * dp(k) / g / dt
        else
          ! If condensation occurs in a restricted portion, put aerosols in the highest layer unaffected by CO2 condensation.
          do m = 1, ntracers
            q(k_scavdn+1,m) = q(k_scavdn+1,m) + scaveff * q(k,m) * dp(k) / dp(k_scavdn+1)
          end do
        end if
        do m = 1, naer
          q(k,m) = q(k,m) * (1 - scaveff)
        end do
        dm(k) = 0
      end do
    end if
  end if

end subroutine coldair