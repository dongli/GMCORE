subroutine tempgr( &
  lat            , &
  ps             , &
  tbot           , &
  qbot           , &
  tm_sfc         , &
  alsp           , &
  polarcap       , &
  rhouch         , &
  co2ice_sfc     , &
  ht_sfc         , &
  irflx_sfc_dn   , &
  vsflx_sfc_dn   , &
  h2osub_sfc     , &
  h2oice_sfc     , &
  rhosoil        , &
  cpsoil         , &
  scond          , &
  stemp          , &
  zin            , &
  tg             , &
  als            )

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  ! Calculate the ground temperature. Where CO2 ice is present on the ground, this also involves
  ! calculation of the mass of CO2 ice on the ground and surface pressure.

  use formula_mod
  use gomars_v1_const_mod
  use gomars_v1_namelist_mod
  use gomars_v1_tracers_mod
  use gomars_v1_lsm_mod

  implicit none

  real(r8), intent(in   ) :: lat
  real(r8), intent(in   ) :: ps
  real(r8), intent(in   ) :: tbot
  real(r8), intent(in   ) :: qbot     (ntracers)
  real(r8), intent(in   ) :: tm_sfc   (ntracers)
  real(r8), intent(in   ) :: alsp
  logical , intent(in   ) :: polarcap
  real(r8), intent(in   ) :: rhouch
  real(r8), intent(inout) :: co2ice_sfc
  real(r8), intent(in   ) :: ht_sfc
  real(r8), intent(in   ) :: irflx_sfc_dn
  real(r8), intent(in   ) :: vsflx_sfc_dn
  real(r8), intent(  out) :: h2osub_sfc
  real(r8), intent(  out) :: h2oice_sfc
  real(r8), intent(in   ) :: rhosoil  (nsoil)
  real(r8), intent(in   ) :: cpsoil   (nsoil)
  real(r8), intent(in   ) :: scond    (2*nsoil+1)
  real(r8), intent(inout) :: stemp    (2*nsoil+1)
  real(r8), intent(in   ) :: zin      (nsoil)
  real(r8), intent(  out) :: tg
  real(r8), intent(  out) :: als

  logical, save :: is_first_call = .true.
  integer k, l
  real(r8) dmgdt
  real(r8) tsat
  real(r8) emg15
  real(r8) emgout
  real(r8) downir
  real(r8) rhoucht
  real(r8) fcdn
  real(r8) tgp
  real(r8) tinp
  real(r8) wflux
  real(r8) qsat
  real(r8) flux(2*nsoil+1)

  if (is_first_call) then
    h2osub_sfc = 0
    ! FIXME: Could surface water vapor be considered as water ice?
    h2oice_sfc = tm_sfc(iMa_vap)
    is_first_call = .false.
  end if

  dmgdt = 0

  ! Set surface albedo.
  als  = alsp
  if (co2ice_sfc > 0) then
    als = merge(alices, alicen, lat < 0)
  else if (albfeed .and. tm_sfc(iMa_vap) > icethresh_kgm2 .and. .not. polarcap) then
    als = icealb
  end if

  tsat = 3182.48_r8 / (23.3494_r8 - log(ps / 100.0_r8))

  if (co2ice_sfc <= 0) then
    co2ice_sfc  = 0
    emg15   = eg15gnd
    emgout  = egognd
    downir  = emg15 * irflx_sfc_dn
    rhoucht = rhouch * tbot

    call newtg(      &
      als          , &
      vsflx_sfc_dn , &
      downir       , &
      rhouch       , &
      rhoucht      , &
      scond        , &
      stemp        , &
      sthick       , &
      ps           , &
      qbot(iMa_vap), &
      h2oice_sfc   , &
      h2osub_sfc   , &
      polarcap     , &
      tg             &
    )

    if (tg < tsat) then
      tg = tsat

      emg15  = eg15gnd
      emgout = egognd

      fcdn = -2 * scond(2) * (stemp(2) - tsat) / sthick(2)
      tgp = dt * ((1 - als) * vsflx_sfc_dn + ht_sfc - emg15 * (stbo * tsat**4) - fcdn) / xlhtc

      ! Check if there is any CO2 ice accumulation.
      if (tgp < 0) then
        dmgdt  = -tgp / dt
        co2ice_sfc = -tgp
      else
        dmgdt  = 0
        ! This term represents the last amounts of ice evaporating resulting in an
        ! increase in Tg. It still depends on TINP until I can figure out what to do with it.
        tinp   = sqrdy / zin(1)
        tg     = tsat + tgp * xlhtc * tinp
        co2ice_sfc = 0
      end if
    end if
  else ! CO2 ice on the ground
    tg = tsat

    ! Only modify if we have water condensation.
    qsat = water_vapor_saturation_mixing_ratio_mars(tg, ps)
    ! See Eq. (1) in Haberle et al. (2019), which used a bulk transfer approach,
    ! but note rhouch contains cpd, so here divides cpd.
    wflux = -rhouch * (qbot(iMa_vap) - qsat) / cpd
    if (wflux < 0) then
      h2oice_sfc = h2oice_sfc - wflux * dt
      h2osub_sfc = h2osub_sfc + wflux * dt
    end if

    if (.not. polarcap .and. h2osub_sfc > tm_sfc(iMa_vap)) then
      h2osub_sfc = tm_sfc(iMa_vap)
    end if

    if (lat < 0) then
      emg15  = eg15co2s
      emgout = egoco2s
    else
      emg15  = eg15co2n
      emgout = egoco2n
    end if

    ! New soil scheme: surface boundary condition with ice on the ground.
    fcdn = -2 * scond(2) * (stemp(2) - tsat) / sthick(2)
    tgp  = -co2ice_sfc + dt * ((1 - als) * vsflx_sfc_dn + ht_sfc - emg15 * (stbo * tsat**4) - fcdn) / xlhtc

    ! Check if there is still CO2 ice left.
    if (tgp < 0) then
      dmgdt      = -(co2ice_sfc + tgp) / dt
      co2ice_sfc = -tgp
    else
      dmgdt      = -co2ice_sfc / dt
      tinp       = sqrdy / zin(1)
      tg         = tsat + tgp * xlhtc * tinp
      co2ice_sfc = 0
    end if
  end if

  ! Calculate fluxes at layer boundaries (positive downward).
  do l = 3, 2 * nsoil - 1, 2
    flux(l) = -scond(l) * (stemp(l+1) - stemp(l-1)) / sthick(l)
  end do
  ! Calculate flux at top and bottom boundaries.
  flux(1        ) = -2 * scond(2) * (stemp(2) - tg) / sthick(2)
  flux(2*nsoil+1) = 0
  ! Update soil temperatures.
  do k = 1, nsoil
    l = 2 * k
    stemp(l) = stemp(l) - dt * (flux(l+1) - flux(l-1)) / (rhosoil(k) * cpsoil(k) * sthick(l))
  end do

  ! Calculate the CO2 condensation temperature at the surface.
  tsat = dewpoint_temperature_mars(ps)
  if (tg < tsat .or. co2ice_sfc > 0) then
    tg = tsat
  end if

end subroutine tempgr