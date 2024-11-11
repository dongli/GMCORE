subroutine tempgr(lat, ps, tbot, q_vap, qcond_vap, alsp, npcflag, rhouch, &
                  co2ice, fa, dnirflux, dnvflux, subflux, gndice, rhosoil, cpsoil, scond, &
                  stemp, sthick, zin, gt, surfalb)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  ! Calculate the ground temperature. Where CO2 ice is present on the ground, this also involves
  ! calculation of the mass of CO2 ice on the ground and surface pressure.

  use gomars_v1_const_mod
  use gomars_v1_namelist_mod

  implicit none

  real(r8), intent(in) :: lat
  real(r8), intent(in) :: ps
  real(r8), intent(in) :: tbot
  real(r8), intent(in) :: q_vap
  real(r8), intent(in) :: qcond_vap
  real(r8), intent(in) :: alsp
  logical , intent(in) :: npcflag
  real(r8), intent(in) :: rhouch
  real(r8), intent(inout) :: co2ice
  real(r8), intent(in) :: fa
  real(r8), intent(in) :: dnirflux
  real(r8), intent(in) :: dnvflux
  real(r8), intent(out) :: subflux
  real(r8), intent(out) :: gndice
  real(r8), intent(in) :: rhosoil(nl)
  real(r8), intent(in) :: cpsoil(nl)
  real(r8), intent(in) :: scond(2*nl+1)
  real(r8), intent(inout) :: stemp(2*nl+1)
  real(r8), intent(in) :: sthick(2*nl+1)
  real(r8), intent(in) :: zin(nl)
  real(r8), intent(out) :: gt
  real(r8), intent(out) :: surfalb

  logical, save :: is_first_call = .true.
  real(r8) dmgdt
  real(r8) als
  real(r8) tsat
  real(r8) emg15
  real(r8) emgout
  real(r8) downir
  real(r8) rhoucht
  real(r8) fcdn
  real(r8) tgp
  real(r8) tinp
  real(r8) wflux
  real(r8) qgnd
  real(r8) flux(nl+1)
  integer k, l

  if (is_first_call) then
    subflux = 0
    gndice = qcond_vap
    is_first_call = .false.
  end if

  dmgdt = 0
  als = alsp

  if (co2ice > 0) then
    als = merge(alices, alicen, lat < 0)
  else if (albfeed .and. qcond_vap > icethresh_kgm2 .and. .not. npcflag) then
    als = icealb
  end if

  surfalb = als

  tsat = 3182.48_r8 / (23.3494_r8 - log(ps / 100.0_r8))

  if (co2ice <= 0) then
    co2ice = 0
    emg15  = eg15gnd
    emgout = egognd

    downir = emg15 * dnirflux
    rhoucht = rhouch * tbot

    ! Call newtg

    if (gt < tsat) then
      gt = tsat

      emg15 = eg15gnd
      emgout = egognd

      fcdn = -2 * scond(2) * (stemp(2) - tsat) / sthick(2)
      tgp = dt * ((1 - als) * dnvflux + fa - emg15 * (stbo * tsat**4) - fcdn) / xlhtc

      ! Check if there is any CO2 ice accumulation.
      if (tgp < 0) then
        dmgdt = -tgp / dt
        co2ice = -tgp
      else
        dmgdt = 0

        ! This term represents the last amounts of ice evaporating resulting in an
        ! increase in Tg. It still depends on TINP until I can figure out what to do with it.
        tinp = sqrdy / zin(1)
        gt = tsat + tgp * xlhtc * tinp
        co2ice = 0
      end if
    end if
  else ! CO2 ice on the ground
    gt = tsat

    ! Only modify if we have water condensation.
    qgnd = (18.0_r8 / 44.0_r8) * 6.11_r8 * exp(22.5_r8 * (1 - (273.16_r8 / gt))) / ps
    wflux = -rhouch * (q_vap - qgnd) / cpd
    if (wflux < 0) then
      gndice = gndice - wflux * dt
      subflux = subflux + wflux * dt
    end if
    if (lat < 0) then
      emg15 = eg15co2s
      emgout = egoco2s
    else
      emg15 = eg15co2n
      emgout = egoco2n
    end if

    ! New soil scheme: surface boundary condition with ice on the ground.
    fcdn = -2 * scond(2) * (stemp(2) - tsat) / sthick(2)
    tgp = -co2ice + dt * ((1 - als) * dnvflux + fa - emg15 * (stbo * tsat**4) - fcdn) / xlhtc

    ! Check if there is still CO2 ice left.
    if (tgp < 0) then
      dmgdt = -(co2ice + tgp) / dt
      co2ice = -tgp
    else
      dmgdt = -co2ice / dt

      tinp = sqrdy / zin(1)
      gt = tsat + tgp * xlhtc * tinp
      co2ice = 0
    end if
  end if

  ! Calculate fluxes at layer boundaries (positive downward).
  do l = 2, nl
    k = 2 * l - 1
    flux(l) = -scond(k) * (stemp(k+1) - stemp(k-1)) / sthick(k)
  end do
  ! Calculate flux at top and bottom boundaries.
  flux(1) = -2 * scond(2) * (stemp(2) - gt) / sthick(2)
  flux(nl+1) = 0
  ! Update soil temperatures.
  do l = 1, nl
    k = 2 * l
    stemp(k) = stemp(k) - dt * (flux(l+1) - flux(l)) / (rhosoil(l) * cpsoil(l) * sthick(k))
  end do

end subroutine tempgr