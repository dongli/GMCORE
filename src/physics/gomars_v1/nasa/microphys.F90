subroutine microphys( &
  ps                , &
  p                 , &
  dp_dry            , &
  dz                , &
  dz_lev            , &
  t                 , &
  t_lev             , &
  tg                , &
  q                 , &
  co2ice_sfc        , &
  taux              , &
  tauy              , &
  ht_pbl            , &
  ptop_pbl          , &
  kh                , &
  tm_sfc            , &
  dstflx_wsl        , &
  dstflx_ddl        , &
  rho               , &
  deposit           , &
  tmflx_sfc_dn      )

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use formula_mod
  use gomars_v1_const_mod
  use gomars_v1_tracers_mod
  use gomars_v1_mp_mod

  implicit none

  real(r8), intent(in   ) :: ps
  real(r8), intent(in   ) :: p           (nlev)
  real(r8), intent(in   ) :: dp_dry      (nlev)
  real(r8), intent(in   ) :: dz          (nlev)
  real(r8), intent(in   ) :: dz_lev      (nlev+1)
  real(r8), intent(in   ) :: t           (nlev)
  real(r8), intent(in   ) :: t_lev       (nlev+1)
  real(r8), intent(in   ) :: tg
  real(r8), intent(in   ) :: q           (nlev,ntracers)
  real(r8), intent(in   ) :: co2ice_sfc
  real(r8), intent(in   ) :: taux
  real(r8), intent(in   ) :: tauy
  real(r8), intent(in   ) :: ht_pbl
  real(r8), intent(in   ) :: ptop_pbl
  real(r8), intent(in   ) :: kh          (nlev+1)
  real(r8), intent(inout) :: tm_sfc      (ntracers)
  real(r8), intent(  out) :: dstflx_wsl
  real(r8), intent(  out) :: dstflx_ddl
  real(r8), intent(  out) :: rho         (nlev)
  real(r8), intent(  out) :: deposit     (ntracers)
  real(r8), intent(  out) :: tmflx_sfc_dn(ntracers)

  integer k, m
  real(r8) Mo
  real(r8) No
  real(r8) ro  (nlev,ntracers)
  real(r8) dens(nlev,ntracers)
  real(r8) rho_lev(nlev+1)

  ! Inject dust into the atmosphere.
  call dust_update(ps, p, dp_dry, t, tg, taux, tauy, &
    ht_pbl, ptop_pbl, co2ice_sfc, dstflx_wsl, dstflx_ddl, q, tm_sfc)

  rho = dry_air_density(t, p)

  do m = 1, ntracers
    select case (m)
    case (1)
      do k = 1, nlev
        Mo = q(k,iMa_dst)
        No = q(k,iNb_dst) + eps ! FIXME: Original code uses 1.0e-50_r8.
        dens(k,m) = rho_aer(m)
        ro  (k,m) = (Mo / No * 0.75_r8 / pi / dens(k,m))**athird * exp(3 * std_aer(m)**2)
        if (Mo < 1.0e-20_r8) ro(k,m) = 1.0e-8_r8
      end do
    case (2)
      do k = 1, nlev
        dens(k,m) = rho_aer(iMa_dst)
        ro  (k,m) = ro(k,iMa_dst) * exp(-3 * std_aer(m)**2)
      end do
    case (3)
      do k = 1, nlev
        Mo = q(k,iMa_cld)
        No = q(k,iNb_cld) + eps ! FIXME: Original code uses 1.0e-50_r8.
        dens(k,m) = q(k,iMa_cld) / Mo * rho_aer(iMa_cld) + &
                    q(k,iMa_cor) / Mo * rho_aer(iMa_dst)
        dens(k,m) = min(max(dens(k,m), rho_ice), rho_dst)
        ro  (k,m) = (Mo / No * 0.75_r8 / pi / dens(k,m))**athird * exp(3 * std_aer(m)**2)
        if (Mo < 1.0e-20_r8) ro(k,m) = 1.0e-8_r8
      end do
    case (4)
      do k = 1, nlev
        dens(k,m) = dens(k,iMa_cld)
        ro  (k,m) = ro  (k,iMa_cld) * exp(-3 * std_aer(m)**2)
      end do
    case (5)
      do k = 1, nlev
        dens(k,m) = dens(k,iMa_cld)
        ro  (k,m) = ro  (k,iMa_cld)
      end do
    end select
  end do

  ! Now compute the microphysical processes.
  ! Always check the order of sedimentation and nucleation condensation.
  ! PBL is called before microphysics, so do sedimentation first.
  call sedim(p, dp_dry, t_lev, rho_lev, dz, dz_lev, kh, q, ro, dens, deposit, tmflx_sfc_dn)

end subroutine microphys

subroutine dust_update(ps, p, dp_dry, t, tg, taux, tauy, &
  ht_pbl, ptop_pbl, co2ice_sfc, dstflx_wsl, dstflx_ddl, q, tm_sfc)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_namelist_mod
  use gomars_v1_tracers_mod

  implicit none

  real(r8), intent(in   ) :: ps
  real(r8), intent(in   ) :: p          (nlev)
  real(r8), intent(in   ) :: dp_dry     (nlev)
  real(r8), intent(in   ) :: t          (nlev)
  real(r8), intent(in   ) :: tg
  real(r8), intent(in   ) :: taux
  real(r8), intent(in   ) :: tauy
  real(r8), intent(in   ) :: ht_pbl
  real(r8), intent(in   ) :: ptop_pbl
  real(r8), intent(in   ) :: co2ice_sfc
  real(r8), intent(  out) :: dstflx_wsl
  real(r8), intent(  out) :: dstflx_ddl
  real(r8), intent(inout) :: q          (nlev,ntracers)
  real(r8), intent(inout) :: tm_sfc     (     ntracers)

  integer k
  real(r8) rho
  real(r8) reff_lift
  real(r8) rm
  real(r8) mp
  real(r8) dm

  if (use_wsl_newman) then
    rho = ps / (rd * tg)
    call wsl_newman(rho, taux, tauy, co2ice_sfc, dstflx_wsl)
    dm = dstflx_wsl * dt
    if (dm > 0) then
      reff_lift = 2.0e-6_r8
      rm = reff_lift * exp(-0.5_r8 * dev_dst**2)
      mp = 4.0_r8 / 3.0_r8 * pi * rm**3 * rho_dst
      q(nlev,iMa_dst) = q(nlev,iMa_dst) + dm / dp_dry(nlev)
      q(nlev,iNb_dst) = q(nlev,iNb_dst) + dm / dp_dry(nlev) / mp
      tm_sfc(iMa_dst) = tm_sfc(iMa_dst) - dm
    end if
  else if (use_wsl_kmh) then

  end if

  if (use_ddl .and. ptop_pbl /= 0) then
    call ddl(ps, ht_pbl, ptop_pbl, dstflx_ddl)
    dm = dstflx_ddl * dt
    if (dm > 0) then
      reff_lift = 2.0e-6_r8
      rm = reff_lift * exp(-0.5_r8 * dev_dst**2)
      mp = 4.0_r8 / 3.0_r8 * pi * rm**3 * rho_dst
      do k = nlev, 1, -1
        if (p(k) < ptop_pbl) exit
        q(k,iMa_dst) = q(k,iMa_dst) + dm / dp_dry(nlev)
        q(k,iNb_dst) = q(k,iNb_dst) + dm / dp_dry(nlev) / mp
      end do
      tm_sfc(iMa_dst) = tm_sfc(iMa_dst) - dm
    end if
  end if

end subroutine dust_update