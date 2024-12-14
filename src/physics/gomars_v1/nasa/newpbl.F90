subroutine newpbl( &
  z0, tg, ht_rad, ps, ts, polarcap, u, v, pt, pt_lev, q, dp_dry, &
  z, dz, z_lev, dz_lev, shr2, ri, km, kh, ustar, tstar, &
  taux, tauy, ht_pbl, rhouch, tm_sfc, h2osub_sfc)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use formula_mod
  use gomars_v1_const_mod
  use gomars_v1_tracers_mod
  use gomars_v1_pbl_mod

  implicit none

  real(r8), intent(in   ) :: z0                       ! Roughness length (m)
  real(r8), intent(in   ) :: tg                       ! Ground temperature (K)
  real(r8), intent(in   ) :: ht_rad(0:nlev  )         ! Radiative heating rate on full levels including TOA (K s-1)
  real(r8), intent(in   ) :: ps                       ! Surface pressure (Pa)
  real(r8), intent(in   ) :: ts                       ! Surface air temperature (K)
  logical , intent(in   ) :: polarcap
  real(r8), intent(inout) :: u     (  nlev  )
  real(r8), intent(inout) :: v     (  nlev  )
  real(r8), intent(inout) :: pt    (  nlev  )
  real(r8), intent(in   ) :: pt_lev(  nlev+1)
  real(r8), intent(inout) :: q     (  nlev,ntracers)
  real(r8), intent(in   ) :: dp_dry(  nlev  )         ! Dry-air weight (Pa)
  real(r8), intent(in   ) :: z     (  nlev  )
  real(r8), intent(in   ) :: dz    (  nlev  )
  real(r8), intent(in   ) :: z_lev (  nlev+1)
  real(r8), intent(in   ) :: dz_lev(  nlev+1)
  real(r8), intent(inout) :: shr2  (  nlev+1)
  real(r8), intent(  out) :: ri    (  nlev+1)
  real(r8), intent(  out) :: km    (  nlev+1)
  real(r8), intent(  out) :: kh    (  nlev+1)
  real(r8), intent(  out) :: ustar
  real(r8), intent(  out) :: tstar
  real(r8), intent(  out) :: taux                     ! Wind stress in x direction (N m-2)
  real(r8), intent(  out) :: tauy                     ! Wind stress in y direction (N m-2)
  real(r8), intent(  out) :: ht_pbl                   ! Heat rate at the surface (K s-1)
  real(r8), intent(  out) :: rhouch
  real(r8), intent(inout) :: tm_sfc(ntracers)         ! Tracer mass on the ground (kg)
  real(r8), intent(inout) :: h2osub_sfc               ! Water ice upward sublimation flux at the surface (kg m-2 s-1)

  integer i, n
  real(r8) lnzz
  real(r8) rhos
  real(r8) cdm
  real(r8) cdh
  real(r8) alpha
  real(r8) qsat
  real(r8) coef
  real(r8) kdf(nlev+1,nvar)
  real(r8) var(nlev  ,nvar)
  real(r8) bnd(       nvar)
  real(r8) rhs(nlev  ,nvar)

  n = nlev

  lnzz = log(z(n) / z0)

  call eddycoef(z_lev, dz_lev, u, v, pt, q, pt_lev, shr2, ri, km, kh)
  call bndcond(u, v, pt, tg, z, lnzz, z0, cdm, cdh, ustar, tstar)

  ! Reduce the sublimation flux by a coefficient. It is a tunable parameter to
  ! avoid the formation of low-lying clouds in summer above the north permanent cap.
  qsat   = water_vapor_saturation_mixing_ratio_mars(tg, ps)
  coef = merge(1.0_r8, 1.0_r8, qsat > q(n,iMa_vap))

  i = 1
  ! Zonal wind
  i = i + 1
  kdf(:,i) = km
  bnd(  i) = ustar * dt / dz(n) * sqrt(cdm)
  rhs(:,i) = 0
  var(:,i) = dp_dry * u
  ! Meridional wind
  i = i + 1
  kdf(:,i) = km
  bnd(  i) = bnd(1)
  rhs(:,i) = 0
  var(:,i) = dp_dry * v
  ! Potential temperature
  i = i + 1
  kdf(:,i) = kh
  bnd(  i) = ustar * dt / dz(n) * cdh
  rhs(:,i) = ht_rad * dt
  rhs(n,i) = rhs(n,i) + bnd(i) * dp_dry(n) * tg
  var(:,i) = dp_dry * pt
  ! Water vapor
  i = i + 1
  kdf(:,i) = kh
  bnd(  i) = ustar * dt / dz(n) * cdh * coef
  rhs(:,i) = 0
  rhs(n,i) = h2osub_sfc
  var(:,i) = dp_dry * q(:,iMa_vap)

  do i = 1, nvar
    call pbl_solve(kdf(:,i), bnd(i), rhs(:,i), var(:,i))
  end do

  u            = var(:,1) / dp_dry
  v            = var(:,2) / dp_dry
  pt           = var(:,3) / dp_dry
  q(:,iMa_vap) = var(:,4) / dp_dry

  ! Diagnose the surface stress and heat fluxes.
  rhos = dry_air_density(ps, ts)
  alpha  = atan2(v(n), u(n))
  taux   =  rhos * ustar**2 * cos(alpha)
  tauy   =  rhos * ustar**2 * sin(alpha)
  ht_pbl = -rhos * cpd * ustar * tstar
  rhouch =  rhos * cpd * ustar * cdh

  ! Update water ice budget on the surface.
  tm_sfc(iMa_vap) = tm_sfc(iMa_vap) - h2osub_sfc
  if (.not. polarcap .and. tm_sfc(iMa_vap) < 0) then
    tm_sfc(iMa_vap) = 0
  end if

contains

  subroutine pbl_solve(kdf, bnd, rhs, var)

    use math_mod

    real(r8), intent(in   ) :: kdf(nlev+1)
    real(r8), intent(in   ) :: bnd         ! Lower boundary condition for coefficient c
    real(r8), intent(in   ) :: rhs(nlev  ) ! Right hand side
    real(r8), intent(inout) :: var(nlev  )

    integer k
    real(r8) a(nlev), b(nlev), c(nlev), d(nlev)

    ! Setup the tridiagonal matrix.
    a(1) = 0
    do k = 2, nlev
      a(k) = -dt * kdf(k  ) / dz(k) / dz_lev(k  )
    end do
    do k = 1, nlev - 1
      c(k) = -dt * kdf(k+1) / dz(k) / dz_lev(k+1)
    end do
    c(nlev) = bnd
    do k = 1, nlev
      b(k) = 1 - a(k) - c(k)
    end do
    d = var + rhs

    call tridiag_thomas(a, b, c, d, var)

  end subroutine pbl_solve

end subroutine newpbl