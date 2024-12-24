subroutine sedim(p, dp_dry, t_lev, rho_lev, dz, dz_lev, kh, q, ro, dens, deposit, tmflx_sfc_dn)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  ! Sedimentation is based on the standard Stokes-Cunningham relationships for
  ! particle fall velocity with the slip correction for the thin Martian
  ! atmosphere.
  !
  !   vf = 2 * g * r^2 * ρ_p / (9 * μ) * (1 + cc * K_n)
  ! 
  ! where K_n is Knudsen number
  !
  !   K_n = λ / r
  !
  ! and λ is the mean free path of gas molecules.
  !
  !   cc = 1.246 + 0.42 * exp(-0.87 / K_n)

  use math_mod
  use gomars_v1_const_mod
  use gomars_v1_mp_mod
  use gomars_v1_tracers_mod

  implicit none

  real(r8), intent(in   ) :: p           (nlev)
  real(r8), intent(in   ) :: dp_dry      (nlev)
  real(r8), intent(in   ) :: t_lev       (nlev)
  real(r8), intent(in   ) :: rho_lev     (nlev+1)
  real(r8), intent(in   ) :: dz          (nlev)
  real(r8), intent(in   ) :: dz_lev      (nlev+1)
  real(r8), intent(in   ) :: kh          (nlev+1)        ! Vertical eddy coefficient from PBL (m2 s-1)
  real(r8), intent(inout) :: q           (nlev,ntracers)
  real(r8), intent(inout) :: ro          (nlev,ntracers) ! Particle radius (m)
  real(r8), intent(inout) :: dens        (nlev,ntracers) ! Particle density (kg m-3)
  real(r8), intent(inout) :: deposit     (     ntracers)
  real(r8), intent(inout) :: tmflx_sfc_dn(     ntracers) ! Tracer mass downward flux at the surface (kg m-2 s-1)

  integer k, m
  real(r8) q_old(nlev)
  real(r8) vt       ! Terminal velocity (m s-1)
  real(r8) dv       ! Gas molecule mean diffusion velocity (m s-1)
  real(r8) mfp      ! Gas molecule mean free path (m)
  real(r8) kn       ! Knudsen number (<<1: continuum regime, ~1: transition regime >>1: free molecular regime)
  real(r8) cc       ! Cunningham slip correction factor
  real(r8) w
  real(r8) lnpr     ! log(dp_try(k-1) / dp_try(k))
  real(r8) theta
  real(r8) sigma
  real(r8) rap
  real(r8) cfl
  real(r8) exp2t
  real(r8) ft(nlev+1)
  real(r8) fs(nlev+1)
  real(r8) as(nlev)
  real(r8) bs(nlev)
  real(r8) cs(nlev)
  real(r8) ds(nlev)
  logical solve_implicit

  deposit = 0

  do m = 1, ntracers
    q_old = dp_dry * q(:,m)
    q(:,m) = q_old
    do k = 2, nlev
      ! Calculate fall velocity on half levels.
      ! - Thermal velocity of CO2 molecules
      vt = sqrt(scale_co2 * t_lev(k))
      ! - CO2 gas diffusion velocity
      dv = 1.59e-6_r8 * t_lev(k)**1.5_r8 / (t_lev(k) + 244)
      ! - CO2 gas mean free path
      mfp = 2 * dv / (rho_lev(k) * vt)
      ! - Knudsen number
      kn = mfp / ro(k,m)
      ! - Cunningham slip correction factor (Alvarez et al., 2024 gave 1.168 + 0.552 * exp(-0.990 / kn)
      cc = 1.246_r8 + 0.42_r8 * exp(-0.87_r8 / kn)
      ! - Fall velocity
      w  = 2 * g * ro(k,m)**2 * dens(k,m) / (9 * dv) * (1 + cc * kn)
      w  = -w * exp(-std_aer(m)**2)
      ! - Courant number
      cfl = w * dt / dz_lev(k)
      ! Correct the fall velocity accounting for mixing.
      lnpr = log(dp_dry(k-1) / dp_dry(k))
      if (kh(k) /= 0) then
        theta = 0.5_r8 * (w * dz_lev(k) / kh(k) + lnpr)
        if (theta /= 0) then
          sigma = abs(1.0_r8 / tanh(theta) - 1.0_r8 / theta)
        else
          sigma = 1
        end if
      else
        sigma = 1
      end if
      if (q_old(k) == 0) then
        rap = 10
        if (q_old(k-1) == 0) then
          rap = 1
        end if
      else
        rap = min(max(q_old(k-1) / q_old(k), 0.1_r8), 10.0_r8)
      end if
      if (.not. (rap > 0.9 .and. rap < 1.1 .or. w * dt > dz(k))) then
        if (w < 0) then
          w = w * exp(sigma * log((exp(-cfl * log(rap)) - 1) / (cfl * (1 - rap))))
        else
          w = 0
        end if
      end if
      ! Calculate fluxes.
      if (kh(k) /= 0) then
        if (theta /= 0) then ! Has turbulent mixing and sedimentation or density gradient.
          exp2t = exp(2 * theta)
          ft(k) = (w + kh(k) / dz_lev(k) * lnpr) / (exp2t - 1)
          fs(k) = ft(k) * exp2t
        else ! Only has turbulent mixing and no density gradient.
          ft(k) = kh(k) / dz_lev(k)
          fs(k) = ft(k)
        end if
      else
        if (w < 0) then ! Only has sedimentation.
          ft(k) = -w
          fs(k) = 0
        else ! Neither has turbulent mixing nor sedimentation.
          ! FIXME: w should never be positive.
          ft(k) = 0
          fs(k) = w
        end if
      end if
    end do
    ! Boundary conditions for the fluxes.
    ft(1) = 0
    fs(1) = 0
    ft(nlev+1) = -w
    fs(nlev+1) = 0
    ! Calculate the coefficients of the continuity equation.
    solve_implicit = .false.
    do k = 1, nlev
      cs(k) = ft(k+1) + fs(k) - dz(k) / dt
      if (cs(k) > 0) then
        solve_implicit = .true.
        exit
      end if
      as(k) = -dz(k) / dt
      bs(k) = -ft(k)
      ds(k) = -fs(k+1)
    end do
    if (solve_implicit) then
      do k = 1, nlev
        as(k) = ft(k)
        bs(k) = -(ft(k+1) + fs(k) + dz(k) / dt)
        cs(k) = fs(k+1)
        ds(k) = -dz(k) / dt * q_old(k)
      end do
      call tridiag_thomas(as, bs, cs, ds, q(:,m))
    else
      q(1,m) = (cs(1) * q_old(1) + ds(1) * q_old(2)) / as(1)
      do k = 2, nlev - 1
        q(k,m) = (bs(k) * q_old(k-1) + cs(k) * q_old(k) + ds(k) * q_old(k+1)) / as(k)
      end do
      q(nlev,m) = (bs(nlev) * q_old(nlev-1) + cs(nlev) * q_old(nlev)) / as(nlev)
    end if
    ! Calculate the tracer mass falling on the ground.
    if (m == iMa_cld) then
      deposit     (iMa_vap) = deposit     (iMa_vap) + q(nlev,m) * ft(nlev+1) * dt
      tmflx_sfc_dn(iMa_vap) = tmflx_sfc_dn(iMa_vap) + q(nlev,m) * ft(nlev+1)
    else
      deposit     (m) = deposit     (m) + q(nlev,m) * ft(nlev+1) * dt
      tmflx_sfc_dn(m) = tmflx_sfc_dn(m) + q(nlev,m) * ft(nlev+1)
    end if
    q(:,m) = q(:,m) / dp_dry
  end do

end subroutine sedim