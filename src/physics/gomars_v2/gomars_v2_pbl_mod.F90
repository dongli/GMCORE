module gomars_v2_pbl_mod

  use gomars_v2_const_mod
  use gomars_v2_namelist_mod
  use gomars_v2_types_mod
  use gomars_v2_objects_mod
  use gomars_v2_sfc_mod

  implicit none

  private

  public gomars_v2_pbl_init
  public gomars_v2_pbl_final
  public gomars_v2_pbl_run

contains

  subroutine gomars_v2_pbl_init()

    integer iblk

    do iblk = 1, size(objects)
      associate (state => objects(iblk)%state)
      state%z0   = 0.01_r8
      state%shr2 = 0
      end associate
    end do

  end subroutine gomars_v2_pbl_init

  subroutine gomars_v2_pbl_final()

  end subroutine gomars_v2_pbl_final

  subroutine gomars_v2_pbl_run(state, tend, dt)

    type(gomars_v2_state_type), intent(inout) :: state
    type(gomars_v2_tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    real(r8) alpha, rhos, c1, c2
    real(r8) c(state%mesh%nlev+1)
    real(r8) A(state%mesh%nlev,3)
    real(r8) b(state%mesh%nlev)
    integer nlev, icol, k

    call calc_eddy_coef(state, dt)
    call gomars_v2_sfc_run(state)

    associate (mesh   => state%mesh   , &
               beta   => pbl_beta     , &
               dz     => state%dz     , & ! in
               u      => state%u      , & ! in
               v      => state%v      , & ! in
               ps     => state%ps     , & ! in
               ts     => state%ts     , & ! in
               rho    => state%rhod   , & ! in
               ustar  => state%ustar  , & ! in
               tstar  => state%tstar  , & ! in
               cdh    => state%cdh    , & ! in
               km     => state%eddy_km, & ! in
               kh     => state%eddy_kh, & ! in
               taux   => state%taux   , & ! out
               tauy   => state%tauy   , & ! out
               hflx   => state%hflx   , & ! out
               rhouch => state%rhouch )   ! out
    nlev = mesh%nlev
    do icol = 1, mesh%ncol
      alpha = atan2(v(icol,nlev), u(icol,nlev))
      ! Calculate surface air density.
      rhos = ps(icol) / rd / ts(icol)
      taux  (icol) =  rhos * ustar(icol) * ustar(icol) * cos(alpha)
      tauy  (icol) =  rhos * ustar(icol) * ustar(icol) * sin(alpha)
      hflx  (icol) = -rhos * ustar(icol) * tstar(icol) * cpd
      rhouch(icol) =  rhos * ustar(icol) * cdh  (icol) * cpd
      ! U
      do k = 2, mesh%nlev ! Interface levels excluding top and bottom
        c(k) = 2 * dt * km(icol,k) / (dz(icol,k-1) + dz(icol,k))
      end do
      do k = 1, mesh%nlev
        A(k,1) = -beta / dz(icol,k) * c(k  )
        A(k,3) = -beta / dz(icol,k) * c(k+1)
        A(k,2) = 1 - A(k,1) - A(k,3)
        c1 = (1 - beta) / dz(icol,k) * c(k  )
        c2 = (1 - beta) / dz(icol,k) * c(k+1)
        b(k) = c1 * u(icol,k-1) + (1 - c1 - c2) * u(icol,k) + c2 * u(icol,k+1)
      end do
    end do
    end associate

  end subroutine gomars_v2_pbl_run

  subroutine calc_eddy_coef(state, dt)

    type(gomars_v2_state_type), intent(inout) :: state
    real(r8), intent(in) :: dt

    real(r8), parameter :: Sm = 0.393_r8
    real(r8), parameter :: Sh = 0.493_r8
    real(r8), parameter :: sqrtGM = sqrt(0.153_r8)
    real(r8) rl2, beta, dptdz, dudz, dvdz, shr, ri, km0, kh0, kmin
    integer icol, k

    associate (mesh   => state%mesh   , &
               z      => state%z      , & ! in
               z_lev  => state%z_lev  , & ! in
               dz     => state%dz     , & ! in
               dz_lev => state%dz_lev , & ! in
               pt     => state%pt     , & ! in
               u      => state%u      , & ! in
               v      => state%v      , & ! in
               shr2   => state%shr2   , & ! inout
               km     => state%eddy_km, & ! out
               kh     => state%eddy_kh)   ! out
    do k = 2, mesh%nlev ! Loop on interfaces excluding top and bottom.
      do icol = 1, mesh%ncol
        ! Calculate mixing length and beta (volume expansion coefficient).
        rl2 = (rl0 * ka * z_lev(icol,k) / (rl0 + ka * z_lev(icol,k)))**2
        beta = (dz(icol,k-1) + dz(icol,k)) / (dz(icol,k-1) * pt(icol,k-1) + dz(icol,k) * pt(icol,k))
        ! Calculate gradient Richardson number.
        dudz = (u(icol,k) - u(icol,k-1)) / dz_lev(icol,k)
        dvdz = (v(icol,k) - v(icol,k-1)) / dz_lev(icol,k)
        dptdz = (pt(icol,k) - pt(icol,k-1)) / dz_lev(icol,k)
        ! Smooth the wind shear used to calculate the gradient Richardson number.
        shr2(icol,k) = shr2(icol,k) - (shr2(icol,k) - dudz**2 - dvdz**2) * dt / 1.0e4_r8
        shr = sqrt(shr2(icol,k))
        ri = beta * g * dptdz / (shr2(icol,k) + 1.0e-9_r8)
        ! Calculate neutral eddy coefficients.
        km0 = Sm * rl2 * shr / sqrtGM
        kh0 = Sh * rl2 * shr / sqrtGM
        ! Calculate eddy mixing coefficients.
        if (ri <= 0) then
          km(icol,k) = km0 * (1 - 15 * ri)**0.25_r8
          kh(icol,k) = kh0 * (1 - 15 * ri)**0.50_r8
        else
          km(icol,k) = km0 * (1 - ri / ric)
          kh(icol,k) = kh0 * (1 - ri / ric)
        end if
        ! Limit the coefficients.
        kmin = merge(0.1_r8, 0.001_r8, z_lev(icol,k) < 300)
        km(icol,k) = max(km(icol,k), kmin)
        kh(icol,k) = max(kh(icol,k), kmin)
      end do
    end do
    end associate

  end subroutine calc_eddy_coef

end module gomars_v2_pbl_mod
