subroutine eddycoef(z_lev, dz_lev, psi, shr2, ri, km, kh)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_pbl_mod

  implicit none

  real(r8), intent(in   ) :: z_lev (nlev+1)
  real(r8), intent(in   ) :: dz_lev(nlev+1)
  real(r8), intent(in   ) :: psi   (2*nlev+1,nvar)
  real(r8), intent(inout) :: shr2  (nlev+1)
  real(r8), intent(  out) :: ri    (nlev+1)  
  real(r8), intent(  out) :: km    (nlev+1)
  real(r8), intent(  out) :: kh    (nlev+1)

  real(r8), parameter :: Sm     = 0.393_r8
  real(r8), parameter :: Sh     = 0.493_r8
  real(r8), parameter :: sqrtGM = sqrt(0.153_r8)

  integer k, l
  real(r8) ml2
  real(r8) dudz
  real(r8) dvdz
  real(r8) dptdz
  real(r8) shr
  real(r8) km0
  real(r8) kh0
  real(r8) kmin

  do l = 3, 2 * nlev - 1, 2 ! Loop on half levels excluding top and bottom.
    k = (l + 1) / 2 ! Half level index
    ! Calculate mixing length and beta (volume expansion coefficient).
    ml2 = (ml0 * ka * z_lev(k) / (ml0 + ka * z_lev(k)))**2
    ! Calculate gradient Richardson number.
    dudz    = (psi(l+1,1) - psi(l-1,1)) / dz_lev(k)
    dvdz    = (psi(l+1,2) - psi(l-1,2)) / dz_lev(k)
    dptdz   = (psi(l+1,3) - psi(l-1,3)) / dz_lev(k)
    ! Smooth the wind shear,
    shr2(k) = shr2(k) - (shr2(k) - dudz**2 - dvdz**2) * dt / 1.0e4_r8
    shr    = sqrt(shr2(k))
    ! NOTE: The potential temperature on half levels is already available in psi.
    ri  (k) = g / psi(l,3) * dptdz / (shr2(k) + 1.0e-9_r8)
    ! Calculate neutral eddy coefficients.
    km0 = Sm * ml2 * shr / sqrtGM
    kh0 = Sh * ml2 * shr / sqrtGM
    ! Calculate eddy mixing coefficients.
    if (ri(k) <= 0) then
      km(k) = km0 * (1 - 15 * ri(k))**0.25_r8
      kh(k) = kh0 * (1 - 15 * ri(k))**0.50_r8
    else
      km(k) = km0 * (1 - ri(k) / ric)
      kh(k) = kh0 * (1 - ri(k) / ric)
    end if
  end do
  ! Limit the coefficients.
  kmin = merge(0.1_r8, 0.001_r8, z_lev(k) < 300)
  km(k) = max(km(k), kmin)
  kh(k) = max(kh(k), kmin)

end subroutine eddycoef