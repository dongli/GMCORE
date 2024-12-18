subroutine sedim(p, t_lev, rho, rho_lev, dz, kh, q, ro, dens, deposit)

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

  use gomars_v1_const_mod

  implicit none

  real(r8), intent(in   ) :: p          (nlev)
  real(r8), intent(in   ) :: t_lev      (nlev)
  real(r8), intent(in   ) :: rho        (nlev)
  real(r8), intent(in   ) :: rho_lev    (nlev+1)
  real(r8), intent(in   ) :: dz         (nlev)
  real(r8), intent(in   ) :: kh         (nlev+1)        ! Vertical eddy coefficient from PBL (m2 s-1)
  real(r8), intent(inout) :: q          (nlev,ntracers)
  real(r8), intent(inout) :: ro         (nlev,ntracers) ! Particle radius (m)
  real(r8), intent(inout) :: dens       (nlev,ntracers) ! Particle density (kg m-3)
  real(r8), intent(inout) :: deposit    (     ntracers)

  integer k, m
  real(r8) vt       ! Terminal velocity (m s-1)
  real(r8) dv       ! Gas molecule mean diffusion velocity (m s-1)
  real(r8) mfp      ! Gas molecule mean free path (m)
  real(r8) kn       ! Knudsen number (<<1: continuum regime, ~1: transition regime >>1: free molecular regime)
  real(r8) cc       ! Cunningham slip correction factor
  real(r8) vf

  deposit = 0

  do m = 1, ntracers
    do k = 2, nlev - 1
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
      vf = 2 * g * ro(k,m)**2 * dens(k,m) / (9 * dv) * (1 + cc * kn)
      vf = -vf * exp(-std_aer(m)**2)
      ! Correct the fall velocity accounting for mixing.
      if (kh(k) /= 0) then

      end if
    end do
  end do

end subroutine sedim