!-----------------------------------------------------------------------
!
!  Version:  2.0
!
!  Date:  January 22nd, 2015
!
!  Change log:
!  v2 - Added sub-cycling of rain sedimentation so as not to violate
!       CFL condition.
!
!  The KESSLER subroutine implements the Kessler (1969) microphysics
!  parameterization as described by Soong and Ogura (1973) and Klemp
!  and Wilhelmson (1978, KW). KESSLER is called at the end of each
!  time step and makes the final adjustments to the potential
!  temperature and moisture variables due to microphysical processes
!  occurring during that time step. KESSLER is called once for each
!  vertical column of grid cells. Increments are computed and added
!  into the respective variables. The Kessler scheme contains three
!  moisture categories: water vapor, cloud water (liquid water that
!  moves with the flow), and rain water (liquid water that falls
!  relative to the surrounding air). There  are no ice categories.
!  Variables in the column are ordered from the surface to the top.
!
!  SUBROUTINE KESSLER(theta, qv, qc, qr, rho, pk, dt, z, nz, rainnc)
!
!  Input variables:
!     theta  - potential temperature (K)
!     qv     - water vapor mixing ratio (gm/gm)
!     qc     - cloud water mixing ratio (gm/gm)
!     qr     - rain  water mixing ratio (gm/gm)
!     rho    - dry air density (not mean state as in KW) (kg/m^3)
!     pk     - Exner function  (not mean state as in KW) (p/p0)**(R/cp)
!     dt     - time step (s)
!     z      - heights of thermodynamic levels in the grid column (m)
!     nz     - number of thermodynamic levels in the column
!     precl  - Precipitation rate (m_water/s)
!
! Output variables:
!     Increments are added into t, qv, qc, qr, and rainnc which are
!     returned to the routine from which KESSLER was called. To obtain
!     the total precip qt, after calling the KESSLER routine, compute:
!
!       qt = sum over surface grid cells of (rainnc * cell area)  (kg)
!       [here, the conversion to kg uses (10^3 kg/m^3)*(10^-3 m/mm) = 1]
!
!
!  Authors: Paul Ullrich
!           University of California, Davis
!           Email: paullrich@ucdavis.edu
!
!           Based on a code by Joseph Klemp
!           (National Center for Atmospheric Research)
!
!  Reference:
!
!    Klemp, J. B., W. C. Skamarock, W. C., and S.-H. Park, 2015:
!    Idealized Global Nonhydrostatic Atmospheric Test Cases on a Reduced
!    Radius Sphere. Journal of Advances in Modeling Earth Systems. 
!    doi:10.1002/2015MS000435
!
!=======================================================================

subroutine kessler(nz, rd, cpd, theta, qv, qc, qr, rho, pk, dt, z, precl)

  implicit none

  integer, intent(in   )                :: nz    ! Number of thermodynamic levels in the column
  real(8), intent(in   )                :: rd    ! Gas constant for dry air (J/kg/K)
  real(8), intent(in   )                :: cpd   ! Specific heat capacity of dry air at constant pressure (J/kg/K)
  real(8), intent(inout), dimension(nz) :: theta ! Potential temperature (K)
  real(8), intent(inout), dimension(nz) :: qv    ! Water vapor mixing ratio (gm/gm)
  real(8), intent(inout), dimension(nz) :: qc    ! Cloud water mixing ratio (gm/gm)
  real(8), intent(inout), dimension(nz) :: qr    ! Rain  water mixing ratio (gm/gm)
  real(8), intent(in   ), dimension(nz) :: rho   ! Dry air density (not mean state as in KW) (kg/m^3)
  real(8), intent(  out)                :: precl ! Precipitation rate (m_water / s)
  real(8), intent(in   ), dimension(nz) :: z     ! Heights of thermo. levels in the grid column (m)
  real(8), intent(in   ), dimension(nz) :: pk    ! Exner function (p/p0)**(R/cp)
  real(8), intent(in   )                :: dt    ! Time step (s)

  real(8), parameter :: lv    = 2.5d6     ! Latent heat of vaporization (J/kg)
  real(8), parameter :: c1    = 17.27d0   ! Constant for saturation vapor pressure (dimensionless)
  real(8), parameter :: psl   = 1000.0d0  ! Pressure at sea level (hPa)
  real(8), parameter :: rhoqr = 1000.0d0  ! Density of liquid water (kg/m^3)

  real(8), dimension(nz) :: r
  real(8), dimension(nz) :: rhalf
  real(8), dimension(nz) :: velqr
  real(8), dimension(nz) :: sed
  real(8), dimension(nz) :: pc
  real(8) c2
  real(8) xk
  real(8) ern
  real(8) qrprod
  real(8) prod
  real(8) qvs
  real(8) dt_max
  real(8) dt0
  integer k, rainsplit, nt

  c2 = 237.3d0 * c1 * lv / cpd
  xk = rd / cpd

  do k = 1, nz
    r    (k) = 0.001d0 * rho(k)
    rhalf(k) = sqrt(rho(1) / rho(k))
    pc   (k) = 3.8d0 / (pk(k)**(1.0d0 / xk) * psl)
    ! Liquid water terminal velocity (m/s) following KW eq. 2.15
    velqr(k) = 36.34d0 * (qr(k) * r(k))**0.1364d0 * rhalf(k)
  end do

  ! Maximum time step size in accordance with CFL condition
  dt_max = dt
  do k = 1, nz - 1
    if (velqr(k) /= 0) then
      dt_max = min(dt_max, 0.8d0 * (z(k+1) - z(k)) / velqr(k))
    end if
  end do

  ! Number of subcycles
  rainsplit = ceiling(dt / dt_max)
  dt0 = dt / real(rainsplit, 8)

  ! Subcycle through rain process
  precl = 0

  do nt = 1, rainsplit
    ! Precipitation rate (m/s)
    precl = precl + rho(1) * qr(1) * velqr(1) / rhoqr
    ! Sedimentation term using upstream differencing
    do k = 1, nz - 1
      sed(k) = dt0 * (r(k+1) * qr(k+1) * velqr(k+1) - r(k) * qr(k) * velqr(k)) / (r(k) * (z(k+1) - z(k)))
    end do
    sed(nz) = -dt0 * qr(nz) * velqr(nz) / (0.5d0 * (z(nz) - z(nz-1)))

    ! Adjustment terms
    do k = 1, nz
      ! Autoconversion and accretion rates following KW eq. 2.13a,b
      qrprod = qc(k) - (qc(k) - dt0 * max(0.001d0 * (qc(k) - 0.001d0), 0.0d0)) / (1 + dt0 * 2.2d0 * qr(k)**0.875d0)
      qc(k) = max(qc(k) - qrprod, 0.0d0)
      qr(k) = max(qr(k) + qrprod + sed(k), 0.0d0)

      ! Saturation vapor mixing ratio (gm/gm) following KW eq. 2.11
      qvs = pc(k) * exp(c1 * (pk(k) * theta(k) - 273.0d0) / (pk(k) * theta(k)- 36.0d0))
      prod = (qv(k) - qvs) / (1 + qvs * c2 / (pk(k) * theta(k) - 36.0d0)**2)

      ! Evaporation rate following KW eq. 2.14a,b
      ern = min(dt0 * (((1.6d0 + 124.9d0 * (r(k) * qr(k))**0.2046d0)  &
            * (r(k) * qr(k))**0.525d0) / (2550000d0 * pc(k)           &
            / (3.8d0 * qvs) + 540000d0)) * (dim(qvs, qv(k))           &
            /(r(k) * qvs)), max(-prod - qc(k), 0.0d0), qr(k))

      ! Saturation adjustment following KW eq. 3.10
      theta(k) = theta(k) + 2500000d0 / (1003.0d0 * pk(k)) * (max(prod, -qc(k)) - ern)
      qv(k) = max(qv(k) - max(prod, -qc(k)) + ern, 0.0d0)
      qc(k) = qc(k) + max(prod, -qc(k))
      qr(k) = qr(k) - ern
    end do

    ! Recalculate liquid water terminal velocity
    if (nt /= rainsplit) then
      do k = 1, nz
        velqr(k)  = 36.34d0 * (qr(k) * r(k))**0.1364 * rhalf(k)
      end do
    end if
  end do

  precl = precl / dble(rainsplit)

end subroutine kessler
