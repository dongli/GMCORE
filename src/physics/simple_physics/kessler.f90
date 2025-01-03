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

subroutine kessler(nz, pt, qv, qc, qr, rho, pk, dt, z, precl)

  use const_mod

  implicit none

  integer, intent(in   )                :: nz    ! Number of thermodynamic levels in the column
  real(8), intent(inout), dimension(nz) :: pt    ! Potential temperature (K)
  real(8), intent(inout), dimension(nz) :: qv    ! Water vapor mixing ratio (gm/gm)
  real(8), intent(inout), dimension(nz) :: qc    ! Cloud water mixing ratio (gm/gm)
  real(8), intent(inout), dimension(nz) :: qr    ! Rain  water mixing ratio (gm/gm)
  real(8), intent(in   ), dimension(nz) :: rho   ! Dry air density (not mean state as in KW) (kg/m^3)
  real(8), intent(in   ), dimension(nz) :: pk    ! Exner function (p/p0)**(R/cp)
  real(8), intent(in   )                :: dt    ! Time step (s)
  real(8), intent(in   ), dimension(nz) :: z     ! Heights of thermo. levels in the grid column (m)
  real(8), intent(  out)                :: precl ! Precipitation rate (m_water / s)

  real(8), parameter :: psl   = 1000.0d0  ! Pressure at sea level (hPa)
  real(8), parameter :: rhoqr = 1000.0d0  ! Density of liquid water (kg/m^3)

  real(8), dimension(nz) :: r
  real(8), dimension(nz) :: rhalf
  real(8), dimension(nz) :: velqr
  real(8), dimension(nz) :: sed
  real(8), dimension(nz) :: pc
  real(8) ern
  real(8) qrprod
  real(8) prod
  real(8) qvs
  real(8) dt_max
  real(8) dt0
  integer k, rainsplit, nt

  do k = 1, nz
    r    (k) = 0.001d0 * rho(k) ! g cm-3
    rhalf(k) = sqrt(rho(nz) / rho(k))
    pc   (k) = 3.8d0 / (pk(k)**(1.0d0 / rd_o_cpd) * psl) ! hPa
    qr   (k) = max(qr(k), 0.0d0)
    ! Liquid water terminal velocity (m/s) following KW eq. 2.15
    velqr(k) = 36.34d0 * (qr(k) * r(k))**0.1364d0 * rhalf(k)
  end do

  ! Maximum time step size in accordance with CFL condition
  dt_max = dt
  do k = nz, 2, -1
    if (abs(velqr(k)) > 1.0d-12) then
      dt_max = min(dt_max, 0.8d0 * (z(k-1) - z(k)) / velqr(k))
    end if
  end do

  ! Number of subcycles
  rainsplit = ceiling(dt / dt_max)
  if (rainsplit < 1) stop 'Kessler: rainsplit < 1'
  dt0 = dt / real(rainsplit, 8)

  ! Subcycle through rain process
  precl = 0

  do nt = 1, rainsplit
    ! Precipitation rate (m/s)
    precl = precl + rho(nz) * qr(nz) * velqr(nz) / rhoqr
    ! Sedimentation term using upstream differencing
    do k = nz, 2, -1
      sed(k) = dt0 * (r(k-1) * qr(k-1) * velqr(k-1) - r(k) * qr(k) * velqr(k)) / (r(k) * (z(k-1) - z(k)))
    end do
    sed(1) = -dt0 * qr(1) * velqr(1) / (0.5d0 * (z(1) - z(2)))

    ! Adjustment terms
    do k = nz, 1, -1
      ! Autoconversion and accretion rates following KW eq. 2.13a,b
      qrprod = qc(k) - (qc(k) - dt0 * max(0.001d0 * (qc(k) - 0.001d0), 0.0d0)) / (1 + dt0 * 2.2d0 * qr(k)**0.875d0)
      qc(k) = max(qc(k) - qrprod, 0.0d0)
      qr(k) = max(qr(k) + qrprod + sed(k), 0.0d0)

      ! Saturation vapor mixing ratio (gm/gm) following KW eq. 2.11
      qvs = pc(k) * exp(17.27d0 * (pk(k) * pt(k) - 273) / (pk(k) * pt(k)- 36))
      prod = (qv(k) - qvs) / (1 + qvs * (4093 * lv / cpd) / (pk(k) * pt(k) - 36)**2)

      ! Evaporation rate following KW eq. 2.14a,b
      ern = min(dt0 * (((1.6d0 + 124.9d0 * (r(k) * qr(k))**0.2046d0) &
            * (r(k) * qr(k))**0.525d0) / (2.55d6 * pc(k)             &
            / (3.8d0 * qvs) + 5.4d5)) * (dim(qvs, qv(k))             &
            / (r(k) * qvs)), max(-prod - qc(k), 0.0d0), qr(k))

      ! Saturation adjustment following KW eq. 3.10
      pt(k) = pt(k) + lv / (cpd * pk(k)) * (max(prod, -qc(k)) - ern)
      qv(k) = max(qv(k) - max(prod, -qc(k)) + ern, 0.0d0)
      qc(k) = qc(k) + max(prod, -qc(k))
      qr(k) = qr(k) - ern
    end do
  end do

  precl = precl / dble(rainsplit)

end subroutine kessler
