module tropical_cyclone_test_mod

! ======================================================================
!
!  Date:  July 29, 2015
!
!  Function for setting up idealized tropical cyclone initial conditions
!
!  SUBROUTINE tropical_cyclone_sample(
!    lon,lat,p,z,zcoords,ptv,rho,qv,u,v,t,phis)
!
!  Given a point specified by:
!      lon    longitude (radians)
!      lat    latitude (radians)
!      p/z    pressure (Pa) / height (m)
!  zcoords    1 if z is specified, 0 if p is specified
!
!  the functions will return:
!        p    pressure if z is specified (Pa)
!        z    geopotential height if p is specified (m)
!        u    zonal wind (m s^-1)
!        v    meridional wind (m s^-1)
!        t    temperature (K)
!      ptv    virtual potential temperature (K)
!     phis    surface geopotential (m^2 s^-2)
!      rho    density (kj m^-3)
!       qv    specific humidity (kg/kg)
!
!  Initial data are currently identical to:
!
!       Reed, K. A., and C. Jablonowski, 2011: An analytic
!       vortex initialization technique for idealized tropical
!       cyclone studies in AGCMs. Mon. Wea. Rev., 139, 689-710.
!
!  Author: Kevin A. Reed
!          Stony Brook University
!          Email: kevin.a.reed@stonybrook.edu
!
! ======================================================================
! Revision history:
!
!   2023-08: Revised by Li Dong following Yi Zhang's work.
! ======================================================================

  use mpi
  use string
  use flogger
  use math_mod
  use const_mod, only: r8, pi, rad, deg, rd, g, cpd, rd_o_cpd, omega, a => radius
  use block_mod
  use formula_mod
  use vert_coord_mod
  use namelist_mod
  use tracer_mod
  use latlon_parallel_mod

  implicit none

  private

  public tropical_cyclone_test_init
  public tropical_cyclone_test_set_ic

  ! ====================================================================
  ! Physical constants
  ! ====================================================================

  real(r8), parameter ::   &
    Lvap     = 2.5d6     , & ! Latent heat of vaporization of water
    Mvap     = 0.608d0   , & ! Ratio of molar mass of dry air/water
    p0       = 100000.0d0    ! Surface pressure (Pa)

  ! ====================================================================
  ! Test case parameters
  ! ====================================================================
  real(r8), parameter ::   &
    rp       = 282000.d0 , & ! Radius for calculation of PS
    dp       = 1115.d0   , & ! Delta P for calculation of PS
    zp       = 7000.d0   , & ! Height for calculation of P
    q0       = 0.021d0   , & ! qv at surface from Jordan
    gamma    = 0.007d0   , & ! lapse rate
    Ts0      = 302.15d0  , & ! Surface temperature (SST)
    p00      = 101500.d0 , & ! global mean surface pressure
    cen_lat  = 10 * rad  , & ! Center latitude of initial vortex
    cen_lon  = 180 * rad , & ! Center longitufe of initial vortex
    zq1      = 3000.d0   , & ! Height 1 for qv calculation
    zq2      = 8000.d0   , & ! Height 2 for qv calculation
    exppr    = 1.5d0     , & ! Exponent for r dependence of p
    exppz    = 2.d0      , & ! Exponent for z dependence of p
    ztrop    = 15000.d0  , & ! Tropopause Height
    qtrop    = 1.d-11    , & ! Tropopause specific humidity
    rfpi     = 1000000.d0, & ! Radius within which to use fixed-point iteration
    constTv  = 0.608d0   , & ! Constant for Virtual Temp Conversion
    deltaz   = 2.d-13    , & ! Small number to ensure convergence in FPI
    T0       = Ts0 * (1 + constTv * q0), & ! Surface temperature
    Ttrop    = T0 - gamma * ztrop          ! Tropopause temperature
  real(r8) exponent          ! rd * gamma / g
  real(r8) ptrop             ! Tropopause pressure

  real(r8), pointer :: gaussx(:) => null()
  real(r8), pointer :: gaussw(:) => null()

  integer :: ngauss = 5
  integer :: max_iter = 200

contains

  subroutine tropical_cyclone_test_init()

    select case (ngauss)
    case (3)
      gaussx => gaussx3
      gaussw => gaussw3
    case (5)
      gaussx => gaussx5
      gaussw => gaussw5
    case (20)
      gaussx => gaussx20
      gaussw => gaussw20
    case default
      if (proc%is_root()) call log_error('supercell_test_init: Invalid ngauss value!', __FILE__, __LINE__)
    end select

    if (idx_qv == 0) call tracer_add('moist', dt_adv, 'qv', 'Water vapor', 'kg kg-1')

  end subroutine tropical_cyclone_test_init

  real(r8) function get_dry_air_pressure(lon, lat, ptop, ztop, z, ngauss, gaussx, gaussw) result(res)

    real(r8), intent(in) :: lon
    real(r8), intent(in) :: lat
    real(r8), intent(in) :: ptop
    real(r8), intent(in) :: ztop
    real(r8), intent(in) :: z
    integer , intent(in) :: ngauss
    real(r8), intent(in) :: gaussx(ngauss)
    real(r8), intent(in) :: gaussw(ngauss)

    integer jgw
    real(r8) z1, z2, a, b, zg
    real(r8) p, ptv, rho, qv

    ! Set vertical height range.
    z1 = z
    z2 = ztop
    ! Transform height into Gaussian quadrature coordinate [-1,1].
    a = 0.5_r8 * (z2 - z1)
    b = 0.5_r8 * (z1 + z2)
    ! Integrate hydrostatic equation to get dry air surface pressure.
    res = 0
    do jgw = 1, ngauss
      zg = a * gaussx(jgw) + b
      call tropical_cyclone_test(lon, lat, p, zg, 1, ptv, rho, qv)
      ! Remove water vapor from integration.
      ! Note: qv is wet mixing ratio or specific humidity here.
      res = res + gaussw(jgw) * rho * (1 - qv)
    end do
    res = a * g * res + ptop

  end function get_dry_air_pressure

  real(r8) function get_pressure(lon, lat, z) result(res)

    real(r8), intent(in   ) :: lon
    real(r8), intent(in   ) :: lat
    real(r8), intent(inout) :: z

    real(r8) ptv, rho, qv

    call tropical_cyclone_test(lon, lat, res, z, 1, ptv, rho, qv)

  end function get_pressure

  real(r8) function get_top_height(lon, lat, ptop) result(res)

    real(r8), intent(in   ) :: lon
    real(r8), intent(in   ) :: lat
    real(r8), intent(inout) :: ptop

    real(r8) ptv, rho, qv

    call tropical_cyclone_test(lon, lat, ptop, res, 0, ptv, rho, qv)

  end function get_top_height

  real(r8) function get_height(lon, lat, pabv, zabv, ps, p) result(res)

    real(r8), intent(in) :: lon
    real(r8), intent(in) :: lat
    real(r8), intent(in) :: pabv ! Dry air pressure of above level (Pa)
    real(r8), intent(in) :: zabv ! Height of above level (m)
    real(r8), intent(in) :: ps   ! Dry air surface pressure (Pa)
    real(r8), intent(in) :: p    ! Dry air pressure (Pa)

    real(r8), parameter :: eps = 1.0e-12_r8
    real(r8) z1, z2, zc
    real(r8) p1, p2, pc
    integer i

    z1 = 0   ; p1 = ps
    z2 = zabv; p2 = pabv
    do i = 1, 1000
      zc = z2 + (z2 - z1) / (p2 - p1) * (p - p2)
      pc = get_dry_air_pressure(lon, lat, pabv, zabv, zc, ngauss, gaussx, gaussw)
      if (abs(pc - p) / p < eps) exit
      if (pc > p) then
        z2 = zc; p2 = pc
      else
        z1 = zc; p1 = pc
      end if
    end do
    if (i == 1001) then
      call log_error('get_height: Iteration did not converge!', __FILE__, __LINE__)
    end if
    res = zc

  end function get_height

  subroutine tropical_cyclone_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    integer i, j, k
    real(r8) ptv, rho, qv
    real(8) time1, time2

    time1 = MPI_Wtime()

    exponent = rd * gamma / g
    ptrop = p00 * (Ttrop / T0)**(1.d0 / exponent)

    associate (mesh   => block%mesh,             &
               lon    => block%mesh%full_lon   , &
               lat    => block%mesh%full_lat   , &
               gzs    => block%static%gzs      , &
               mgs    => block%dstate(1)%mgs   , &
               mg_lev => block%dstate(1)%mg_lev, &
               p      => block%dstate(1)%p     , &
               p_lev  => block%dstate(1)%p_lev , &
               z      => block%dstate(1)%gz    , &
               z_lev  => block%dstate(1)%gz_lev, &
               u      => block%aux      %u     , &
               u_lon  => block%dstate(1)%u_lon , &
               v      => block%aux      %v     , &
               v_lat  => block%dstate(1)%v_lat , &
               t      => block%aux%t           , &
               pt     => block%dstate(1)%pt    , &
               q      => tracers(block%id)%q   )
    if (proc%is_root()) call log_notice('Calculate model top height and surface dry air pressure.')
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        z_lev%d(i,j,1) = get_top_height(lon(i), lat(j), ptop)
        mgs%d(i,j) = get_dry_air_pressure(lon(i), lat(j), ptop, z_lev%d(i,j,1), 0.0_r8, 20, gaussx20, gaussw20)
      end do
    end do
    call fill_halo(mgs)
    if (proc%is_root()) call log_notice('Calculate dry air pressure at half levels.')
    do k = mesh%half_kds, mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          mg_lev%d(i,j,k) = vert_coord_calc_mg_lev(k, mgs%d(i,j))
        end do
      end do
    end do
    if (proc%is_root()) call log_notice('Calculate height at half levels.')
    do k = mesh%half_kds + 1, mesh%half_kde ! Skip the first half level, since ztop has already been calculated.
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          z_lev%d(i,j,k) = get_height(lon(i), lat(j), mg_lev%d(i,j,k-1), z_lev%d(i,j,k-1), mgs%d(i,j), mg_lev%d(i,j,k))
        end do
      end do
    end do
    where (z_lev%d(:,:,mesh%half_kde) /= 0) z_lev%d(:,:,mesh%half_kde) = 0
    if (proc%is_root()) call log_notice('Calculate height at full levels.')
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          z%d(i,j,k) = 0.5_r8 * (z_lev%d(i,j,k) + z_lev%d(i,j,k+1))
        end do
      end do
    end do
    if (proc%is_root()) call log_notice('Calculate variables on model full levels.')
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          call tropical_cyclone_test(lon(i), lat(j), p%d(i,j,k), z%d(i,j,k), 1, &
            ptv, rho, qv, u%d(i,j,k), v%d(i,j,k), t%d(i,j,k))
          q%d(i,j,k,idx_qv) = dry_mixing_ratio(qv, qv)
          pt%d(i,j,k) = modified_potential_temperature(t%d(i,j,k), p%d(i,j,k), q%d(i,j,k,idx_qv))
        end do
      end do
    end do
    call fill_halo(p)
    call fill_halo(u)
    call fill_halo(v)
    call fill_halo(q, idx_qv)
    call fill_halo(pt)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%half_ids, mesh%half_ide
          u_lon%d(i,j,k) = 0.5_r8 * (u%d(i,j,k) + u%d(i+1,j,k))
        end do
      end do
    end do
    call fill_halo(u_lon)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          v_lat%d(i,j,k) = 0.5_r8 * (v%d(i,j,k) + v%d(i,j+1,k))
        end do
      end do
    end do
    call fill_halo(v_lat)
    if (proc%is_root()) call log_notice('Calculate variables on model half levels.')
    do k = mesh%half_kds, mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          call tropical_cyclone_test(lon(i), lat(j), p_lev%d(i,j,k), z_lev%d(i,j,k), 1, ptv, rho, qv)
        end do
      end do
    end do
    call fill_halo(p_lev)
    z_lev%d = z_lev%d * g
    call fill_halo(z_lev)
    z%d = z%d * g
    call fill_halo(z)
    end associate

    time2 = MPI_Wtime()
    if (proc%is_root()) call log_notice('Time to set tropical_cyclone IC: ' // to_str(time2 - time1, 5) // ' seconds')

  end subroutine tropical_cyclone_test_set_ic

  subroutine tropical_cyclone_test(lon, lat, p, z, zcoords, ptv, rho, qv, u, v, t, phis)

    real(r8), intent(in   )           :: lon
    real(r8), intent(in   )           :: lat
    real(r8), intent(inout)           :: p
    real(r8), intent(inout)           :: z
    integer , intent(in   )           :: zcoords ! 1 if z coordinates are specified
                                                 ! 0 if p coordinates are specified
    real(r8), intent(  out)           :: ptv
    real(r8), intent(  out)           :: rho
    real(r8), intent(  out)           :: qv
    real(r8), intent(  out), optional :: u
    real(r8), intent(  out), optional :: v
    real(r8), intent(  out), optional :: t
    real(r8), intent(  out), optional :: phis

    real(r8)  ps, d1, d2, d, vfac, ufac, gr, f, zn, exner
    integer   n

    !------------------------------------------------
    !   Define Great circle distance (gr) and
    !   Coriolis parameter (f)
    !------------------------------------------------
    ! Coriolis parameter
    f  = 2 * omega * sin(cen_lat)
    ! Great circle radius
    gr = a * acos(sin(cen_lat) * sin(lat) + cos(cen_lat) * cos(lat) * cos(lon - cen_lon))

    ps = p00 - dp * exp(-(gr / rp)**exppr)

    if (zcoords == 1) then
      ! Eqn (8):
      if (z > ztrop) then
        p = ptrop * exp(-(g * (z - ztrop))/(rd * Ttrop))
      else
        p = (p00 - dp * exp(-(gr / rp)**exppr) &
            * exp(-(z / zp)**exppz)) &
            * ((T0 - gamma * z) / T0)**(1 / exponent)
      end if
    else
      ! Eqn (24):
      if (ps >= p .and. p >= ptrop) then
        z = (T0 / gamma) * (1 - (p / ps)**exponent)
      else if (ptrop >= p) then
        z = ztrop + rd * Ttrop / g * log(ptrop / p)
      else
        call log_error('p < ps!', __FILE__, __LINE__)
      end if
      ! If inside a certain distance of the center of the storm
      ! perform a Fixed-point iteration to calculate the height
      ! more accurately
      if (gr < rfpi) then
        do n = 1, 21
          zn = z - fpif(p, gr, z) / fpidfdz(gr, z)
          if (n == 21) then
            write(*, *) 'FPI did not converge after 20 interations in qv & T!!!'
          else if (abs(zn - z) / abs(zn) > deltaz) then
            z = zn
          else
            exit
          end if
        end do
        z = zn
      end if
    end if

    if (present(u) .and. present(v)) then
      d1 = sin(cen_lat) * cos(lat) - &
           cos(cen_lat) * sin(lat) * cos(lon - cen_lon)
      d2 = cos(cen_lat) * sin(lon - cen_lon)
      d  = max(1.0e-25_r8, sqrt(d1**2 + d2**2))
      ufac = d1 / d
      vfac = d2 / d

      if (z > ztrop) then
        u = 0
        v = 0
      else
        v = vfac * (-f * gr / 2.d0 + sqrt((f * gr / 2.d0)**2 &
            - exppr * (gr / rp)**exppr * rd * (T0 - gamma * z) &
            / (exppz * z * rd * (T0 - gamma * z) / (g * zp**exppz) &
            + (1 - p00 / dp * exp((gr / rp)**exppr) * exp((z / zp)**exppz)))))
        u = ufac * (-f * gr / 2.d0 + sqrt((f * gr / 2.d0)**2 &
            - exppr * (gr / rp)**exppr * rd * (T0 - gamma * z) &
            / (exppz * z * rd * (T0 - gamma * z) / (g * zp**exppz) &
            + (1 - p00 / dp * exp((gr / rp)**exppr) * exp((z / zp)**exppz)))))
      end if
    end if

    if (z > ztrop) then
      qv = qtrop
    else
      qv = q0 * exp(-z / zq1) * exp(-(z / zq2)**exppz)
    end if

    exner = (p / p0)**rd_o_cpd
    if (z > ztrop) then
      ptv = Ttrop * (1 + constTv * qv) / exner
    else
      ptv = (T0 - gamma * z) / (1 + constTv * qv) &
          / (1 + exppz * rd * (T0 - gamma * z) * z &
          / (g * zp**exppz * (1 - p00 / dp * exp((gr / rp)**exppr) &
          * exp((z / zp)**exppz)))) * (1 + constTv * qv) / exner
    end if

    if (present(t)) then
      t = ptv * exner / (1 + constTv * qv)
    end if

    if (present(phis)) phis = 0

    rho = p / (rd * exner * ptv)

  end subroutine tropical_cyclone_test

!-----------------------------------------------------------------------
!    First function for fixed point iterations
!-----------------------------------------------------------------------
  real(r8) function fpif(p, gr, z) result(res)

    real(r8), intent(in) :: p, gr, z

    if (z >= 0 .and. z <= ztrop) then
      res = p - (p00 - dp * exp(-(gr / rp)**exppr) &
          * exp(-(z / zp)**exppz)) &
          * ((T0 - gamma * z) / T0)**(g / (rd * gamma))
    else
      res = p - (ptrop * exp(-g * (z - ztrop) / (rd * Ttrop)))
    end if

  end function fpif

!-----------------------------------------------------------------------
!    Second function for fixed point iterations
!-----------------------------------------------------------------------
  real(r8) function fpidfdz(gr, z) result(res)

    real(r8), intent(in) :: gr, z

    res = -exppz * dp * z / (zp**2) * exp(-(gr / rp)**exppr) &
        * exp(-(z / zp)**exppz) * ((T0 - gamma * z) / T0)**(g / (rd * gamma)) &
        + g / (rd * T0) * (p00 - dp * exp(-(gr / rp)**exppr) * exp(-(z / zp)**exppz)) &
        * ((T0 - gamma * z) / T0)**(g / (rd * gamma) - 1)

  end function fpidfdz

end module tropical_cyclone_test_mod
