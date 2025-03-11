#ifdef HAS_LAPACK
module supercell_test_mod

!=======================================================================
!
!  Date:  April 22, 2016
!
!  Functions for setting up idealized initial conditions for the
!  Klemp et al. supercell test.  Before sampling the result,
!  supercell_init() must be called.
!
!  SUBROUTINE supercell_test(
!    lon,lat,p,z,zcoords,u,v,t,ptv,ps,rho,q,pert)
!
!  Given a point specified by:
!      lon    longitude (radians)
!      lat    latitude (radians)
!      p/z    pressure (Pa) / height (m)
!  zcoords    1 if z is specified, 0 if p is specified
!     pert    1 if thermal perturbation included, 0 if not
!
!  the functions will return:
!        p    pressure if z is specified (Pa)
!        z    geopotential height if p is specified (m)
!        u    zonal wind (m s^-1)
!        v    meridional wind (m s^-1)
!        t    temperature (K)
!      ptv    virtual potential temperature (K)
!       ps    surface pressure (Pa)
!      rho    density (kj m^-3)
!        q    water vapor mixing ratio (kg/kg)
!
!  Author: Paul Ullrich
!          University of California, Davis
!          Email: paullrich@ucdavis.edu
!
!          Based on a code by Joseph Klemp
!          (National Center for Atmospheric Research)
!
!=======================================================================

  use flogger
  use const_mod
  use namelist_mod
  use math_mod
  use formula_mod
  use vert_coord_mod
  use tracer_mod
  use process_mod
  use latlon_parallel_mod

  implicit none

  private

  public supercell_test_set_params
  public supercell_test_init
  public supercell_test_set_ic

!=======================================================================
!    Test case parameters
!=======================================================================
  integer(4), parameter ::         &
       nz         = 30         ,   & ! number of vertical levels in init
       nphi       = 16               ! number of meridional points in init

  real(8), parameter ::            &
       z1         = 0          ,   & ! lower sample altitude
       z2         = 50000            ! upper sample altitude

  real(8), parameter ::            &
       X          = 120        ,   & ! Earth reduction factor
       theta0     = 300        ,   & ! theta at the equatorial surface
       theta_tr   = 343        ,   & ! theta at the tropopause
       z_tr       = 12000      ,   & ! altitude at the tropopause
       t_tr       = 213        ,   & ! temperature at the tropopause
       pseq       = 100000           ! surface pressure at equator (Pa)

  real(8), parameter ::            &
       us         = 30         ,   & ! maximum zonal wind velocity
       uc         = 15         ,   & ! coordinate reference velocity
       zs         = 5000       ,   & ! lower altitude of maximum velocity
       zt         = 1000             ! transition distance of velocity

  real(8), parameter ::            &
       pert_dtheta = 3         ,   & ! perturbation magnitude
       pert_lonc   = 0         ,   & ! perturbation longitude
       pert_latc   = 0         ,   & ! perturbation latitude
       pert_rh     = 10000     ,   & ! perturbation horiz. halfwidth
       pert_zc     = 1500      ,   & ! perturbation center altitude
       pert_rz     = 1500            ! perturbation vert. halfwidth

!-----------------------------------------------------------------------
!    Coefficients computed from initialization
!-----------------------------------------------------------------------
  integer(4)                  :: initialized = 0

  real(8), dimension(nphi   ) :: phicoord
  real(8), dimension(     nz) :: zcoord
  real(8), dimension(nphi,nz) :: ptvyz
  real(8), dimension(nphi,nz) :: exneryz
  real(8), dimension(     nz) :: qveq

  real(8), pointer :: gaussx(:) => null()
  real(8), pointer :: gaussw(:) => null()

  integer :: ngauss = 5
  integer :: pert = 1
  integer :: max_iter = 200

  namelist /supercell_test_control/ ngauss, pert, max_iter

contains

  subroutine supercell_test_set_params()

    omega = 0
    radius = radius / X

  end subroutine supercell_test_set_params

  subroutine supercell_test_init(file_path)

    character(*), intent(in) :: file_path

    integer ignore

    open(11, file=file_path, status='old', action='read')
    read(11, nml=supercell_test_control, iostat=ignore)
    close(11)

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
    if (idx_qc == 0) call tracer_add('moist', dt_adv, 'qc', 'Cloud water', 'kg kg-1')
    if (idx_qr == 0) call tracer_add('moist', dt_adv, 'qr', 'Rain water' , 'kg kg-1')

    call supercell_init()

  end subroutine supercell_test_init

  real(8) function get_dry_air_pressure(pert, lon, lat, ptop, ztop, z, ngauss, gaussx, gaussw) result(res)

    integer, intent(in) :: pert
    real(8), intent(in) :: lon
    real(8), intent(in) :: lat
    real(8), intent(in) :: ptop
    real(8), intent(in) :: ztop
    real(8), intent(in) :: z
    integer, intent(in) :: ngauss
    real(8), intent(in) :: gaussx(ngauss)
    real(8), intent(in) :: gaussw(ngauss)

    integer jgw
    real(8) z1, z2, a, b, zg
    real(8) p, ptv, rho, qv

    ! Set vertical height range.
    z1 = z
    z2 = ztop
    ! Transform height into Gaussian quadrature coordinate [-1,1].
    a = 0.5d0 * (z2 - z1)
    b = 0.5d0 * (z1 + z2)
    ! Integrate hydrostatic equation to get dry air surface pressure.
    res = 0
    do jgw = 1, ngauss
      zg = a * gaussx(jgw) + b
      call supercell_test(pert, lon, lat, p, zg, 1, ptv, rho, qv)
      ! Remove water vapor from integration.
      ! Note: qv is wet mixing ratio or specific humidity here.
      res = res + gaussw(jgw) * rho * (1 - qv)
    end do
    res = a * g * res + ptop

  end function get_dry_air_pressure

  real(8) function get_pressure(pert, lon, lat, z) result(res)

    integer, intent(in   ) :: pert
    real(8), intent(in   ) :: lon
    real(8), intent(in   ) :: lat
    real(8), intent(inout) :: z

    real(8) u, v, t, ptv, ps, rho, qv

    call supercell_test(pert, lon, lat, res, z, 1, u, v, t, ptv, ps, rho, qv)

  end function get_pressure

  real(8) function get_top_height(pert, lon, lat, ptop) result(res)

    integer, intent(in   ) :: pert
    real(8), intent(in   ) :: lon
    real(8), intent(in   ) :: lat
    real(8), intent(inout) :: ptop

    real(8) ptv, rho, qv

    call supercell_test(pert, lon, lat, ptop, res, 0, ptv, rho, qv)

  end function get_top_height

  real(8) function get_height(pert, lon, lat, pabv, zabv, ps, p) result(res)

    integer, intent(in) :: pert
    real(8), intent(in) :: lon
    real(8), intent(in) :: lat
    real(8), intent(in) :: pabv ! Dry air pressure of above level (Pa)
    real(8), intent(in) :: zabv ! Height of above level (m)
    real(8), intent(in) :: ps   ! Dry air surface pressure (Pa)
    real(8), intent(in) :: p    ! Dry air pressure (Pa)

    real(8), parameter :: eps = 1.0d-12
    real(8) z1, z2, zc
    real(8) p1, p2, pc
    integer i

    z1 = 0   ; p1 = ps
    z2 = zabv; p2 = pabv
    do i = 1, max_iter
      zc = z2 + (z2 - z1) / (p2 - p1) * (p - p2)
      pc = get_dry_air_pressure(pert, lon, lat, pabv, zabv, zc, ngauss, gaussx, gaussw)
      if (abs(pc - p) / p < eps) exit
      if (pc > p) then
        z2 = zc; p2 = pc
      else
        z1 = zc; p1 = pc
      end if
    end do
    if (i == max_iter + 1) then
      call log_error('get_height: Iteration did not converge!', __FILE__, __LINE__)
    end if
    res = zc

  end function get_height

  subroutine supercell_test_set_ic(block)

    use mpi
    use string
    use block_mod

    type(block_type), intent(inout), target :: block

    integer i, j, k
    real(8) ptv, rho, qv
    real(8) time1, time2

    time1 = MPI_Wtime()

    associate (mesh   => block%mesh            , &
               lon    => block%mesh%full_lon   , &
               lat    => block%mesh%full_lat   , &
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
               t      => block%aux      %t     , &
               pt     => block%dstate(1)%pt    , &
               q      => tracers(block%id)%q   )
    if (proc%is_root()) call log_notice('Calculate model top height and surface dry air pressure.')
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        z_lev%d(i,j,1) = get_top_height(0, lon(i), lat(j), ptop)
        mgs%d(i,j) = get_dry_air_pressure(0, lon(i), lat(j), ptop, z_lev%d(i,j,1), 0.0_r8, 20, gaussx20, gaussw20)
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
          z_lev%d(i,j,k) = get_height(0, lon(i), lat(j), mg_lev%d(i,j,k-1), z_lev%d(i,j,k-1), mgs%d(i,j), mg_lev%d(i,j,k))
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
          call supercell_test(pert, lon(i), lat(j), p%d(i,j,k), z%d(i,j,k), 1, &
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
          call supercell_test(pert, lon(i), lat(j), p_lev%d(i,j,k), z_lev%d(i,j,k), 1, ptv, rho, qv)
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
    if (proc%is_root()) call log_notice('Time to set supercell IC: ' // to_str(time2 - time1, 5) // ' seconds')

  end subroutine supercell_test_set_ic

!=======================================================================
!    Generate the supercell initial conditions
!=======================================================================
  subroutine supercell_init()

    ! d/dphi and int(dphi) operators
    real(8), dimension(nphi,nphi) :: ddphi, intphi
    ! Buffer matrices for computing SVD of d/dphi operator
    real(8), dimension(nphi,nphi) :: ddphibak
    real(8), dimension(nphi,nphi) :: svdpu, svdpvt
    real(8), dimension(nphi)      :: svdps
    real(8), dimension(5*nphi)    :: pwork

    ! d/dz and int(dz) operators
    real(8), dimension(nz,nz) :: ddz, intz
    ! Buffer matrices for computing SVD of d/dz operator
    real(8), dimension(nz,nz) :: ddzbak
    real(8), dimension(nz,nz) :: svdzu, svdzvt
    real(8), dimension(nz)    :: svdzs
    real(8), dimension(5*nz)  :: zwork

    ! Buffer data for calculation of SVD
    integer(4) lwork, info

    ! Sampled values of ueq**2 and d/dz(ueq**2)
    real(8), dimension(nphi,nz) :: ueq2, dueq2

    ! Buffer matrices for iteration
    real(8), dimension(nphi,nz) :: phicoordmat, dztheta, rhs, irhs

    ! Buffer for sampled potential temperature at equator
    real(8), dimension(nz) :: thetaeq

    ! Buffer for computed equatorial Exner pressure and relative humidity
    real(8), dimension(nz) :: exnereq, h

    ! Variables for calculation of equatorial profile
    real(8) exnereqs, p, t, qvs, qv

    ! Error metric
    real(8) err

    ! Loop indices
    integer(4) i, k, iter

    ! Chebyshev nodes in the phi direction
    do i = 1, nphi
      phicoord(i) = - cos(dble(i - 1) * pi / dble(nphi - 1))
      phicoord(i) = 0.25d0 * pi * (phicoord(i) + 1.0d0)
    end do

    ! Matrix of phis
    do k = 1, nz
      phicoordmat(:,k) = phicoord
    end do

    ! Chebyshev nodes in the z direction
    do k = 1, nz
      zcoord(k) = - cos(dble(k - 1) * pi / dble(nz - 1))
      zcoord(k) = z1 + 0.5d0 * (z2 - z1) * (zcoord(k) + 1.0d0)
    end do

    ! Compute the d/dphi operator
    do i = 1, nphi
      call diff_lagrangian_polynomial_coeffs( &
        nphi, phicoord, ddphi(:,i), phicoord(i))
    end do

    ! Zero derivative at pole
    ddphi(:,nphi) = 0.0d0

    ! Compute the d/dz operator
    do k = 1, nz
      call diff_lagrangian_polynomial_coeffs( &
        nz, zcoord, ddz(:,k), zcoord(k))
    end do

    ! Compute the int(dphi) operator via pseudoinverse
    lwork = 5 * nphi

    ddphibak = ddphi
    call DGESVD('A', 'A', &
       nphi, nphi, ddphibak, nphi, &
       svdps, svdpu, nphi, svdpvt, nphi, &
       pwork, lwork, info)

    if (info /= 0) then
      call log_error('Unable to compute SVD of d/dphi matrix')
    end if

    do i = 1, nphi
      if (abs(svdps(i)) <= 1.0d-12) then
        call DSCAL(nphi, 0.0d0, svdpu(1,i), 1)
      else
        call DSCAL(nphi, 1.0d0 / svdps(i), svdpu(1,i), 1)
      end if
    end do
    call DGEMM('T', 'T', &
      nphi, nphi, nphi, 1.0d0, svdpvt, nphi, svdpu, nphi, 0.0d0, &
      intphi, nphi)

    ! Compute the int(dz) operator via pseudoinverse
    lwork = 5*nz

    ddzbak = ddz
    call DGESVD('A', 'A', &
       nz, nz, ddzbak, nz, &
       svdzs, svdzu, nz, svdzvt, nz, &
       zwork, lwork, info)

    if (info /= 0) then
      call log_error('Unable to compute SVD of d/dz matrix')
    end if

    do i = 1, nz
      if (abs(svdzs(i)) <= 1.0d-12) then
        call DSCAL(nz, 0.0d0, svdzu(1,i), 1)
      else
        call DSCAL(nz, 1.0d0 / svdzs(i), svdzu(1,i), 1)
      end if
    end do
    call DGEMM('T', 'T', &
      nz, nz, nz, 1.0d0, svdzvt, nz, svdzu, nz, 0.0d0, &
      intz, nz)

    ! Sample the equatorial velocity field and its derivative
    do k = 1, nz
      ueq2(1,k) = zonal_velocity(zcoord(k), 0.0d0)
      ueq2(1,k) = ueq2(1,k)**2
    end do
    do k = 1, nz
      dueq2(1,k) = dot_product(ddz(:,k), ueq2(1,:))
    end do
    do i = 2, nphi
      ueq2 (i,:) = ueq2 (1,:)
      dueq2(i,:) = dueq2(1,:)
    end do

    ! Initialize potential temperature at equator
    do k = 1, nz
      thetaeq(k) = equator_theta(zcoord(k))
      H(k) = equator_relative_humidity(zcoord(k))
    end do
    ptvyz(1,:) = thetaeq

    ! Exner pressure at the equatorial surface
    exnereqs = (pseq / p0)**rd_o_cpd

    ! Iterate on equatorial profile
    do iter = 1, 12
      ! Calculate Exner pressure in equatorial column (p0 at surface)
      rhs(1,:) = - g / cpd / ptvyz(1,:)
      do k = 1, nz
        exnereq(k) = dot_product(intz(:,k), rhs(1,:))
      end do
      do k = 2, nz
        exnereq(k) = exnereq(k) + (exnereqs - exnereq(1))
      end do
      exnereq(1) = exnereqs

      ! Calculate new pressure and temperature
      do k = 1, nz
        p = p0 * exnereq(k)**(1.0d0 / rd_o_cpd)
        T = thetaeq(k) * exnereq(k)

        qvs = saturation_mixing_ratio(p, T)
        qveq(k) = qvs * H(k)

        ptvyz(1,k) = thetaeq(k) * (1.d0 + 0.61d0 * qveq(k))
      end do
    end do

    ! Iterate on remainder of domain
    do iter = 1, 12
      ! Compute d/dz(theta)
      do i = 1, nphi
        do k = 1, nz
          dztheta(i,k) = dot_product(ddz(:,k), ptvyz(i,:))
        end do
      end do

      ! Compute rhs
      rhs = sin(2.0d0*phicoordmat)/(2.0d0*g) &
            * (ueq2 * dztheta - ptvyz * dueq2)

      ! Integrate
      do k = 1, nz
        do i = 1, nphi
          irhs(i,k) = dot_product(intphi(:,i), rhs(:,k))
        end do
      end do

      ! Apply boundary conditions (fixed Dirichlet condition at equator)
      do i = 2, nphi
        irhs(i,:) = irhs(i,:) + (ptvyz(1,:) - irhs(1,:))
      end do
      irhs(1,:) = ptvyz(1,:)

      ! Update iteration
      ptvyz = irhs
    end do

    ! Calculate pressure through remainder of domain
    rhs = - ueq2 * sin(phicoordmat) * cos(phicoordmat) / cpd / ptvyz

    do k = 1, nz
      do i = 1, nphi
        exneryz(i,k) = dot_product(intphi(:,i), rhs(:,k))
      end do
      do i = 2, nphi
        exneryz(i,k) = exneryz(i,k) + (exnereq(k) - exneryz(1,k))
      end do

      exneryz(1,k) = exnereq(k)
    end do

    ! Initialization successful
    initialized = 1

  end subroutine supercell_init

!-----------------------------------------------------------------------
!    Evaluate the supercell initial conditions
!-----------------------------------------------------------------------
  subroutine supercell_test(pert, lon, lat, p, z, zcoords, ptv, rho, qv, u, v, t, ps)

    integer, intent(in   )           :: pert    ! 1 if perturbation should be included
                                                ! 0 if no perturbation should be included
    real(8), intent(in   )           :: lon     ! Longitude (radians)
    real(8), intent(in   )           :: lat     ! Latitude (radians)
    real(8), intent(inout)           :: p       ! Pressure (Pa)
    real(8), intent(inout)           :: z       ! Altitude (m)
    integer, intent(in   )           :: zcoords ! 1 if z coordinates are specified
                                                ! 0 if p coordinates are specified
    real(8), intent(  out)           :: ptv     ! Virtual potential Temperature (K)
    real(8), intent(  out)           :: rho     ! Density (kg m^-3)
    real(8), intent(  out)           :: qv      ! Water vapor mixing ratio (kg/kg)
    real(8), intent(  out), optional :: u       ! Zonal wind (m s^-1)
    real(8), intent(  out), optional :: v       ! Meridional wind (m s^-1)
    real(8), intent(  out), optional :: t       ! Temperature (K)
    real(8), intent(  out), optional :: ps      ! Surface Pressure (Pa)

    ! Check that we are initialized
    if (initialized /= 1) then
      call log_error('supercell_init() has not been called', __FILE__, __LINE__)
    end if

    !------------------------------------------------
    !   Begin sampling
    !------------------------------------------------

    ! Sample surface pressure
    if (present(ps)) then
      call supercell_z(pert, lon, lat, 0.d0, ps, ptv, rho, qv)
    end if

    ! Calculate dependent variables
    if (zcoords == 1) then
      call supercell_z(pert, lon, lat, z, p, ptv, rho, qv)
    else
      call supercell_p(pert, lon, lat, p, z, ptv, rho, qv)
    end if

    ! Sample the zonal velocity
    if (present(u)) then
      u = zonal_velocity(z, lat)
    end if

    ! Zero meridional velocity
    if (present(v)) then
      v = 0
    end if

    ! Temperature
    if (present(t)) then
      t = ptv / (1.d0 + 0.61d0 * qv) * (p / p0)**rd_o_cpd
    end if

  end subroutine supercell_test

!-----------------------------------------------------------------------
!    Calculate pointwise pressure and temperature
!-----------------------------------------------------------------------
  subroutine supercell_z(pert, lon, lat, z, p, ptv, rho, qv)

    integer, intent(in ) :: pert  ! 1 if perturbation should be included
                                  ! 0 if no perturbation should be included
    real(8), intent(in ) :: lon   ! Longitude (radians)
    real(8), intent(in ) :: lat   ! Latitude (radians)
    real(8), intent(in ) :: z     ! Altitude (m)
    real(8), intent(out) :: p
    real(8), intent(out) :: ptv
    real(8), intent(out) :: rho
    real(8), intent(out) :: qv

    ! Northern hemisphere latitude
    real(8) abs_lat

    ! Pointwise Exner pressure
    real(8) exner

    ! Assembled variable values in a column
    real(8), dimension(nz) :: varcol

    ! Coefficients for computing a polynomial fit in each coordinate
    real(8), dimension(nphi) :: fitphi
    real(8), dimension(nz)   :: fitz

    integer k

    ! Northern hemisphere latitude
    abs_lat = abs(lat)

    ! Perform fit
    call lagrangian_polynomial_coeffs(nz, zcoord, fitz, z)
    call lagrangian_polynomial_coeffs(nphi, phicoord, fitphi, abs_lat)

    ! Obtain exner pressure of background state
    do k = 1, nz
      varcol(k) = dot_product(fitphi, exneryz(:,k))
    end do
    exner = dot_product(fitz, varcol)
    p = p0 * exner**(1.0d0 / rd_o_cpd)

    ! Sample the initialized fit at this point for theta_v
    do k = 1, nz
      varcol(k) = dot_product(fitphi, ptvyz(:,k))
    end do
    ptv = dot_product(fitz, varcol)

    ! Sample water vapor mixing ratio
    qv = dot_product(fitz, qveq)

    ! Fixed density
    rho = p / (rd * exner * ptv)

    ! Modified virtual potential temperature
    if (pert /= 0) then
      ptv = ptv + thermal_perturbation(lon, lat, z) * (1 + 0.61d0 * qv)
    end if

    ! Updated pressure
    p = p0 * (rho * rd * ptv / p0)**(cpd / (cpd - rd))

  end subroutine supercell_z

!-----------------------------------------------------------------------
!    Calculate pointwise z and temperature given pressure
!-----------------------------------------------------------------------
  subroutine supercell_p(pert, lon, lat, p, z, ptv, rho, qv)

    integer, intent(in ) :: pert   ! 1 if perturbation should be included
                                   ! 0 if no perturbation should be included
    real(8), intent(in ) :: lon    ! Longitude (radians)
    real(8), intent(in ) :: lat    ! Latitude (radians)
    real(8), intent(in ) :: p      ! Pressure (Pa)
    real(8), intent(out) :: z
    real(8), intent(out) :: ptv
    real(8), intent(out) :: rho
    real(8), intent(out) :: qv

    ! Bounding interval and sampled values
    real(8) za, zb, zc, pa, pb, pc

    integer iter

    za = z1
    zb = z2

    call supercell_z(pert, lon, lat, za, pa, ptv, rho, qv)
    call supercell_z(pert, lon, lat, zb, pb, ptv, rho, qv)

    if (pa < p) then
      write(*,*) 'Requested pressure out of range on bottom, adjust sample interval'
      write(*,*) pa, p
      stop
    end if
    if (pb > p) then
      write(*,*) 'Requested pressure out of range on top, adjust sample interval'
      write(*,*) pb, p
      stop
    end if

    ! Iterate using fixed point method
    do iter = 1, 1000
      zc = (za * (pb - p) - zb * (pa - p)) / (pb - pa)

      call supercell_z(pert, lon, lat, zc, pc, ptv, rho, qv)

      if (abs((pc - p) / p) < 1.d-12) exit

      if (pc > p) then
        za = zc
        pa = pc
      else
        zb = zc
        pb = pc
      end if
    end do

    if (iter == 1001) then
      call log_error('Iteration failed to converge', __FILE__, __LINE__)
    end if

    z = zc

  end subroutine supercell_p

!-----------------------------------------------------------------------
!    Calculate pointwise z and temperature given pressure
!-----------------------------------------------------------------------
  pure real(8) function thermal_perturbation(lon, lat, z) result(res)

    real(8), intent(in) :: lon ! Longitude (radians)
    real(8), intent(in) :: lat ! Latitude (radians)
    real(8), intent(in) :: z   ! Altitude (m)

    ! Great circle radius from the perturbation centerpoint
    real(8) gr

    ! Approximately spherical radius from the perturbation centerpoint
    real(8) Rtheta

    gr = radius * acos(sin(pert_latc * rad) * sin(lat) + &
         (cos(pert_latc * rad) * cos(lat) * cos(lon - pert_lonc * rad)))

    Rtheta = sqrt((gr / pert_rh)**2 + ((z - pert_zc) / pert_rz)**2)

    if (Rtheta <= 1) then
      res = pert_dtheta * cos(0.5d0 * pi * Rtheta)**2
    else
      res = 0
    end if

  end function thermal_perturbation

!-----------------------------------------------------------------------
!    Calculate the reference zonal velocity
!-----------------------------------------------------------------------
  pure real(8) function zonal_velocity(z, lat) result(res)

    real(8), intent(in) :: z, lat

    if (z <= zs - zt) then
      res = us * (z / zs) - uc
    elseif (abs(z - zs) <= zt) then
      res = (-4.0d0/5.0d0 + 3.0d0*z/zs - 5.0d0/4.0d0*(z**2)/(zs**2)) * us - uc
    else
      res = us - uc
    end if

    res = res * cos(lat)

  end function zonal_velocity

!-----------------------------------------------------------------------
!    Calculate pointwise theta at the equator at the given altitude
!-----------------------------------------------------------------------
  pure real(8) function equator_theta(z) result(res)

    real(8), intent(in) :: z

    if (z <= z_tr) then
      res = theta0 + (theta_tr - theta0) * (z / z_tr)**1.25d0
    else
      res = theta_tr * exp(g / cpd / t_tr * (z - z_tr))
    end if

  end function equator_theta

!-----------------------------------------------------------------------
!    Calculate pointwise relative humidity (in %) at the equator at the
!    given altitude
!-----------------------------------------------------------------------
  pure real(8) function equator_relative_humidity(z) result(res)

    real(8), intent(in) :: z

    if (z <= z_tr) then
      res = 1.0d0 - 0.75d0 * (z / z_tr)**1.25d0
    else
      res = 0.25d0
    end if

  end function equator_relative_humidity

!-----------------------------------------------------------------------
!    Calculate saturation mixing ratio (in kg/kg) in terms of pressure
!    (in Pa) and temperature (in K)
!-----------------------------------------------------------------------
  pure real(8) function saturation_mixing_ratio(p, t) result(res)

    real(8), intent(in) :: p ! Pressure in Pa
    real(8), intent(in) :: t ! Temperature in K

    res = 380.0d0 / p * exp(17.27d0 * (t - 273.d0) / (t - 36.0d0))

    if (res > 0.014) res = 0.014

  end function saturation_mixing_ratio

!-----------------------------------------------------------------------
!    Calculate coefficients for a Lagrangian polynomial
!-----------------------------------------------------------------------
  subroutine lagrangian_polynomial_coeffs(npts, x, coeffs, xs)

    ! Number of points to fit
    integer(4), intent(in) :: npts

    ! Sample points to fit
    real(8), dimension(npts), intent(in) :: x

    ! Computed coefficients
    real(8), dimension(npts), intent(out) :: coeffs

    ! Point at which sample is taken
    real(8), intent(in) :: xs

    integer(4) i, j

    ! Compute the lagrangian polynomial coefficients
    do i = 1, npts
      coeffs(i) = 1.0d0
      do j = 1, npts
        if (i == j) cycle
        coeffs(i) = coeffs(i) * (xs - x(j)) / (x(i) - x(j))
      end do
    end do

  end subroutine lagrangian_polynomial_coeffs

!-----------------------------------------------------------------------
!    Calculate coefficients of the derivative of a Lagrangian polynomial
!-----------------------------------------------------------------------
  subroutine diff_lagrangian_polynomial_coeffs(npts, x, coeffs, xs)

    ! Number of points to fit
    integer(4), intent(in) :: npts

    ! Sample points to fit
    real(8), dimension(npts), intent(in) :: x

    ! Computed coefficients
    real(8), dimension(npts), intent(out) :: coeffs

    ! Point at which sample is taken
    real(8), intent(in) :: xs

    integer(4) i, j, imatch

    ! Buffer sum
    real(8) coeffsum, differential

    ! Check if xs is equivalent to one of the values of x
    imatch = -1
    do i = 1, npts
      if (abs(xs - x(i)) < 1.0d-14) then
        imatch = i
        exit
      end if
    end do

    ! Equivalence detected; special treatment required
    if (imatch /= -1) then
      do i = 1, npts
        coeffs(i) = 1.0d0
        coeffsum = 0.0d0

        do j = 1, npts
          if (j == i .or. j == imatch) then
            cycle
          end if

          coeffs(i) = coeffs(i) * (xs - x(j)) / (x(i) - x(j))
          coeffsum = coeffsum + 1.0 / (xs - x(j))
        end do

        if (i /= imatch) then
          coeffs(i) = coeffs(i)                   &
            * (1.0 + (xs - x(imatch)) * coeffsum) &
            / (x(i) - x(imatch))
        else
          coeffs(i) = coeffs(i) * coeffsum
        end if
      end do

    ! No equivalence; simply differentiate lagrangian fit
    else
      call lagrangian_polynomial_coeffs(npts, x, coeffs, xs)

      do i = 1, npts
        differential = 0.0d0
        do j = 1, npts
          if (i == j) cycle
          differential = differential + 1.0 / (xs - x(j))
        end do
        coeffs(i) = coeffs(i) * differential
      end do
    end if

  end subroutine diff_lagrangian_polynomial_coeffs

end module supercell_test_mod
#endif