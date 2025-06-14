module dcmip31_test_mod

  use flogger
  use namelist_mod
  use const_mod, only: r8, pi, rd, rd_o_cpd, cpd, g, p0, radius, omega
  use latlon_parallel_mod
  use block_mod
  use formula_mod
  use operators_mod

  implicit none

  private

  public dcmip31_test_set_params
  public dcmip31_test_set_ic

  real(r8), parameter :: X    = 125.0_r8
  real(r8), parameter :: u0   = 20.0_r8      ! m s-1
  real(r8), parameter :: teq  = 300.0_r8     ! K
  real(r8), parameter :: peq  = 1.0e5_r8     ! Pa
  real(r8), parameter :: d    = 5000.0_r8    ! m
  real(r8), parameter :: lonc = 2 * pi / 3
  real(r8), parameter :: latc = 0.0_r8
  real(r8), parameter :: dpt  = 1.0_r8       ! K
  real(r8), parameter :: lz   = 20000.0_r8   ! m
  real(r8), parameter :: N    = 0.01_r8      ! s-1
  real(r8), parameter :: N2   = N**2
  real(r8)            :: t0

contains

  subroutine dcmip31_test_set_params()

    omega = 0.0_r8
    radius = radius / X
    t0 = g**2 / N2 / cpd

  end subroutine dcmip31_test_set_params

  subroutine dcmip31_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    integer i, j, k
    real(r8) ts, cos_2lat, r

    associate (mesh   => block%mesh            , &
               u      => block%dstate(1)%u_lon , &
               v      => block%dstate(1)%v_lat , &
               w      => block%dstate(1)%w_lev , &
               mgs    => block%dstate(1)%mgs   , &
               mg_lev => block%dstate(1)%mg_lev, &
               mg     => block%dstate(1)%mg    , &
               pt     => block%dstate(1)%pt    , &
               gz_lev => block%dstate(1)%gz_lev, &
               gz     => block%dstate(1)%gz    , &
               gzs    => block%static%gzs)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%half_ids, mesh%half_ide
          u%d(i,j,k) = u0 * mesh%full_cos_lat(j)
        end do
      end do
    end do
    call fill_halo(u)

    v  %d = 0
    gzs%d = 0

    if (nonhydrostatic) w%d = 0

    do j = mesh%full_jds, mesh%full_jde
      cos_2lat = cos(2 * mesh%full_lat(j))
      ts = t0 + (teq - t0) * exp(-u0 * N2 / (4 * g**2) * (u0 + 2 * omega * radius) * (cos_2lat - 1))
      do i = mesh%full_ids, mesh%full_ide
        mgs%d(i,j) = peq * exp(u0 / (4 * t0 * rd) * (u0 + 2 * omega * radius) * (cos_2lat - 1)) * &
                      (ts / teq)**(1.0_r8 / rd_o_cpd)
      end do
    end do
    call fill_halo(mgs)

    call calc_mg(block, block%dstate(1))

    do k = mesh%half_kds, mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde
        cos_2lat = cos(2 * mesh%full_lat(j))
        ts = t0 + (teq - t0) * exp(-u0 * N2 / (4 * g**2) * (u0 + 2 * omega * radius) * (cos_2lat - 1))
        do i = mesh%full_ids, mesh%full_ide
          gz_lev%d(i,j,k) = -g**2 / n2 * log(ts / t0 * ((mg_lev%d(i,j,k) / mgs%d(i,j))**rd_o_cpd - 1) + 1)
        end do
      end do
    end do
    call fill_halo(gz_lev)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          gz%d(i,j,k) = 0.5_r8 * (gz_lev%d(i,j,k) + gz_lev%d(i,j,k+1))
        end do
      end do
    end do
    call fill_halo(gz)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        cos_2lat = cos(2 * mesh%full_lat(j))
        ts = t0 + (teq - t0) * exp(-u0 * N2 / (4 * g**2) * (u0 + 2 * omega * radius) * (cos_2lat - 1))
        do i = mesh%full_ids, mesh%full_ide
          pt%d(i,j,k) = ts * (p0 / mgs%d(i,j))**rd_o_cpd / (ts / t0 * ((mg%d(i,j,k) / mgs%d(i,j))**rd_o_cpd - 1) + 1)
          ! Perturbation
          r = radius * acos(sin(latc) * mesh%full_sin_lat(j) + cos(latc) * mesh%full_cos_lat(j) * cos(mesh%full_lon(i) - lonc))
          pt%d(i,j,k) = pt%d(i,j,k) + dpt * d**2 / (d**2 + r**2) * sin(2 * pi * gz%d(i,j,k) / g / lz)
        end do
      end do
    end do
    call fill_halo(pt)

    init_hydrostatic_gz = .true.
    end associate

  end subroutine dcmip31_test_set_ic

end module dcmip31_test_mod
