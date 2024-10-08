module baroclinic_wave_test_mod

  use const_mod
  use formula_mod
  use latlon_parallel_mod
  use block_mod
  use operators_mod

  implicit none

  private

  public baroclinic_wave_test_set_ic

  real(r8), parameter :: alpha = 0.0_r8
  real(r8), parameter :: u0    = 35.0_r8   ! m s-1
  real(r8), parameter :: t0    = 288.0_r8  ! K
  real(r8), parameter :: gamma = 0.005_r8  ! K m-1
  real(r8), parameter :: dt    = 4.8e5_r8  ! K
  real(r8), parameter :: eta0  = 0.252
  real(r8), parameter :: etat  = 0.2       ! Tropopause level
  real(r8), parameter :: lonc  = pi / 9.0
  real(r8), parameter :: latc  = pi2 / 9.0
  real(r8), parameter :: up    = 1.0       ! m s-1

contains

  subroutine baroclinic_wave_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    real(r8) etav, eta, tbar, gzbar, sin_lat, cos_lat, half_lon, r
    integer i, j, k

    associate (mesh   => block%mesh            , &
               u      => block%dstate(1)%u_lon , &
               v      => block%dstate(1)%v_lat , &
               mgs    => block%dstate(1)%mgs   , &
               mg     => block%dstate(1)%mg    , &
               t      => block%aux%t           , &
               pt     => block%dstate(1)%pt    , &
               gz_lev => block%dstate(1)%gz_lev, &
               gz     => block%dstate(1)%gz    , &
               gzs    => block%static%gzs)
    mgs%d = 1.0e5_r8
    v  %d = 0

    call calc_mg(block, block%dstate(1))

    do k = mesh%full_kds, mesh%full_kde
      eta = mesh%full_lev(k)
      etav = (eta - eta0) * pi05
      do j = mesh%full_jds, mesh%full_jde
        sin_lat = mesh%full_sin_lat(j)
        cos_lat = mesh%full_cos_lat(j)
        do i = mesh%half_ids, mesh%half_ide
          half_lon = mesh%half_lon(i)
          r = 10.0 * acos(sin(latc) * sin_lat + cos(latc) * cos_lat * cos(half_lon - lonc))
          u%d(i,j,k) = u0 * cos(etav)**(1.5d0) * sin(2 * mesh%full_lat(j))**2 + &
                       up * merge(exp(-r**2), 0.0_r8, r**2 <= 50) ! FIXME: Why we set a limit on r?
        end do
      end do
    end do
    call fill_halo(u)

    do k = mesh%full_kds, mesh%full_kde
      eta = mesh%full_lev(k)
      etav = (eta - eta0) * pi / 2d0
      if (etat <= eta .and. eta <= 1) then
        tbar = t0 * eta**(Rd * gamma / g)
      else
        tbar = t0 * eta**(Rd * gamma / g) + dt * (etat - eta)**5
      end if
      do j = mesh%full_jds, mesh%full_jde
        sin_lat = mesh%full_sin_lat(j)
        cos_lat = mesh%full_cos_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          t%d(i,j,k) = tbar + 3.0d0 / 4.0d0 * eta * pi * u0 / rd * sin(etav) * sqrt(cos(etav)) * (             &
              (-2 * sin_lat**6 * (cos_lat**2 + 1.0d0 / 3.0d0) + 10.0d0 / 63.0d0) * 2 * u0 * cos(etav)**1.5d0 + &
              (8.0d0 / 5.0d0 * cos_lat**3 * (sin_lat**2 + 2.0d0 / 3.0d0) - pi / 4.0d0) * radius * omega        &
            )
          pt%d(i,j,k) = modified_potential_temperature(t%d(i,j,k), mg%d(i,j,k), 0.0_r8)
        end do
      end do
    end do
    call fill_halo(pt)

    do k = mesh%half_kds, mesh%half_kde
      eta = merge(1.0d-12, mesh%half_lev(k), mesh%half_lev(k) == 0)
      etav = (eta - eta0) * pi / 2d0
      if (etat <= eta .and. eta <= 1) then
        gzbar = t0 * g / gamma * (1 - eta**(Rd * gamma / g))
      else
        gzbar = t0 * g / gamma * (1 - eta**(Rd * gamma / g)) - Rd * dt * (       &
            (log(eta / etat) + 137.0d0 / 60.0d0) * etat**5 - 5 * etat**4 * eta + &
            5 * etat**3 * eta**2 - 10.0d0 / 3.0d0 * etat**2 * eta**3 +           &
            5.0d0 / 4.0d0 * etat * eta**4 - 1.0d0 / 5.0d0 * eta**5               &
          )
      end if
      do j = mesh%full_jds, mesh%full_jde
        sin_lat = mesh%full_sin_lat(j)
        cos_lat = mesh%full_cos_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          gz_lev%d(i,j,k) = gzbar + u0 * cos(etav)**1.5d0 * (                                            &
            (-2 * sin_lat**6 * (cos_lat**2 + 1.0d0 / 3.0d0) + 10.0d0 / 63.0d0) * u0 * cos(etav)**1.5d0 + &
            (8.0d0 / 5.0d0 * cos_lat**3 * (sin_lat**2 + 2.0d0 / 3.0d0) - pi / 4.0d0) * radius * omega    &
          )
        end do
      end do
    end do
    call fill_halo(gz_lev)

    do k = mesh%full_kds, mesh%full_kde
      eta = mesh%full_lev(k)
      etav = (eta - eta0) * pi / 2d0
      if (etat <= eta .and. eta <= 1) then
        gzbar = t0 * g / gamma * (1 - eta**(Rd * gamma / g))
      else
        gzbar = t0 * g / gamma * (1 - eta**(Rd * gamma / g)) - Rd * dt * (   &
            (log(eta / etat) + 137d0 / 60d0) * etat**5 - 5 * etat**4 * eta + &
            5 * etat**3 * eta**2 - 10d0 / 3d0 * etat**2 * eta**3 +           &
            5d0 / 4d0 * etat * eta**4 - 1d0 / 5d0 * eta**5                   &
          )
      end if
      do j = mesh%full_jds, mesh%full_jde
        sin_lat = mesh%full_sin_lat(j)
        cos_lat = mesh%full_cos_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          gz%d(i,j,k) = gzbar + u0 * cos(etav)**1.5d0 * (                                        &
            (-2 * sin_lat**6 * (cos_lat**2 + 1d0 / 3d0) + 10d0 / 63d0) * u0 * cos(etav)**1.5d0 + &
            (8d0 / 5d0 * cos_lat**3 * (sin_lat**2 + 2d0 / 3d0) - pi / 4d0) * radius * omega      &
          )
        end do
      end do
    end do
    call fill_halo(gz)

    etav = (1 - eta0) * pi / 2
    do j = mesh%full_jds, mesh%full_jde
      sin_lat = mesh%full_sin_lat(j)
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_ids, mesh%full_ide
        gzs%d(i,j) = u0 * cos(etav)**1.5d0 * (                                                 &
          (-2 * sin_lat**6 * (cos_lat**2 + 1d0 / 3d0) + 10d0 / 63d0) * u0 * cos(etav)**1.5d0 + &
          (8d0 / 5d0 * cos_lat**3 * (sin_lat**2 + 2d0 / 3d0) - pi / 4d0) * radius * omega      &
        )
      end do
    end do
    call fill_halo(gzs)
    end associate

  end subroutine baroclinic_wave_test_set_ic

end module baroclinic_wave_test_mod
