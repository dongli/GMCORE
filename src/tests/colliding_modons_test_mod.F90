module colliding_modons_test_mod

  use const_mod, only: r8, pi, radius, omega, g
  use namelist_mod
  use sphere_geometry_mod
  use formula_mod
  use block_mod
  use operators_mod
  use latlon_parallel_mod

  implicit none

  private

  public colliding_modons_test_set_params
  public colliding_modons_test_set_ic

  real(r8), parameter :: z0   = 10000.0_r8  ! m
  real(r8), parameter :: T0   = 300.0_r8    ! K
  real(r8), parameter :: p0   = 1.0e5_r8    ! Pa
  real(r8), parameter :: U0   = 40.0_r8     ! m/s
  real(r8), parameter :: r0   = 500.0e3_r8  ! m
  real(r8), parameter :: lon1 = pi * 0.5_r8 ! rad
  real(r8), parameter :: lat1 = 0           ! rad
  real(r8), parameter :: lon2 = pi * 1.5_r8 ! rad
  real(r8), parameter :: lat2 = 0           ! rad

contains

  subroutine colliding_modons_test_set_params()

    omega = 0

  end subroutine colliding_modons_test_set_params

  subroutine colliding_modons_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    integer i, j, k
    real(r8) r, m1, m2

    associate (mesh   => block%mesh           , &
               dstate => block%dstate(1)      , &
               gz     => block%dstate(1)%gz   , &
               mgs    => block%dstate(1)%mgs  , &
               mg     => block%dstate(1)%mg   , &
               u_lon  => block%dstate(1)%u_lon, &
               v_lat  => block%dstate(1)%v_lat, &
               t      => block%aux      %t    , &
               pt     => block%dstate(1)%pt   )
    if (baroclinic) then
      mgs%d = p0
      call calc_mg (block, dstate)
      call calc_dmg(block, dstate)

      t%d = T0
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            pt%d(i,j,k) = modified_potential_temperature(t%d(i,j,k), mg%d(i,j,k), 0.0_r8)
          end do
        end do
      end do
      call fill_halo(pt)
    else
      gz%d = g * z0
    end if

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          r = great_circle(radius, lon1, lat1, mesh%half_lon(i), mesh%full_lat(j))
          m1 = U0 * exp(-(r / r0)**2)
          r = great_circle(radius, lon2, lat2, mesh%half_lon(i), mesh%full_lat(j))
          m2 = U0 * exp(-(r / r0)**2)
          u_lon%d(i,j,k) = m1 - m2
        end do
      end do
    end do
    call fill_halo(u_lon)

    v_lat%d = 0
    end associate

  end subroutine colliding_modons_test_set_ic

end module colliding_modons_test_mod