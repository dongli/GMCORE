module rossby_haurwitz_wave_test_mod

  use flogger
  use const_mod
  use namelist_mod
  use block_mod
  use latlon_parallel_mod

  implicit none

  private

  public rossby_haurwitz_wave_test_set_ic

  real(r8), parameter :: R = 4.0d0
  real(r8), parameter :: omg = 7.848d-6
  real(r8)            :: gz0

contains

  ! u = a ω (cosφ + R cosᴿ⁻¹φ sin²φ cosRλ - cosᴿ⁺¹φ cosRλ)
  !
  ! v = - a ω R cosᴿ⁻¹φ sinφ sinRλ
  !
  ! gz = gz0 + a² A(φ) + a² B(φ) cosRλ + a² C(φ) cos2Rλ
  !
  ! A(φ) = 1/2 ω (2 Ω + ω) cos²φ + 1/4 ω² cos²ᴿφ ((R + 1) cos²φ + (2 R² - R - 2) - 2 R² cos⁻²φ)
  ! B(φ) = 2 (Ω + ω) ω cosᴿφ ((R² + 2 R + 2) - (R + 1)² cos²φ) / (R + 1) / (R + 2)
  ! C(φ) = 1/4 ω² cos²ᴿφ ((R + 1) cos²φ - (R + 2))


  subroutine rossby_haurwitz_wave_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    real(r8) lon, cos_lat, sin_lat
    real(r8) a, b, c
    integer i, j

    gz0 = 8.0d3 * g

    associate (mesh   => block%mesh           , &
               u      => block%dstate(1)%u_lon, &
               v      => block%dstate(1)%v_lat, &
               gz     => block%dstate(1)%gz   , &
               gzs    => block%static%gzs)
    gzs%d = 0.0

    do j = mesh%full_jds, mesh%full_jde
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = mesh%half_ids, mesh%half_ide
        lon = mesh%half_lon(i)
        a = cos_lat
        b = R * cos_lat**(R - 1) * sin_lat**2 * cos(R * lon)
        c = cos_lat**(R + 1) * cos(R * lon)
        u%d(i,j,1) = radius * omg * (a + b - c)
      end do
    end do
    call fill_halo(u)

    do j = mesh%half_jds, mesh%half_jde
      cos_lat = mesh%half_cos_lat(j)
      sin_lat = mesh%half_sin_lat(j)
      do i = mesh%full_ids, mesh%full_ide
        lon = mesh%full_lon(i)
        a = R * cos_lat**(R - 1) * sin_lat * sin(R * lon)
        v%d(i,j,1) = - radius * omg * a
      end do
    end do
    call fill_halo(v)

    do j = mesh%full_jds, mesh%full_jde
      cos_lat = mesh%full_cos_lat(j)
      a = 0.5 * omg * (2 * omega + omg) * cos_lat**2 + &
        0.25 * omg**2 * ((R + 1) * cos_lat**(2 * R + 2) + (2 * R**2 - R - 2) * cos_lat**(2 * R) - 2 * R**2 * cos_lat**(2 * R - 2))
      b = 2 * (omega + omg) * omg * cos_lat**R * &
        (R**2 + 2 * R + 2 - (R + 1)**2 * cos_lat**2) / (R + 1) / (R + 2)
      c = 0.25 * omg**2 * cos_lat**(2 * R) * ((R + 1) * cos_lat**2 - R - 2)
      do i = mesh%full_ids, mesh%full_ide
        lon = mesh%full_lon(i)
        gz%d(i,j,1) = gz0 + radius**2 * (a + b * cos(R * lon) + c * cos(2 * R * lon))
      end do
    end do
    call fill_halo(gz)

    init_hydrostatic_gz = .true.
    end associate

  end subroutine rossby_haurwitz_wave_test_set_ic

end module rossby_haurwitz_wave_test_mod
