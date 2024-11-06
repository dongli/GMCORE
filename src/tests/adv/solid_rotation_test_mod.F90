module solid_rotation_test_mod

  use const_mod
  use namelist_mod, dt => dt_adv
  use sphere_geometry_mod
  use latlon_parallel_mod
  use block_mod
  use adv_mod
  use tracer_mod

  implicit none

  public solid_rotation_test_init
  public solid_rotation_test_set_ic
  public solid_rotation_test_set_uv

  real(8), parameter :: period = 12 * 86400
  real(8), parameter :: h0     = 1000 ! m
  real(8), parameter :: lon0   = 3 * pi / 2.0_r8
  real(8), parameter :: lat0   = 0
  real(8), parameter :: alpha  = 90.0_r8 * rad
  real(8) u0

contains

  subroutine solid_rotation_test_init()

    u0 = pi2 * radius / period

    call tracer_add('adv', dt, 'q1', 'Cosine bell tracer')
    call tracer_add('adv', dt, 'q2', 'Slotted cylinder tracer')
    call tracer_add('adv', dt, 'q3', 'Gaussian hill tracer')
    call tracer_add('adv', dt, 'q4', 'Constant one tracer')

  end subroutine solid_rotation_test_init

  subroutine solid_rotation_test_set_ic()

    integer iblk, i, j
    real(8) lon, lat, r, r0, qmin, qmax

    r0 = radius / 3.0_r8

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)              , &
                 mesh  => blocks(iblk)%mesh         , &
                 m     => blocks(iblk)%dstate(1)%dmg, &
                 q     => tracers(iblk)%q           )
      ! Background
      m%d(:,:,:) = 1
      ! Cosine bell
      do j = mesh%full_jds, mesh%full_jde
        lat = mesh%full_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          lon = mesh%full_lon(i)
          r = great_circle(radius, lon0, lat0, lon, lat)
          if (r < r0) then
            q%d(i,j,1,1) = h0 / 2.0_r8 * (1 + cos(pi * r / r0))
          else
            q%d(i,j,1,1) = 0
          end if
        end do
      end do
      call fill_halo(q, 1)
      ! Slotted cylinders
      qmax = 1.0_r8; qmin = 0.1_r8; r = radius * 0.5_r8
      do j = mesh%full_jds, mesh%full_jde
        lat = mesh%full_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          lon = mesh%full_lon(i)
          r0 = great_circle(radius, lon0, lat0, lon, lat)
          if (r0 <= r .and. abs(lon - lon0) >= r / radius / 6.0_r8) then
            q%d(i,j,1,2) = qmax
          else if (r0 <= r .and. abs(lon - lon0) < r / radius / 6.0_r8 .and. lat - lat0 < -5.0_r8 / 12.0_r8 * (r / radius)) then
            q%d(i,j,1,2) = qmax
          else
            q%d(i,j,1,2) = qmin
          end if
        end do
      end do
      call fill_halo(q, 2)
      ! Gaussian hill
      do j = mesh%full_jds, mesh%full_jde
        lat = mesh%full_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          lon = mesh%full_lon(i)
          r = great_circle(radius, lon0, lat0, lon, lat)
          q%d(i,j,1,3) = exp(-5 * r**2 / radius**2)
        end do
      end do
      call fill_halo(q, 3)
      ! Constant one
      q%d(:,:,:,4) = 1
      end associate
    end do

  end subroutine solid_rotation_test_set_ic

  subroutine solid_rotation_test_set_uv(time_in_seconds, itime)

    real(r8), intent(in) :: time_in_seconds
    integer, intent(in) :: itime

    integer iblk, i, j
    real(r8) lon, lat

    do iblk = 1, size(blocks)
      associate (block   => blocks(iblk)                    , &
                 mesh    => blocks(iblk)%mesh               , &
                 u       => blocks(iblk)%dstate(itime)%u_lon, &
                 v       => blocks(iblk)%dstate(itime)%v_lat)
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        lat = mesh%full_lat(j)
        do i = mesh%half_ids, mesh%half_ide
          lon = mesh%half_lon(i)
          u%d(i,j,1) = u0 * (cos(lat) * cos(alpha) + sin(lat) * cos(lon) * sin(alpha))
        end do
      end do
      call fill_halo(u)
      do j = mesh%half_jds, mesh%half_jde
        lat = mesh%half_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          lon = mesh%full_lon(i)
          v%d(i,j,1) = -u0 * sin(lon) * sin(alpha)
        end do
      end do
      call fill_halo(v)
      end associate
    end do

  end subroutine solid_rotation_test_set_uv

end module solid_rotation_test_mod
