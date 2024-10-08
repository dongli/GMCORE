! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================
! Description:
!
!   This module is a test case with cross pole wind field.
!
! References:
!
!   - MacDonald, A., J. R. Bates, 1989: Semi-Lagrangian Integration of a
!     Gridpoint Shallow Water Model on the Sphere. MWR, 117: 130-137.
!
! Authors:

!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
!   - Jianghao Li
! ==============================================================================


module cross_pole_flow_test_mod

  use flogger
  use const_mod
  use latlon_parallel_mod
  use block_mod

  implicit none

  private

  public cross_pole_flow_test_set_ic

  real, parameter :: v0 = 20 ! m s-1
  real, parameter :: gz0 = 5.7684e4 ! m2 s-2

contains

  subroutine cross_pole_flow_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    real(r8) cos_lat, sin_lat, cos_lon, sin_lon
    integer i, j

    associate (mesh   => block%mesh           , &
               u      => block%dstate(1)%u_lon, &
               v      => block%dstate(1)%v_lat, &
               gz     => block%dstate(1)%gz   , &
               gzs    => block%static%gzs)
    gzs%d = 0

    do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = mesh%half_ids, mesh%half_ide
        sin_lon = mesh%half_sin_lon(i)
        u%d(i,j,1) = -v0 * sin_lon * sin_lat * (4.0 * cos_lat**2 - 1.0)
      end do
    end do
    call fill_halo(u)

    do j = mesh%half_jds, mesh%half_jde
      sin_lat = mesh%half_sin_lat(j)
      do i = mesh%full_ids, mesh%full_ide
        cos_lon = mesh%full_cos_lon(i)
        v%d(i,j,1) = v0 * sin_lat**2 * cos_lon
      end do
    end do
    call fill_halo(v)

    do j = mesh%full_jds, mesh%full_jde
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = mesh%full_ids, mesh%full_ide
        sin_lon = mesh%full_sin_lon(i)
        gz%d(i,j,1) = gz0 + 2 * radius * omega * v0 * sin_lat**3 * cos_lat * sin_lon
      end do
    end do
    call fill_halo(gz)
    end associate

  end subroutine cross_pole_flow_test_set_ic

end module cross_pole_flow_test_mod
