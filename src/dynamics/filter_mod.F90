! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module filter_mod

  use const_mod
  use perf_mod
  use filter_types_mod
  use latlon_field_types_mod
  use latlon_parallel_mod

  implicit none

  private

  public filter_type
  public filter_run
  public filter_run_vector

  interface filter_run
    module procedure filter_run_2d
    module procedure filter_run_3d
    module procedure filter_run_4d
  end interface filter_run

contains

  subroutine filter_run_2d(filter, x, y)

    type(filter_type), intent(in), target :: filter
    type(latlon_field2d_type), intent(inout) :: x
    type(latlon_field2d_type), intent(inout), optional :: y

    real(r8) tmp(x%mesh%full_ids:x%mesh%full_ide)
    real(r8), pointer :: wgt(:,:)
    integer, pointer :: ngrid(:)
    integer is, ie, js, je, i, j, n, hn

    call perf_start('filter_run_2d')

    call fill_halo(x, south_halo=.false., north_halo=.false.)

    is = x%mesh%full_ids
    ie = x%mesh%full_ide
    select case (x%loc)
    case ('cell', 'lon')
      js = x%mesh%full_jds
      je = x%mesh%full_jde
      wgt => filter%wgt_lon
      ngrid => filter%ngrid_lon
    end select

    do j = js, je
      if (ngrid(j) >= 3) then
        n  = ngrid(j)
        hn = (n - 1) / 2
        do i = is, ie
          tmp(i) = sum(wgt(:n,j) * x%d(i-hn:i+hn,j))
        end do
        if (present(y)) then
          y%d(is:ie,j) = tmp(is:ie)
        else
          x%d(is:ie,j) = tmp(is:ie)
        end if
      else if (present(y)) then
        y%d(is:ie,j) = x%d(is:ie,j)
      end if
    end do

    call perf_stop('filter_run_2d')

  end subroutine filter_run_2d

  subroutine filter_run_3d(filter, x, y)

    type(filter_type), intent(in), target :: filter
    type(latlon_field3d_type), intent(inout) :: x
    type(latlon_field3d_type), intent(inout), optional :: y

    real(r8) tmp(x%mesh%full_ids:x%mesh%full_ide)
    real(r8), pointer :: wgt(:,:)
    integer, pointer :: ngrid(:)
    integer is, ie, js, je, ks, ke, i, j, k, n, hn

    call perf_start('filter_run_3d')

    call fill_halo(x, south_halo=.false., north_halo=.false.)

    is = x%mesh%full_ids
    ie = x%mesh%full_ide
    select case (x%loc)
    case ('cell', 'lon')
      js = x%mesh%full_jds
      je = x%mesh%full_jde
      ks = x%mesh%full_kds
      ke = x%mesh%full_kde
      wgt => filter%wgt_lon
      ngrid => filter%ngrid_lon
    case ('lat', 'vtx')
      js = x%mesh%half_jds
      je = x%mesh%half_jde
      ks = x%mesh%full_kds
      ke = x%mesh%full_kde
      wgt => filter%wgt_lat
      ngrid => filter%ngrid_lat
    case ('lev')
      js = x%mesh%full_jds
      je = x%mesh%full_jde
      ks = x%mesh%half_kds
      ke = x%mesh%half_kde
      wgt => filter%wgt_lon
      ngrid => filter%ngrid_lon
    end select

    do k = ks, ke
      do j = js, je
        if (ngrid(j) >= 3) then
          n  = ngrid(j)
          hn = (n - 1) / 2
          do i = is, ie
            tmp(i) = sum(wgt(:n,j) * x%d(i-hn:i+hn,j,k))
          end do
          if (present(y)) then
            y%d(is:ie,j,k) = tmp(is:ie)
          else
            x%d(is:ie,j,k) = tmp(is:ie)
          end if
        else if (present(y)) then
          y%d(is:ie,j,k) = x%d(is:ie,j,k)
        end if
      end do
    end do

    call perf_stop('filter_run_3d')

  end subroutine filter_run_3d

  subroutine filter_run_4d(filter, x, i4, y)

    type(filter_type), intent(in), target :: filter
    type(latlon_field4d_type), intent(inout) :: x
    integer, intent(in) :: i4
    type(latlon_field4d_type), intent(inout), optional :: y

    real(r8) tmp(x%mesh%full_ids:x%mesh%full_ide)
    real(r8), pointer :: wgt(:,:)
    integer, pointer :: ngrid(:)
    integer is, ie, js, je, ks, ke, i, j, k, n, hn

    call perf_start('filter_run_4d')

    is = x%mesh%full_ids
    ie = x%mesh%full_ide
    select case (x%loc)
    case ('cell', 'lon')
      js = x%mesh%full_jds
      je = x%mesh%full_jde
      ks = x%mesh%full_kds
      ke = x%mesh%full_kde
      wgt => filter%wgt_lon
      ngrid => filter%ngrid_lon
    end select

    do k = ks, ke
      do j = js, je
        if (ngrid(j) >= 3) then
          n  = ngrid(j)
          hn = (n - 1) / 2
          do i = is, ie
            tmp(i) = sum(wgt(:n,j) * x%d(i-hn:i+hn,j,k,i4))
          end do
          if (present(y)) then
            y%d(is:ie,j,k,i4) = tmp(is:ie)
          else
            x%d(is:ie,j,k,i4) = tmp(is:ie)
          end if
        else if (present(y)) then
          y%d(is:ie,j,k,i4) = x%d(is:ie,j,k,i4)
        end if
      end do
    end do

    call perf_stop('filter_run_4d')

  end subroutine filter_run_4d

  subroutine filter_run_vector(filter, x_lon, y_lat, x_lon_save, y_lat_save)

    type(filter_type        ), intent(in   ) :: filter
    type(latlon_field3d_type), intent(inout) :: x_lon
    type(latlon_field3d_type), intent(inout) :: y_lat
    type(latlon_field3d_type), intent(inout) :: x_lon_save
    type(latlon_field3d_type), intent(inout) :: y_lat_save

    real(r8) xs (filter%mesh%full_ims:filter%mesh%full_ime)
    real(r8) ys (filter%mesh%full_ims:filter%mesh%full_ime)
    real(r8) tmp(filter%mesh%full_ids:filter%mesh%full_ide)
    real(r8) s, y_lon, x_lat
    integer i, j, k, n, hn

    call fill_halo(x_lon, south_halo=.false.)
    call fill_halo(y_lat, north_halo=.false.)

    x_lon_save%d = x_lon%d
    y_lat_save%d = y_lat%d

    associate (mesh => filter%mesh)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        if (filter%ngrid_lon(j) >= 3) then
          n  = filter%ngrid_lon(j)
          hn = (n - 1) / 2
          if (abs(mesh%full_lat(j)) < 70) then
            do i = mesh%half_ids, mesh%half_ide
              tmp(i) = sum(filter%wgt_lon(:n,j) * x_lon%d(i-hn:i+hn,j,k))
            end do
            x_lon%d(mesh%half_ids:mesh%half_ide,j,k) = tmp
          else
            s = sign(1.0_r8, mesh%full_lat(j))
            ! Transform onto polar plane.
            do i = mesh%half_ids - hn, mesh%half_ide + hn
              y_lon = mesh%tg_wgt_lon(1,j) * (y_lat_save%d(i,j-1,k) + y_lat_save%d(i+1,j-1,k)) + &
                      mesh%tg_wgt_lon(2,j) * (y_lat_save%d(i,j  ,k) + y_lat_save%d(i+1,j  ,k))
              xs(i) = s * (-x_lon_save%d(i,j,k) * mesh%half_sin_lon(i) / mesh%full_sin_lat(j) - y_lon * mesh%half_cos_lon(i) / mesh%full_sin_lat(j)**2)
              ys(i) = s * ( x_lon_save%d(i,j,k) * mesh%half_cos_lon(i) / mesh%full_sin_lat(j) - y_lon * mesh%half_sin_lon(i) / mesh%full_sin_lat(j)**2)
            end do
            do i = mesh%half_ids, mesh%half_ide
              tmp(i) = sum(filter%wgt_lon(:n,j) * xs(i-hn:i+hn))
            end do
            xs(mesh%half_ids:mesh%half_ide) = tmp
            do i = mesh%half_ids, mesh%half_ide
              tmp(i) = sum(filter%wgt_lon(:n,j) * ys(i-hn:i+hn))
            end do
            ys(mesh%half_ids:mesh%half_ide) = tmp
            ! Transform back.
            do i = mesh%half_ids, mesh%half_ide
              x_lon%d(i,j,k) = -s * mesh%full_sin_lat(j) * (mesh%half_sin_lon(i) * xs(i) - mesh%half_cos_lon(i) * ys(i))
            end do
          end if
        end if
      end do
      do j = mesh%half_jds, mesh%half_jde
        if (filter%ngrid_lat(j) >= 3) then
          n  = filter%ngrid_lat(j)
          hn = (n - 1) / 2
          if (abs(mesh%half_lat(j)) < 70) then
            do i = mesh%full_ids, mesh%full_ide
              tmp(i) = sum(filter%wgt_lat(:n,j) * y_lat%d(i-hn:i+hn,j,k))
            end do
            y_lat%d(mesh%full_ids:mesh%full_ide,j,k) = tmp
          else
            s = sign(1.0_r8, mesh%half_lat(j))
            ! Transform onto polar plane.
            do i = mesh%full_ids - hn, mesh%full_ide + hn
              x_lat = mesh%tg_wgt_lat(1,j) * (x_lon_save%d(i-1,j  ,k) + x_lon_save%d(i,j  ,k)) + &
                      mesh%tg_wgt_lat(2,j) * (x_lon_save%d(i-1,j+1,k) + x_lon_save%d(i,j+1,k))
              xs(i) = s * (-x_lat * mesh%full_sin_lon(i) / mesh%half_sin_lat(j) - y_lat_save%d(i,j,k) * mesh%full_cos_lon(i) / mesh%half_sin_lat(j)**2)
              ys(i) = s * ( x_lat * mesh%full_cos_lon(i) / mesh%half_sin_lat(j) - y_lat_save%d(i,j,k) * mesh%full_sin_lon(i) / mesh%half_sin_lat(j)**2)
            end do
            do i = mesh%full_ids, mesh%full_ide
              tmp(i) = sum(filter%wgt_lat(:n,j) * xs(i-hn:i+hn))
            end do
            xs(mesh%full_ids:mesh%full_ide) = tmp
            do i = mesh%full_ids, mesh%full_ide
              tmp(i) = sum(filter%wgt_lat(:n,j) * ys(i-hn:i+hn))
            end do
            ys(mesh%full_ids:mesh%full_ide) = tmp
            ! Transform back.
            do i = mesh%full_ids, mesh%full_ide
              y_lat%d(i,j,k) = -s * mesh%half_sin_lat(j)**2 * (mesh%full_cos_lon(i) * xs(i) + mesh%full_sin_lon(i) * ys(i))
            end do
          end if
        end if
      end do
    end do
    end associate

  end subroutine filter_run_vector

end module filter_mod
