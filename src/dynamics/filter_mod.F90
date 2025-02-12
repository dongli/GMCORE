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
  use math_mod
  use filter_types_mod
  use latlon_field_types_mod
  use latlon_parallel_mod

  implicit none

  private

  public filter_type
  public filter_run

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

    if (.not. associated(x%mesh, filter%mesh)) then
      call log_error('x%mesh is not associated with filter%mesh', __FILE__, __LINE__, pid=proc%id_model)
    end if

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

    if (.not. associated(x%mesh, filter%mesh)) then
      call log_error('x%mesh is not associated with filter%mesh', __FILE__, __LINE__, pid=proc%id_model)
    end if

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

    if (.not. associated(x%mesh, filter%mesh)) then
      call log_error('x%mesh is not associated with filter%mesh', __FILE__, __LINE__, pid=proc%id_model)
    end if

    call fill_halo(x, i4, south_halo=.false., north_halo=.false.)

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

end module filter_mod
