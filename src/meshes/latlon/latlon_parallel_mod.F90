! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module latlon_parallel_mod

  use mpi
  use flogger
  use const_mod
  use perf_mod
  use latlon_mesh_mod
  use latlon_halo_mod
  use latlon_parallel_types_mod
  use latlon_parallel_zonal_mod
  use latlon_parallel_global_mod
  use latlon_field_types_mod

  implicit none

  private

  public proc
  public fill_halo
  public zonal_sum
  public zonal_max
  public zonal_avg
  public global_sum
  public global_max
  public global_min

  interface fill_halo
    module procedure fill_halo_2d
    module procedure fill_halo_3d
    module procedure fill_halo_4d
  end interface fill_halo

contains

  subroutine fill_halo_2d(field, west_halo, east_halo, south_halo, north_halo)

    type(latlon_field2d_type), intent(inout) :: field
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    logical west_halo_opt, east_halo_opt, south_halo_opt, north_halo_opt
    integer t1, t2, ierr

    call perf_start('fill_halo_2d')

    west_halo_opt  = .true. ; if (present(west_halo )) west_halo_opt  = west_halo
    east_halo_opt  = .true. ; if (present(east_halo )) east_halo_opt  = east_halo
    south_halo_opt = .true. ; if (present(south_halo)) south_halo_opt = south_halo
    north_halo_opt = .true. ; if (present(north_halo)) north_halo_opt = north_halo

    t1 = merge(1, 2, field%full_lon)
    t2 = merge(1, 2, field%full_lat)

    if (west_halo_opt) then
      call MPI_ISEND(field%d, 1, field%halo(east)%send_type_2d(t1,t2), field%halo(east)%proc_id, 21, &
                     proc%comm_model, field%send_req(east), ierr)
      call MPI_IRECV(field%d, 1, field%halo(west)%recv_type_2d(t1,t2), field%halo(west)%proc_id, 21, &
                     proc%comm_model, field%recv_req(west), ierr)
    end if

    if (east_halo_opt) then
      call MPI_ISEND(field%d, 1, field%halo(west)%send_type_2d(t1,t2), field%halo(west)%proc_id, 22, &
                     proc%comm_model, field%send_req(west), ierr)
      call MPI_IRECV(field%d, 1, field%halo(east)%recv_type_2d(t1,t2), field%halo(east)%proc_id, 22, &
                     proc%comm_model, field%recv_req(east), ierr)
    end if

    if (south_halo_opt) then
      if (.not. proc%at_north_pole) then
        call MPI_ISEND(field%d, 1, field%halo(north)%send_type_2d(t1,t2), field%halo(north)%proc_id, 23, &
                       proc%comm_model, field%send_req(north), ierr)
      end if
      if (proc%at_south_pole) then
        call MPI_ISEND(field%d, 1, field%halo(south)%send_type_2d(t1,t2), field%halo(south)%proc_id, 23, &
                       proc%comm_model, field%send_req(south), ierr)
      end if
      call MPI_IRECV(field%d, 1, field%halo(south)%recv_type_2d(t1,t2), field%halo(south)%proc_id, 23, &
                     proc%comm_model, field%recv_req(south), ierr)
    end if

    if (north_halo_opt) then
      if (.not. proc%at_south_pole) then
        call MPI_ISEND(field%d, 1, field%halo(south)%send_type_2d(t1,t2), field%halo(south)%proc_id, 24, &
                       proc%comm_model, field%send_req(south), ierr)
      end if
      if (proc%at_north_pole) then
        call MPI_ISEND(field%d, 1, field%halo(north)%send_type_2d(t1,t2), field%halo(north)%proc_id, 24, &
                       proc%comm_model, field%send_req(north), ierr)
      end if
      call MPI_IRECV(field%d, 1, field%halo(north)%recv_type_2d(t1,t2), field%halo(north)%proc_id, 24, &
                     proc%comm_model, field%recv_req(north), ierr)
    end if

    if (field%halo_diagonal .and. field%halo(south_west)%initialized) then
      if (.not. proc%at_north_pole) then
        call MPI_ISEND(field%d, 1, field%halo(north_east)%send_type_2d(t1,t2), field%halo(north_east)%proc_id, 25, &
                       proc%comm_model, field%send_req(north_east), ierr)
      end if
      if (.not. proc%at_south_pole) then
        call MPI_IRECV(field%d, 1, field%halo(south_west)%recv_type_2d(t1,t2), field%halo(south_west)%proc_id, 25, &
                       proc%comm_model, field%recv_req(south_west), ierr)
      end if

      if (.not. proc%at_north_pole) then
        call MPI_ISEND(field%d, 1, field%halo(north_west)%send_type_2d(t1,t2), field%halo(north_west)%proc_id, 26, &
                       proc%comm_model, field%send_req(north_west), ierr)
      end if
      if (.not. proc%at_south_pole) then
        call MPI_IRECV(field%d, 1, field%halo(south_east)%recv_type_2d(t1,t2), field%halo(south_east)%proc_id, 26, &
                       proc%comm_model, field%recv_req(south_east), ierr)
      end if
    end if

    if (field%halo_diagonal .and. field%halo(north_west)%initialized) then
      if (.not. proc%at_south_pole) then
        call MPI_ISEND(field%d, 1, field%halo(south_east)%send_type_2d(t1,t2), field%halo(south_east)%proc_id, 27, &
                       proc%comm_model, field%send_req(south_east), ierr)
      end if
      if (.not. proc%at_north_pole) then
        call MPI_IRECV(field%d, 1, field%halo(north_west)%recv_type_2d(t1,t2), field%halo(north_west)%proc_id, 27, &
                       proc%comm_model, field%recv_req(north_west), ierr)
      end if

      if (.not. proc%at_south_pole) then
        call MPI_ISEND(field%d, 1, field%halo(south_west)%send_type_2d(t1,t2), field%halo(south_west)%proc_id, 28, &
                       proc%comm_model, field%send_req(south_west), ierr)
      end if
      if (.not. proc%at_north_pole) then
        call MPI_IRECV(field%d, 1, field%halo(north_east)%recv_type_2d(t1,t2), field%halo(north_east)%proc_id, 28, &
                       proc%comm_model, field%recv_req(north_east), ierr)
      end if
    end if

    call perf_stop('fill_halo_2d')

    call wait_halo_2d(field)

  end subroutine fill_halo_2d

  subroutine wait_halo_2d(field)

    type(latlon_field2d_type), intent(inout) :: field

    integer j, js, je, nx, mx, hx, hy, ierr
    real(r8) tmp(size(field%d,1),field%halo(1)%lat_hw)

    call perf_start('wait_halo_2d')

    hx = field%halo(1)%lon_hw
    hy = field%halo(1)%lat_hw
    if (field%full_lon) then
      nx = field%mesh%full_nlon
      mx = field%mesh%full_nlon / 2
    else
      nx = field%mesh%half_nlon
      mx = field%mesh%half_nlon / 2
    end if
    if (field%full_lat) then
      js = field%mesh%full_jms
      je = field%mesh%full_jme
    else
      js = field%mesh%half_jms
      je = field%mesh%half_jme
    end if

    if (field%recv_req(south) /= MPI_REQUEST_NULL) then
      call MPI_WAIT(field%recv_req(south), MPI_STATUS_IGNORE, ierr)
      if (proc%at_south_pole .and. field%halo_cross_pole) then
        ! Reverse array order.
        tmp = field%d(:,js:js+hy-1)
        if (field%halo(south)%proc_id == proc%id_model) then ! 1D decompostion, also reverse in lon
          do j = js, js + hy - 1
            field%d(   1:mx,j) = tmp(hx+1+mx:hx+nx,hy+js-j)
            field%d(mx+1:nx,j) = tmp(hx+1   :hx+mx,hy+js-j)
          end do
        else
          do j = js, js + hy - 1
            field%d(:,j) = tmp(:,hy+js-j)
          end do
        end if
      end if
      field%recv_req(south) = MPI_REQUEST_NULL
    end if

    if (field%recv_req(north) /= MPI_REQUEST_NULL) then
      call MPI_WAIT(field%recv_req(north), MPI_STATUS_IGNORE, ierr)
      if (proc%at_north_pole .and. field%halo_cross_pole) then
        ! Reverse array order.
        tmp = field%d(:,je-hy+1:je)
        if (field%halo(north)%proc_id == proc%id_model) then ! 1D decompostion, also reverse in lon
          do j = je - hy + 1, je
            field%d(   1:mx,j) = tmp(hx+1+mx:hx+nx,je+1-j)
            field%d(mx+1:nx,j) = tmp(hx+1   :hx+mx,je+1-j)
          end do
        else
          do j = je - hy + 1, je
            field%d(:,j) = tmp(:,je+1-j)
          end do
        end if
      end if
      field%recv_req(north) = MPI_REQUEST_NULL
    end if

    do j = 1, size(field%send_req)
      call MPI_WAIT(field%send_req(j), MPI_STATUS_IGNORE, ierr); field%send_req(j) = MPI_REQUEST_NULL
      call MPI_WAIT(field%recv_req(j), MPI_STATUS_IGNORE, ierr); field%recv_req(j) = MPI_REQUEST_NULL
    end do

    call perf_stop('wait_halo_2d')

  end subroutine wait_halo_2d

  subroutine fill_halo_3d(field, west_halo, east_halo, south_halo, north_halo)

    type(latlon_field3d_type), intent(inout) :: field
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    logical west_halo_opt, east_halo_opt, south_halo_opt, north_halo_opt
    integer t1, t2, t3, ierr

    call perf_start('fill_halo_3d')

    west_halo_opt  = .true. ; if (present(west_halo )) west_halo_opt  = west_halo
    east_halo_opt  = .true. ; if (present(east_halo )) east_halo_opt  = east_halo
    south_halo_opt = .true. ; if (present(south_halo)) south_halo_opt = south_halo
    north_halo_opt = .true. ; if (present(north_halo)) north_halo_opt = north_halo

    t1 = merge(1, 2, field%full_lon)
    t2 = merge(1, 2, field%full_lat)
    t3 = merge(1, 2, field%full_lev)

    if (west_halo_opt) then
      call MPI_ISEND(field%d, 1, field%halo(east)%send_type_3d(t1,t2,t3), field%halo(east)%proc_id, 31, &
                     proc%comm_model, field%send_req(east), ierr)
      call MPI_IRECV(field%d, 1, field%halo(west)%recv_type_3d(t1,t2,t3), field%halo(west)%proc_id, 31, &
                     proc%comm_model, field%recv_req(west), ierr)
    end if

    if (east_halo_opt) then
      call MPI_ISEND(field%d, 1, field%halo(west)%send_type_3d(t1,t2,t3), field%halo(west)%proc_id, 32, &
                     proc%comm_model, field%send_req(west), ierr)
      call MPI_IRECV(field%d, 1, field%halo(east)%recv_type_3d(t1,t2,t3), field%halo(east)%proc_id, 32, &
                     proc%comm_model, field%recv_req(east), ierr)
    end if

    if (south_halo_opt) then
      if (.not. proc%at_north_pole) then
        call MPI_ISEND(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 33, &
                       proc%comm_model, field%send_req(north), ierr)
      end if
      if (proc%at_south_pole) then
        call MPI_ISEND(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 33, &
                       proc%comm_model, field%send_req(south), ierr)
      end if
      call MPI_IRECV(field%d, 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 33, &
                     proc%comm_model, field%recv_req(south), ierr)
    end if

    if (north_halo_opt) then
      if (.not. proc%at_south_pole) then
        call MPI_ISEND(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 34, &
                       proc%comm_model, field%send_req(south), ierr)
      end if
      if (proc%at_north_pole) then
        call MPI_ISEND(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 34, &
                       proc%comm_model, field%send_req(north), ierr)
      end if
      call MPI_IRECV(field%d, 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 34, &
                     proc%comm_model, field%recv_req(north), ierr)
    end if

    if (field%halo_diagonal .and. field%halo(south_west)%initialized) then
      if (.not. proc%at_north_pole) then
        call MPI_ISEND(field%d, 1, field%halo(north_east)%send_type_3d(t1,t2,t3), field%halo(north_east)%proc_id, 37, &
                       proc%comm_model, field%send_req(north_east), ierr)
      end if
      if (.not. proc%at_south_pole) then
        call MPI_IRECV(field%d, 1, field%halo(south_west)%recv_type_3d(t1,t2,t3), field%halo(south_west)%proc_id, 37, &
                       proc%comm_model, field%recv_req(south_west), ierr)
      end if

      if (.not. proc%at_north_pole) then
        call MPI_ISEND(field%d, 1, field%halo(north_west)%send_type_3d(t1,t2,t3), field%halo(north_west)%proc_id, 38, &
                       proc%comm_model, field%send_req(north_west), ierr)
      end if
      if (.not. proc%at_south_pole) then
        call MPI_IRECV(field%d, 1, field%halo(south_east)%recv_type_3d(t1,t2,t3), field%halo(south_east)%proc_id, 38, &
                       proc%comm_model, field%recv_req(south_east), ierr)
      end if
    end if

    if (field%halo_diagonal .and. field%halo(north_west)%initialized) then
      if (.not. proc%at_south_pole) then
        call MPI_ISEND(field%d, 1, field%halo(south_east)%send_type_3d(t1,t2,t3), field%halo(south_east)%proc_id, 39, &
                       proc%comm_model, field%send_req(south_east), ierr)
      end if
      if (.not. proc%at_north_pole) then
        call MPI_IRECV(field%d, 1, field%halo(north_west)%recv_type_3d(t1,t2,t3), field%halo(north_west)%proc_id, 39, &
                       proc%comm_model, field%recv_req(north_west), ierr)
      end if

      if (.not. proc%at_south_pole) then
        call MPI_ISEND(field%d, 1, field%halo(south_west)%send_type_3d(t1,t2,t3), field%halo(south_west)%proc_id, 40, &
                       proc%comm_model, field%send_req(south_west), ierr)
      end if
      if (.not. proc%at_north_pole) then
        call MPI_IRECV(field%d, 1, field%halo(north_east)%recv_type_3d(t1,t2,t3), field%halo(north_east)%proc_id, 40, &
                       proc%comm_model, field%recv_req(north_east), ierr)
      end if
    end if

    call perf_stop('fill_halo_3d')

    call wait_halo_3d(field)

  end subroutine fill_halo_3d

  subroutine wait_halo_3d(field)

    type(latlon_field3d_type), intent(inout) :: field

    integer j, js, je, nx, mx, hx, hy, ierr
    real(r8) tmp(size(field%d,1),field%halo(1)%lat_hw,size(field%d,3))

    call perf_start('wait_halo_3d')

    hx = field%halo(1)%lon_hw
    hy = field%halo(1)%lat_hw
    if (field%full_lon) then
      nx = field%mesh%full_nlon
      mx = field%mesh%full_nlon / 2
    else
      nx = field%mesh%half_nlon
      mx = field%mesh%half_nlon / 2
    end if
    if (field%full_lat) then
      js = field%mesh%full_jms
      je = field%mesh%full_jme
    else
      js = field%mesh%half_jms
      je = field%mesh%half_jme
    end if

    if (field%recv_req(south) /= MPI_REQUEST_NULL) then
      call MPI_WAIT(field%recv_req(south), MPI_STATUS_IGNORE, ierr)
      if (proc%at_south_pole .and. field%halo_cross_pole) then
        ! Reverse array order.
        tmp = field%d(:,js:js+hy-1,:)
        if (field%halo(south)%proc_id == proc%id_model) then ! 1D decompostion, also reverse in lon
          do j = js, js + hy - 1
            field%d(   1:mx,j,:) = tmp(hx+1+mx:hx+nx,hy+js-j,:)
            field%d(mx+1:nx,j,:) = tmp(hx+1   :hx+mx,hy+js-j,:)
          end do
        else
          do j = js, js + hy - 1
            field%d(:,j,:) = tmp(:,hy+js-j,:)
          end do
        end if
      end if
      field%recv_req(south) = MPI_REQUEST_NULL
    end if

    if (field%recv_req(north) /= MPI_REQUEST_NULL) then
      call MPI_WAIT(field%recv_req(north), MPI_STATUS_IGNORE, ierr)
      if (proc%at_north_pole .and. field%halo_cross_pole) then
        ! Reverse array order.
        tmp = field%d(:,je-hy+1:je,:)
        if (field%halo(north)%proc_id == proc%id_model) then ! 1D decompostion, also reverse in lon
          do j = je - hy + 1, je
            field%d(   1:mx,j,:) = tmp(hx+1+mx:hx+nx,je+1-j,:)
            field%d(mx+1:nx,j,:) = tmp(hx+1   :hx+mx,je+1-j,:)
          end do
        else
          do j = je - hy + 1, je
            field%d(:,j,:) = tmp(:,je+1-j,:)
          end do
        end if
      end if
      field%recv_req(north) = MPI_REQUEST_NULL
    end if

    do j = 1, size(field%send_req)
      call MPI_WAIT(field%send_req(j), MPI_STATUS_IGNORE, ierr); field%send_req(j) = MPI_REQUEST_NULL
      call MPI_WAIT(field%recv_req(j), MPI_STATUS_IGNORE, ierr); field%recv_req(j) = MPI_REQUEST_NULL
    end do

    call perf_stop('wait_halo_3d')

  end subroutine wait_halo_3d

  subroutine fill_halo_4d(field, i4, west_halo, east_halo, south_halo, north_halo)

    type(latlon_field4d_type), intent(inout) :: field
    integer, intent(in) :: i4
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    logical west_halo_opt, east_halo_opt, south_halo_opt, north_halo_opt
    integer t1, t2, t3, ierr

    call perf_start('fill_halo_4d')

    west_halo_opt  = .true. ; if (present(west_halo )) west_halo_opt  = west_halo
    east_halo_opt  = .true. ; if (present(east_halo )) east_halo_opt  = east_halo
    south_halo_opt = .true. ; if (present(south_halo)) south_halo_opt = south_halo
    north_halo_opt = .true. ; if (present(north_halo)) north_halo_opt = north_halo

    t1 = merge(1, 2, field%full_lon)
    t2 = merge(1, 2, field%full_lat)
    t3 = merge(1, 2, field%full_lev)

    if (west_halo_opt) then
      call MPI_ISEND(field%d(:,:,:,i4), 1, field%halo(east)%send_type_3d(t1,t2,t3), field%halo(east)%proc_id, 41, &
                     proc%comm_model, field%send_req(east,i4), ierr)
      call MPI_IRECV(field%d(:,:,:,i4), 1, field%halo(west)%recv_type_3d(t1,t2,t3), field%halo(west)%proc_id, 41, &
                     proc%comm_model, field%recv_req(west,i4), ierr)
    end if

    if (east_halo_opt) then
      call MPI_ISEND(field%d(:,:,:,i4), 1, field%halo(west)%send_type_3d(t1,t2,t3), field%halo(west)%proc_id, 42, &
                     proc%comm_model, field%send_req(west,i4), ierr)
      call MPI_IRECV(field%d(:,:,:,i4), 1, field%halo(east)%recv_type_3d(t1,t2,t3), field%halo(east)%proc_id, 42, &
                     proc%comm_model, field%recv_req(east,i4), ierr)
    end if

    if (south_halo_opt) then
      if (.not. proc%at_north_pole) then
        call MPI_ISEND(field%d(:,:,:,i4), 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 43, &
                       proc%comm_model, field%send_req(north,i4), ierr)
      end if
      if (proc%at_south_pole) then
        call MPI_ISEND(field%d(:,:,:,i4), 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 43, &
                       proc%comm_model, field%send_req(south,i4), ierr)
      end if
      call MPI_IRECV(field%d(:,:,:,i4), 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 43, &
                     proc%comm_model, field%recv_req(south,i4), ierr)
    end if

    if (north_halo_opt) then
      if (.not. proc%at_south_pole) then
        call MPI_ISEND(field%d(:,:,:,i4), 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 44, &
                       proc%comm_model, field%send_req(south,i4), ierr)
      end if
      if (proc%at_north_pole) then
        call MPI_ISEND(field%d(:,:,:,i4), 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 44, &
                       proc%comm_model, field%send_req(north,i4), ierr)
      end if
      call MPI_IRECV(field%d(:,:,:,i4), 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 44, &
                     proc%comm_model, field%recv_req(north,i4), ierr)
    end if

    if (field%halo_diagonal .and. field%halo(south_west)%initialized) then
      if (.not. proc%at_north_pole) then
        call MPI_ISEND(field%d(:,:,:,i4), 1, field%halo(north_east)%send_type_3d(t1,t2,t3), field%halo(north_east)%proc_id, 37, &
                       proc%comm_model, field%send_req(north_east,i4), ierr)
      end if
      if (.not. proc%at_south_pole) then
        call MPI_IRECV(field%d(:,:,:,i4), 1, field%halo(south_west)%recv_type_3d(t1,t2,t3), field%halo(south_west)%proc_id, 37, &
                       proc%comm_model, field%recv_req(south_west,i4), ierr)
      end if

      if (.not. proc%at_north_pole) then
        call MPI_ISEND(field%d(:,:,:,i4), 1, field%halo(north_west)%send_type_3d(t1,t2,t3), field%halo(north_west)%proc_id, 38, &
                       proc%comm_model, field%send_req(north_west,i4), ierr)
      end if
      if (.not. proc%at_south_pole) then
        call MPI_IRECV(field%d(:,:,:,i4), 1, field%halo(south_east)%recv_type_3d(t1,t2,t3), field%halo(south_east)%proc_id, 38, &
                       proc%comm_model, field%recv_req(south_east,i4), ierr)
      end if
    end if

    if (field%halo_diagonal .and. field%halo(north_west)%initialized) then
      if (.not. proc%at_south_pole) then
        call MPI_ISEND(field%d(:,:,:,i4), 1, field%halo(south_east)%send_type_3d(t1,t2,t3), field%halo(south_east)%proc_id, 39, &
                       proc%comm_model, field%send_req(south_east,i4), ierr)
      end if
      if (.not. proc%at_north_pole) then
        call MPI_IRECV(field%d(:,:,:,i4), 1, field%halo(north_west)%recv_type_3d(t1,t2,t3), field%halo(north_west)%proc_id, 39, &
                       proc%comm_model, field%recv_req(north_west,i4), ierr)
      end if

      if (.not. proc%at_south_pole) then
        call MPI_ISEND(field%d(:,:,:,i4), 1, field%halo(south_west)%send_type_3d(t1,t2,t3), field%halo(south_west)%proc_id, 40, &
                       proc%comm_model, field%send_req(south_west,i4), ierr)
      end if
      if (.not. proc%at_north_pole) then
        call MPI_IRECV(field%d(:,:,:,i4), 1, field%halo(north_east)%recv_type_3d(t1,t2,t3), field%halo(north_east)%proc_id, 40, &
                       proc%comm_model, field%recv_req(north_east,i4), ierr)
      end if
    end if

    call perf_stop('fill_halo_4d')

    call wait_halo_4d(field, i4)

  end subroutine fill_halo_4d

  subroutine wait_halo_4d(field, i4)

    type(latlon_field4d_type), intent(inout) :: field
    integer, intent(in) :: i4

    integer j, js, je, nx, mx, hx, hy, ierr
    real(r8) tmp(size(field%d,1),field%halo(1)%lat_hw,size(field%d,3))

    call perf_start('wait_halo_4d')

    hx = field%halo(1)%lon_hw
    hy = field%halo(1)%lat_hw
    if (field%full_lon) then
      nx = field%mesh%full_nlon
      mx = field%mesh%full_nlon / 2
    else
      nx = field%mesh%half_nlon
      mx = field%mesh%half_nlon / 2
    end if
    if (field%full_lat) then
      js = field%mesh%full_jms
      je = field%mesh%full_jme
    else
      js = field%mesh%half_jms
      je = field%mesh%half_jme
    end if

    if (field%recv_req(south,i4) /= MPI_REQUEST_NULL) then
      call MPI_WAIT(field%recv_req(south,i4), MPI_STATUS_IGNORE, ierr)
      if (proc%at_south_pole .and. field%halo_cross_pole) then
        ! Reverse array order.
        tmp = field%d(:,js:js+hy-1,:,i4)
        if (field%halo(south)%proc_id == proc%id_model) then ! 1D decompostion, also reverse in lon
          do j = js, js + hy - 1
            field%d(   1:mx,j,:,i4) = tmp(hx+1+mx:hx+nx,hy+js-j,:)
            field%d(mx+1:nx,j,:,i4) = tmp(hx+1   :hx+mx,hy+js-j,:)
          end do
        else
          do j = js, js + hy - 1
            field%d(:,j,:,i4) = tmp(:,hy+js-j,:)
          end do
        end if
      end if
      field%recv_req(south,i4) = MPI_REQUEST_NULL
    end if

    if (field%recv_req(north,i4) /= MPI_REQUEST_NULL) then
      call MPI_WAIT(field%recv_req(north,i4), MPI_STATUS_IGNORE, ierr)
      if (proc%at_north_pole .and. field%halo_cross_pole) then
        ! Reverse array order.
        tmp = field%d(:,je-hy+1:je,:,i4)
        if (field%halo(north)%proc_id == proc%id_model) then ! 1D decompostion, also reverse in lon
          do j = je - hy + 1, je
            field%d(   1:mx,j,:,i4) = tmp(hx+1+mx:hx+nx,je+1-j,:)
            field%d(mx+1:nx,j,:,i4) = tmp(hx+1   :hx+mx,je+1-j,:)
          end do
        else
          do j = je - hy + 1, je
            field%d(:,j,:,i4) = tmp(:,je+1-j,:)
          end do
        end if
      end if
      field%recv_req(north,i4) = MPI_REQUEST_NULL
    end if

    call perf_stop('wait_halo_4d')

  end subroutine wait_halo_4d

end module latlon_parallel_mod
