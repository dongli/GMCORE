! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module latlon_parallel_zonal_mod

  use mpi
  use latlon_mesh_mod
  use latlon_parallel_types_mod
  use latlon_field_types_mod

  implicit none

  private

  public zonal_sum
  public zonal_max
  public zonal_avg

  interface zonal_sum
    module procedure zonal_sum_0d_i4
    module procedure zonal_sum_0d_r4
    module procedure zonal_sum_0d_r8
    module procedure zonal_sum_1d_r4
    module procedure zonal_sum_1d_r8
  end interface zonal_sum

  interface zonal_max
    module procedure zonal_max_r4
    module procedure zonal_max_r8
  end interface zonal_max

  interface zonal_avg
    module procedure zonal_avg_3d_r8
    module procedure zonal_avg_4d_r8
  end interface zonal_avg

  interface gather_zonal_array
    module procedure gather_zonal_array_1d_i4
    module procedure gather_zonal_array_1d_r4
    module procedure gather_zonal_array_1d_r8
    module procedure gather_zonal_array_2d_r4
    module procedure gather_zonal_array_2d_r8
  end interface gather_zonal_array

  interface scatter_zonal_array
    module procedure scatter_zonal_array_1d_r4
    module procedure scatter_zonal_array_1d_r8
    module procedure scatter_zonal_array_2d_r4
    module procedure scatter_zonal_array_2d_r8
  end interface scatter_zonal_array

contains

  subroutine zonal_sum_0d_i4(zonal_circle, work, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    integer, intent(in) :: work(:)
    integer, intent(out) :: value

#ifdef ENSURE_ORDER
    integer allvalue(global_mesh%full_nlon)
#endif
    integer ierr

#ifdef ENSURE_ORDER
    if (zonal_circle%np == 1) then
      value = sum(work)
    else
      call gather_zonal_array(zonal_circle, work, allvalue)
      if (zonal_circle%id == 0) value = sum(allvalue)
      call MPI_BCAST(value, 1, MPI_INT, 0, zonal_circle%comm, ierr)
    end if
#else
    call MPI_ALLREDUCE(sum(work), value, 1, MPI_INT, MPI_SUM, zonal_circle%comm, ierr)
#endif

  end subroutine zonal_sum_0d_i4

  subroutine zonal_sum_0d_r4(zonal_circle, work, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: work(:)
    real(4), intent(out) :: value

#ifdef ENSURE_ORDER
    real(4) allvalue(global_mesh%full_nlon)
#endif
    integer ierr

#ifdef ENSURE_ORDER
    if (zonal_circle%np == 1) then
      value = sum(work)
    else
      call gather_zonal_array(zonal_circle, work, allvalue)
      if (zonal_circle%id == 0) value = sum(allvalue)
      call MPI_BCAST(value, 1, MPI_REAL, 0, zonal_circle%comm, ierr)
    end if
#else
    call MPI_ALLREDUCE(sum(work), value, 1, MPI_REAL, MPI_SUM, zonal_circle%comm, ierr)
#endif

  end subroutine zonal_sum_0d_r4

  subroutine zonal_sum_0d_r8(zonal_circle, work, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: work(:)
    real(8), intent(out) :: value

#ifdef ENSURE_ORDER
    real(8) allvalue(global_mesh%full_nlon)
#endif
    integer ierr

#ifdef ENSURE_ORDER
    if (zonal_circle%np == 1) then
      value = sum(work)
    else
      call gather_zonal_array(zonal_circle, work, allvalue)
      if (zonal_circle%id == 0) value = sum(allvalue)
      call MPI_BCAST(value, 1, MPI_DOUBLE, 0, zonal_circle%comm, ierr)
    end if
#else
    call MPI_ALLREDUCE(sum(work), value, 1, MPI_DOUBLE, MPI_SUM, zonal_circle%comm, ierr)
#endif

  end subroutine zonal_sum_0d_r8

  subroutine zonal_sum_1d_r4(zonal_circle, work, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: work(:,:)
    real(4), intent(out) :: value(:)

#ifdef ENSURE_ORDER
    real(4) allvalue(global_mesh%full_nlon,size(value))
#endif
    integer ierr

#ifdef ENSURE_ORDER
    if (zonal_circle%np == 1) then
      value = sum(work, dim=1)
    else
      call gather_zonal_array(zonal_circle, work, allvalue)
      if (zonal_circle%id == 0) value = sum(allvalue, dim=1)
      call MPI_BCAST(value, size(value), MPI_REAL, 0, zonal_circle%comm, ierr)
    end if
#else
    call MPI_ALLREDUCE(sum(work, dim=1), value, size(value), MPI_REAL, MPI_SUM, zonal_circle%comm, ierr)
#endif

  end subroutine zonal_sum_1d_r4

  subroutine zonal_sum_1d_r8(zonal_circle, work, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: work(:,:)
    real(8), intent(out) :: value(:)

#ifdef ENSURE_ORDER
    real(8) allvalue(global_mesh%full_nlon,size(value))
#endif
    integer ierr

#ifdef ENSURE_ORDER
    if (zonal_circle%np == 1) then
      value = sum(work, dim=1)
    else
      call gather_zonal_array(zonal_circle, work, allvalue)
      if (zonal_circle%id == 0) value = sum(allvalue, dim=1)
      call MPI_BCAST(value, size(value), MPI_DOUBLE, 0, zonal_circle%comm, ierr)
    end if
#else
    call MPI_ALLREDUCE(sum(work, dim=1), value, size(value), MPI_DOUBLE, MPI_SUM, zonal_circle%comm, ierr)
#endif

  end subroutine zonal_sum_1d_r8

  subroutine zonal_max_r4(zonal_circle, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(inout) :: value

    integer ierr
    real(4) res

    call MPI_ALLREDUCE(value, res, 1, MPI_REAL, MPI_MAX, zonal_circle%comm, ierr)
    value = res

  end subroutine zonal_max_r4

  subroutine zonal_max_r8(zonal_circle, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(inout) :: value

    integer ierr
    real(8) res

    call MPI_ALLREDUCE(value, res, 1, MPI_DOUBLE, MPI_MAX, zonal_circle%comm, ierr)
    value = res

  end subroutine zonal_max_r8

  subroutine gather_zonal_array_1d_i4(zonal_circle, local_array, array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    integer, intent(in) :: local_array(:)
    integer, intent(out) :: array(:)

    integer ierr, i

    if (zonal_circle%id == 0) then
      array(1:size(local_array)) = local_array
      do i = 2, zonal_circle%np
        call MPI_RECV(array, 1, zonal_circle%recv_type_i4(i,0), i - 1, 30, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
      end do
    else
      call MPI_SEND(local_array, size(local_array), MPI_INT, 0, 30, zonal_circle%comm, ierr)
    end if

  end subroutine gather_zonal_array_1d_i4

  subroutine gather_zonal_array_1d_r4(zonal_circle, local_array, array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: local_array(:)
    real(4), intent(out) :: array(:)

    integer ierr, i

    if (zonal_circle%id == 0) then
      array(1:size(local_array)) = local_array
      do i = 2, zonal_circle%np
        call MPI_RECV(array, 1, zonal_circle%recv_type_r4(i,0), i - 1, 30, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
      end do
    else
      call MPI_SEND(local_array, size(local_array), MPI_REAL, 0, 30, zonal_circle%comm, ierr)
    end if

  end subroutine gather_zonal_array_1d_r4

  subroutine gather_zonal_array_1d_r8(zonal_circle, local_array, array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: local_array(:)
    real(8), intent(out) :: array(:)

    integer ierr, i

    if (zonal_circle%id == 0) then
      array(1:size(local_array)) = local_array
      do i = 2, zonal_circle%np
        call MPI_RECV(array, 1, zonal_circle%recv_type_r8(i,0), i - 1, 30, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
      end do
    else
      call MPI_SEND(local_array, size(local_array), MPI_DOUBLE, 0, 30, zonal_circle%comm, ierr)
    end if

  end subroutine gather_zonal_array_1d_r8

  subroutine gather_zonal_array_2d_r4(zonal_circle, local_array, array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: local_array(:,:)
    real(4), intent(inout) :: array(:,:)

    integer ierr, i, k

    if (zonal_circle%id == 0) then
      k = merge(1, 2, size(local_array, 2) == global_mesh%full_nlev)
      array(1:size(local_array, 1),:) = local_array
      do i = 2, zonal_circle%np
        call MPI_RECV(array, 1, zonal_circle%recv_type_r4(i,k), i - 1, 31, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
      end do
    else
      call MPI_SEND(local_array, size(local_array), MPI_REAL, 0, 31, zonal_circle%comm, ierr)
    end if

  end subroutine gather_zonal_array_2d_r4

  subroutine gather_zonal_array_2d_r8(zonal_circle, local_array, array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: local_array(:,:)
    real(8), intent(inout) :: array(:,:)

    integer ierr, i, k

    if (zonal_circle%id == 0) then
      k = merge(1, 2, size(local_array, 2) == global_mesh%full_nlev)
      array(1:size(local_array, 1),:) = local_array
      do i = 2, zonal_circle%np
        call MPI_RECV(array, 1, zonal_circle%recv_type_r8(i,k), i - 1, 31, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
      end do
    else
      call MPI_SEND(local_array, size(local_array), MPI_DOUBLE, 0, 31, zonal_circle%comm, ierr)
    end if

  end subroutine gather_zonal_array_2d_r8

  subroutine scatter_zonal_array_1d_r4(zonal_circle, array, local_array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: array(:)
    real(4), intent(out) :: local_array(:)

    integer ierr, i

    if (zonal_circle%id == 0) then
      local_array = array(1:size(local_array))
      do i = 2, zonal_circle%np
        call MPI_SEND(array, 1, zonal_circle%recv_type_r4(i,0), i - 1, 32, zonal_circle%comm, ierr)
      end do
    else
      call MPI_RECV(local_array, size(local_array), MPI_REAL, 0, 32, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
    end if

  end subroutine scatter_zonal_array_1d_r4

  subroutine scatter_zonal_array_1d_r8(zonal_circle, array, local_array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: array(:)
    real(8), intent(out) :: local_array(:)

    integer ierr, i

    if (zonal_circle%id == 0) then
      local_array = array(1:size(local_array))
      do i = 2, zonal_circle%np
        call MPI_SEND(array, 1, zonal_circle%recv_type_r8(i,0), i - 1, 32, zonal_circle%comm, ierr)
      end do
    else
      call MPI_RECV(local_array, size(local_array), MPI_DOUBLE, 0, 32, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
    end if

  end subroutine scatter_zonal_array_1d_r8

  subroutine scatter_zonal_array_2d_r4(zonal_circle, array, local_array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: array(:,:)
    real(4), intent(out) :: local_array(:,:)

    integer ierr, i, k

    if (zonal_circle%id == 0) then
      k = merge(1, 2, size(local_array, 2) == global_mesh%full_nlev)
      local_array = array(1:size(local_array, 1),:)
      do i = 2, zonal_circle%np
        call MPI_SEND(array, 1, zonal_circle%recv_type_r4(i,k), i - 1, 33, zonal_circle%comm, ierr)
      end do
    else
      call MPI_RECV(local_array, size(local_array), MPI_REAL, 0, 33, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
    end if

  end subroutine scatter_zonal_array_2d_r4

  subroutine scatter_zonal_array_2d_r8(zonal_circle, array, local_array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: array(:,:)
    real(8), intent(out) :: local_array(:,:)

    integer ierr, i, k

    if (zonal_circle%id == 0) then
      k = merge(1, 2, size(local_array, 2) == global_mesh%full_nlev)
      local_array = array(1:size(local_array, 1),:)
      do i = 2, zonal_circle%np
        call MPI_SEND(array, 1, zonal_circle%recv_type_r8(i,k), i - 1, 33, zonal_circle%comm, ierr)
      end do
    else
      call MPI_RECV(local_array, size(local_array), MPI_DOUBLE, 0, 33, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
    end if

  end subroutine scatter_zonal_array_2d_r8

  subroutine zonal_avg_3d_r8(zonal_circle, f, j)

    type(zonal_circle_type), intent(in) :: zonal_circle
    type(latlon_field3d_type), intent(inout) :: f
    integer, intent(in) :: j

    real(8) work(f%mesh%full_ids:f%mesh%full_ide,f%mesh%full_kds:f%mesh%full_kde)
    real(8) pole(f%mesh%full_kds:f%mesh%full_kde)
    integer i, k

    do k = f%mesh%full_kds, f%mesh%full_kde
      do i = f%mesh%full_ids, f%mesh%full_ide
        work(i,k) = f%d(i,j,k)
      end do
    end do
    call zonal_sum(zonal_circle, work, pole)
    pole = pole / global_mesh%full_nlon
    do k = f%mesh%full_kds, f%mesh%full_kde
      do i = f%mesh%full_ids, f%mesh%full_ide
        f%d(i,j,k) = pole(k)
      end do
    end do

  end subroutine zonal_avg_3d_r8

  subroutine zonal_avg_4d_r8(zonal_circle, f, j)

    type(zonal_circle_type), intent(in) :: zonal_circle
    type(latlon_field4d_type), intent(inout) :: f
    integer, intent(in) :: j

    real(8) work(f%mesh%full_ids:f%mesh%full_ide,f%mesh%full_kds:f%mesh%full_kde)
    real(8) pole(f%mesh%full_kds:f%mesh%full_kde)
    integer i, k, l

    do l = 1, f%dim4_size
      do k = f%mesh%full_kds, f%mesh%full_kde
        do i = f%mesh%full_ids, f%mesh%full_ide
          work(i,k) = f%d(i,j,k,l)
        end do
      end do
      call zonal_sum(zonal_circle, work, pole)
      pole = pole / global_mesh%full_nlon
      do k = f%mesh%full_kds, f%mesh%full_kde
        do i = f%mesh%full_ids, f%mesh%full_ide
          f%d(i,j,k,l) = pole(k)
        end do
      end do
    end do

  end subroutine zonal_avg_4d_r8

end module latlon_parallel_zonal_mod
