module latlon_parallel_global_mod

  use mpi

  implicit none

  private

  public global_sum
  public global_max
  public global_min

  interface global_sum
    module procedure global_sum_0d_r4
    module procedure global_sum_0d_r8
  end interface global_sum

  interface global_max
    module procedure global_max_0d_r4
    module procedure global_max_0d_r8
    module procedure global_max_0d_i4
  end interface global_max

  interface global_min
    module procedure global_min_0d_r4
    module procedure global_min_0d_r8
    module procedure global_min_0d_i4
  end interface global_min

contains

  real(4) function global_sum_0d_r4(comm, value) result(res)

    integer, intent(in) :: comm
    real(4), intent(in) :: value

    integer ierr

    call MPI_ALLREDUCE(value, res, 1, MPI_REAL, MPI_SUM, comm, ierr)

  end function global_sum_0d_r4

  real(8) function global_sum_0d_r8(comm, value) result(res)

    integer, intent(in) :: comm
    real(8), intent(in) :: value

    integer ierr

    call MPI_ALLREDUCE(value, res, 1, MPI_DOUBLE, MPI_SUM, comm, ierr)

  end function global_sum_0d_r8

  real(4) function global_max_0d_r4(comm, value) result(res)

    integer, intent(in) :: comm
    real(4), intent(in) :: value

    integer ierr

    call MPI_ALLREDUCE(value, res, 1, MPI_REAL, MPI_MAX, comm, ierr)

  end function global_max_0d_r4

  real(8) function global_max_0d_r8(comm, value) result(res)

    integer, intent(in) :: comm
    real(8), intent(in) :: value

    integer ierr

    call MPI_ALLREDUCE(value, res, 1, MPI_DOUBLE, MPI_MAX, comm, ierr)

  end function global_max_0d_r8

  integer function global_max_0d_i4(comm, value) result(res)

    integer, intent(in) :: comm
    integer(4), intent(in) :: value

    integer ierr

    call MPI_ALLREDUCE(value, res, 1, MPI_INT, MPI_MAX, comm, ierr)

  end function global_max_0d_i4

  real(4) function global_min_0d_r4(comm, value) result(res)

    integer, intent(in) :: comm
    real(4), intent(in) :: value

    integer ierr

    call MPI_ALLREDUCE(value, res, 1, MPI_REAL, MPI_MIN, comm, ierr)

  end function global_min_0d_r4

  real(8) function global_min_0d_r8(comm, value) result(res)

    integer, intent(in) :: comm
    real(8), intent(in) :: value

    integer ierr

    call MPI_ALLREDUCE(value, res, 1, MPI_DOUBLE, MPI_MIN, comm, ierr)

  end function global_min_0d_r8

  integer function global_min_0d_i4(comm, value) result(res)

    integer, intent(in) :: comm
    integer(4), intent(in) :: value

    integer ierr

    call MPI_ALLREDUCE(value, res, 1, MPI_INT, MPI_MIN, comm, ierr)

  end function global_min_0d_i4

end module latlon_parallel_global_mod