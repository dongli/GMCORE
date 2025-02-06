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
!   This module setups the halo exchange pattern, and defines the MPI data types
!   for latter halo filling.
!
! Author:
!
!   - Li Dong <dongli@lasg.iap.ac.cn>
! ==============================================================================

module latlon_halo_mod

  use mpi
  use flogger
  use const_mod
  use latlon_mesh_mod
  use latlon_parallel_types_mod

  implicit none

  private

  public latlon_halo_type

  integer, parameter :: cross_proc_halo = 1
  integer, parameter :: cross_comm_halo = 2
  integer, parameter :: inner_halo = 3
  integer, parameter :: nest_halo = 4

  type latlon_halo_type
    logical :: initialized = .false.
    integer :: comm     = MPI_COMM_NULL
    integer :: host_id  = MPI_PROC_NULL
    integer :: proc_id  = MPI_PROC_NULL
    integer :: iblk     = 0
    integer :: orient   = 0
    integer :: dtype    = 0
    integer :: type     = 0
    integer :: lon_hw   = 0
    integer :: lat_hw   = 0
    ! (1,1): full_lon,full_lat (1,2): full_lon,half_lat
    ! (2,1): half_lon,full_lat (2,2): half_lon,half_lat
    integer :: send_type_2d(2,2) = MPI_DATATYPE_NULL
    integer :: recv_type_2d(2,2) = MPI_DATATYPE_NULL
    ! (1,1,1): full_lon,full_lat,full_lev (1,2,1): full_lon,half_lat,full_lev
    ! (2,1,1): half_lon,full_lat,full_lev (2,2,1): half_lon,half_lat,full_lev
    ! (1,1,2): full_lon,full_lat,half_lev (1,2,2): full_lon,half_lat,half_lev
    ! (2,1,2): half_lon,full_lat,half_lev (2,2,2): half_lon,half_lat,half_lev
    integer :: send_type_3d(2,2,2) = MPI_DATATYPE_NULL
    integer :: recv_type_3d(2,2,2) = MPI_DATATYPE_NULL
  contains
    procedure :: init => latlon_halo_init
    procedure :: init_nest => latlon_halo_init_nest
    procedure :: clear => latlon_halo_clear
    final latlon_halo_final
  end type latlon_halo_type

contains

  subroutine latlon_halo_init(this, mesh, orient, dtype, host_id, ngb_id, iblk, lon_hw, lat_hw)

    class(latlon_halo_type), intent(out) :: this
    type(latlon_mesh_type), intent(in) :: mesh
    integer, intent(in) :: orient
    integer, intent(in) :: dtype
    integer, intent(in) :: host_id
    integer, intent(in), optional :: ngb_id
    integer, intent(in), optional :: iblk
    integer, intent(in), optional :: lon_hw
    integer, intent(in), optional :: lat_hw

    integer full_ihs, full_ihe
    integer full_jhs, full_jhe
    integer half_ihs, half_ihe
    integer half_jhs, half_jhe
    integer array_size(3,2,2)
    integer send_subarray_size(3,2,2)
    integer recv_subarray_size(3,2,2)
    integer send_subarray_start(3,2,2)
    integer recv_subarray_start(3,2,2)
    integer nlev(2)
    integer i, j, k, ierr

    if (present(ngb_id)) then
      this%proc_id = ngb_id
    else if (present(iblk)) then
      call log_error('Handle internal halo!', __FILE__, __LINE__)
    end if

    this%host_id = host_id
    this%dtype = dtype

    this%lon_hw = mesh%lon_hw; if (present(lon_hw)) this%lon_hw = lon_hw
    this%lat_hw = mesh%lat_hw; if (present(lat_hw)) this%lat_hw = lat_hw
    ! Calculate the start and end indices of halo for MPI.
    ! NOTE: MPI array index starts from zero.
    select case (orient)
    case (south, north)
      full_ihs = mesh%lon_hw
      full_ihe = full_ihs + mesh%full_nlon - 1
    case (west)
      full_ihs = 0
      full_ihe = mesh%lon_hw - 1
    case (east)
      full_ihs = mesh%full_nlon + mesh%lon_hw
      full_ihe = full_ihs + mesh%lon_hw - 1
    case (south_west, north_west)
      ! NOTE: mesh%lon_hw may not equal to this%lon_hw.
      full_ihs = mesh%lon_hw - this%lon_hw
      full_ihe = full_ihs + this%lon_hw - 1
    case (south_east, north_east)
      full_ihs = mesh%full_nlon + mesh%lon_hw
      full_ihe = full_ihs + this%lon_hw - 1
    end select
    half_ihs = full_ihs
    half_ihe = full_ihe
    select case (orient)
    case (west, east)
      full_jhs = mesh%lat_hw
      full_jhe = full_jhs + mesh%full_nlat - 1
    case (south, south_west, south_east)
      ! NOET: mesh%lat_hw equals to this%lat_hw.
      full_jhs = 0
      full_jhe = mesh%lat_hw - 1
    case (north, north_west, north_east)
      full_jhs = mesh%full_nlat + mesh%lat_hw
      full_jhe = full_jhs + mesh%lat_hw - 1
    end select
    !
    !   South Pole           Internal           North Pole
    !
    !    U     V             U     V             U     V
    !                                          _____       
    !        _____               _____        |\\\\\|_____ 
    !       |\\\\\|             |\\\\\|       |\\6\\|\\\\\|
    !  _____|\\5\\|        _____|\\5\\|       |\\\\\|\\5\\|
    ! |\\\\\|\\\\\|       |\\\\\|\\\\\|       |\\\\\|\\\\\| 
    ! |\\5\\|\\\\\|       |\\5\\|\\\\\|       |\\5\\|\\\\\|
    ! |\\\\\|\\4\\|       |\\\\\|\\4\\|       |\\\\\|\\4\\|
    ! |\\\\\|\\\\\|       |\\\\\|\\\\\|       |@@@@@|\\\\\|
    ! |\\4\\|     |       |\\4\\|     |       |@@4@@|     |
    ! |\\\\\|  3  |       |\\\\\|  3  |       |@@@@@|  3  |
    ! |     |_____|       |     |_____|       |     |_____|
    ! |  3  |     |       |  3  |     |       |  3  |     |
    ! |_____|  2  |       |_____|  2  |       |_____|  2  |
    ! |@@@@@|_____|       |     |_____|       |     |_____|
    ! |@@2@@|\\\\\|       |  2  |\\\\\|       |  2  |\\\\\|
    ! |@@@@@|\\1\\|       |_____|\\1\\|       |_____|\\1\\|
    ! |\\\\\|\\\\\|       |\\\\\|\\\\\|       |\\\\\|\\\\\|
    ! |\\1\\|\\\\\|       |\\1\\|\\\\\|       |\\1\\|\\\\\|
    ! |\\\\\|\\0\\|       |\\\\\|\\0\\|       |\\\\\|\\0\\|
    ! |\\\\\|\\\\\|       |\\\\\|\\\\\|       |\\\\\|\\\\\|
    ! |\\0\\|             |\\0\\|             |\\0\\|      
    ! |\\\\\|             |\\\\\|             |\\\\\|      
    !                                         
    half_jhs = full_jhs + merge(-1, 0, mesh%has_north_pole() .and. (orient == north .or.  orient == north_west .or.  orient == north_east))
    half_jhe = full_jhe + merge(-1, 0, mesh%has_north_pole() .and. (orient /= south .and. orient /= south_west .and. orient /= south_east))

    !                          wx                          nx                          wx
    !                          |                           |                           |
    !                  |---------------|---------------------------------------|---------------|
    !                  |_______________|_______________________________________|_______________|__
    !                  |XXXXXXX|XXXXXXX|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|XXXXXXX|XXXXXXX|  |
    !         wy + ny -|XXXXXXX|XXXXXXX|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|XXXXXXX|XXXXXXX|  |- wy
    !                  |XXXXXXX|XXXXXXX|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|XXXXXXX|XXXXXXX|  |
    !                  |_______|_______|_______________________________________|_______________|__|
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !     wy + ny - 1 -|///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |_______|_______|_______|_______|_______|_______|_______|_______|_______|  |
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |- ny
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |_______|_______|_______|_______|_______|_______|_______|_______|_______|  |
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !              wy -|///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |_______|_______|_______|_______|_______|_______|_______|_______|_______|__|
    !                  |XXXXXXX|XXXXXXX|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|XXXXXXX|XXXXXXX|  |
    !               0 -|XXXXXXX|XXXXXXX|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|XXXXXXX|XXXXXXX|  |- wy
    !                  |XXXXXXX|XXXXXXX|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|XXXXXXX|XXXXXXX|  |
    !                  |_______|_______|_______|_______|_______|_______|_______|_______|_______|__|
    !                      |               |                       |       |       |
    !                      0               wx                      | wx + nx - 1   |
    !                                                              |            wx + nx
    !                                                              nx

    this%orient = orient
    this%type = cross_proc_halo
    nlev = [mesh%full_kme-mesh%full_kms+1,mesh%half_kme-mesh%half_kms+1]

    do k = 1, 2 ! From full level to half level
      array_size(:,1,1) = [mesh%full_nlon+2*mesh%lon_hw,mesh%full_nlat+2*mesh%lat_hw,nlev(k)]
      array_size(:,2,1) = [mesh%half_nlon+2*mesh%lon_hw,mesh%full_nlat+2*mesh%lat_hw,nlev(k)]
      array_size(:,1,2) = [mesh%full_nlon+2*mesh%lon_hw,mesh%half_nlat+2*mesh%lat_hw,nlev(k)]
      array_size(:,2,2) = [mesh%half_nlon+2*mesh%lon_hw,mesh%half_nlat+2*mesh%lat_hw,nlev(k)]
      send_subarray_size(:,1,1) = [full_ihe-full_ihs+1,full_jhe-full_jhs+1,nlev(k)]
      recv_subarray_size(:,1,1) = send_subarray_size(:,1,1)
      send_subarray_size(:,2,1) = [half_ihe-half_ihs+1,full_jhe-full_jhs+1,nlev(k)]
      recv_subarray_size(:,2,1) = send_subarray_size(:,2,1)
      send_subarray_size(:,1,2) = [full_ihe-full_ihs+1,half_jhe-half_jhs+1,nlev(k)]
      recv_subarray_size(:,1,2) = send_subarray_size(:,1,2)
      send_subarray_size(:,2,2) = [half_ihe-half_ihs+1,half_jhe-half_jhs+1,nlev(k)]
      recv_subarray_size(:,2,2) = send_subarray_size(:,2,2)

      recv_subarray_start(:,1,1) = [full_ihs,full_jhs,0]
      recv_subarray_start(:,2,1) = [half_ihs,full_jhs,0]
      recv_subarray_start(:,1,2) = [full_ihs,half_jhs,0]
      recv_subarray_start(:,2,2) = [half_ihs,half_jhs,0]
      select case (orient)
      case (west)
          send_subarray_start(:,1,1) = [full_ihe+1          ,full_jhs              ,0]
          send_subarray_start(:,2,1) = [half_ihe+1          ,full_jhs              ,0]
          send_subarray_start(:,1,2) = [full_ihe+1          ,half_jhs              ,0]
          send_subarray_start(:,2,2) = [half_ihe+1          ,half_jhs              ,0]
      case (east)
          send_subarray_start(:,1,1) = [full_ihs-this%lon_hw,full_jhs              ,0]
          send_subarray_start(:,2,1) = [half_ihs-this%lon_hw,full_jhs              ,0]
          send_subarray_start(:,1,2) = [full_ihs-this%lon_hw,half_jhs              ,0]
          send_subarray_start(:,2,2) = [half_ihs-this%lon_hw,half_jhs              ,0]
      case (south)
        if (mesh%has_south_pole()) then
          send_subarray_start(:,1,1) = [full_ihs            ,full_jhe+2            ,0] ! Skip pole grid.
          send_subarray_start(:,2,1) = [half_ihs            ,full_jhe+2            ,0] ! Skip pole grid.
          send_subarray_start(:,1,2) = [full_ihs            ,half_jhe+1            ,0]
          send_subarray_start(:,2,2) = [half_ihs            ,half_jhe+1            ,0]
        else
          send_subarray_start(:,1,1) = [full_ihs            ,full_jhe+1            ,0]
          send_subarray_start(:,2,1) = [half_ihs            ,full_jhe+1            ,0]
          send_subarray_start(:,1,2) = [full_ihs            ,half_jhe+1            ,0]
          send_subarray_start(:,2,2) = [half_ihs            ,half_jhe+1            ,0]
        end if
      case (north)
        if (mesh%has_north_pole()) then
          send_subarray_start(:,1,1) = [full_ihs            ,full_jhs-this%lat_hw-1,0] ! Skip pole grid.
          send_subarray_start(:,2,1) = [half_ihs            ,full_jhs-this%lat_hw-1,0] ! Skip pole grid.
          send_subarray_start(:,1,2) = [full_ihs            ,half_jhs-this%lat_hw  ,0]
          send_subarray_start(:,2,2) = [half_ihs            ,half_jhs-this%lat_hw  ,0]
        else
          send_subarray_start(:,1,1) = [full_ihs            ,full_jhs-this%lat_hw  ,0]
          send_subarray_start(:,2,1) = [half_ihs            ,full_jhs-this%lat_hw  ,0]
          send_subarray_start(:,1,2) = [full_ihs            ,half_jhs-this%lat_hw  ,0]
          send_subarray_start(:,2,2) = [half_ihs            ,half_jhs-this%lat_hw  ,0]
        end if
      case (south_west)
        if (mesh%has_south_pole()) then
          send_subarray_start(:,1,1) = [full_ihe+1          ,full_jhe+2            ,0] ! Skip pole grid.
          send_subarray_start(:,2,1) = [half_ihe+1          ,full_jhe+2            ,0] ! Skip pole grid.
          send_subarray_start(:,1,2) = [full_ihe+1          ,half_jhe+1            ,0]
          send_subarray_start(:,2,2) = [half_ihe+1          ,half_jhe+1            ,0]
        else
          send_subarray_start(:,1,1) = [full_ihe+1          ,full_jhe+1            ,0]
          send_subarray_start(:,2,1) = [half_ihe+1          ,full_jhe+1            ,0]
          send_subarray_start(:,1,2) = [full_ihe+1          ,half_jhe+1            ,0]
          send_subarray_start(:,2,2) = [half_ihe+1          ,half_jhe+1            ,0]
        end if
      case (south_east)
        if (mesh%has_south_pole()) then
          send_subarray_start(:,1,1) = [full_ihs-this%lon_hw,full_jhe+2            ,0] ! Skip pole grid.
          send_subarray_start(:,2,1) = [half_ihs-this%lon_hw,full_jhe+2            ,0] ! Skip pole grid.
          send_subarray_start(:,1,2) = [full_ihs-this%lon_hw,half_jhe+1            ,0]
          send_subarray_start(:,2,2) = [half_ihs-this%lon_hw,half_jhe+1            ,0]
        else
          send_subarray_start(:,1,1) = [full_ihs-this%lon_hw,full_jhe+1            ,0]
          send_subarray_start(:,2,1) = [half_ihs-this%lon_hw,full_jhe+1            ,0]
          send_subarray_start(:,1,2) = [full_ihs-this%lon_hw,half_jhe+1            ,0]
          send_subarray_start(:,2,2) = [half_ihs-this%lon_hw,half_jhe+1            ,0]
        end if
      case (north_west)
        if (mesh%has_north_pole()) then
          send_subarray_start(:,1,1) = [full_ihe+1          ,full_jhs-this%lat_hw-1,0] ! Skip pole grid.
          send_subarray_start(:,2,1) = [half_ihe+1          ,full_jhs-this%lat_hw-1,0] ! Skip pole grid.
          send_subarray_start(:,1,2) = [full_ihe+1          ,half_jhs-this%lat_hw  ,0]
          send_subarray_start(:,2,2) = [half_ihe+1          ,half_jhs-this%lat_hw  ,0]
        else
          send_subarray_start(:,1,1) = [full_ihe+1          ,full_jhs-this%lat_hw  ,0]
          send_subarray_start(:,2,1) = [half_ihe+1          ,full_jhs-this%lat_hw  ,0]
          send_subarray_start(:,1,2) = [full_ihe+1          ,half_jhs-this%lat_hw  ,0]
          send_subarray_start(:,2,2) = [half_ihe+1          ,half_jhs-this%lat_hw  ,0]
        end if
      case (north_east)
        if (mesh%has_north_pole()) then
          send_subarray_start(:,1,1) = [full_ihs-this%lon_hw,full_jhs-this%lat_hw-1,0] ! Skip pole grid.
          send_subarray_start(:,2,1) = [half_ihs-this%lon_hw,full_jhs-this%lat_hw-1,0] ! Skip pole grid.
          send_subarray_start(:,1,2) = [full_ihs-this%lon_hw,half_jhs-this%lat_hw  ,0]
          send_subarray_start(:,2,2) = [half_ihs-this%lon_hw,half_jhs-this%lat_hw  ,0]
        else
          send_subarray_start(:,1,1) = [full_ihs-this%lon_hw,full_jhs-this%lat_hw  ,0]
          send_subarray_start(:,2,1) = [half_ihs-this%lon_hw,full_jhs-this%lat_hw  ,0]
          send_subarray_start(:,1,2) = [full_ihs-this%lon_hw,half_jhs-this%lat_hw  ,0]
          send_subarray_start(:,2,2) = [half_ihs-this%lon_hw,half_jhs-this%lat_hw  ,0]
        end if
      end select
      do j = 1, 2
        do i = 1, 2
          call MPI_TYPE_CREATE_SUBARRAY(3, array_size(:,i,j), send_subarray_size(:,i,j), &
                                        send_subarray_start(:,i,j), MPI_ORDER_FORTRAN, dtype, &
                                        this%send_type_3d(i,j,k), ierr)
          call MPI_TYPE_COMMIT(this%send_type_3d(i,j,k), ierr)
          call MPI_TYPE_CREATE_SUBARRAY(3, array_size(:,i,j), recv_subarray_size(:,i,j), &
                                        recv_subarray_start(:,i,j), MPI_ORDER_FORTRAN, dtype, &
                                        this%recv_type_3d(i,j,k), ierr)
          call MPI_TYPE_COMMIT(this%recv_type_3d(i,j,k), ierr)
        end do
      end do
    end do

    do j = 1, 2
      do i = 1, 2
        call MPI_TYPE_CREATE_SUBARRAY(2, array_size(1:2,i,j), send_subarray_size(1:2,i,j), &
                                      send_subarray_start(1:2,i,j), MPI_ORDER_FORTRAN, dtype, &
                                      this%send_type_2d(i,j), ierr)
        call MPI_TYPE_COMMIT(this%send_type_2d(i,j), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(2, array_size(1:2,i,j), recv_subarray_size(1:2,i,j), &
                                      recv_subarray_start(1:2,i,j), MPI_ORDER_FORTRAN, dtype, &
                                      this%recv_type_2d(i,j), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_2d(i,j), ierr)
      end do
    end do

    this%initialized = .true.

  end subroutine latlon_halo_init

  subroutine latlon_halo_init_nest(this, parent_mesh, parent_proc_id)

    class(latlon_halo_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: parent_mesh
    integer, intent(in) :: parent_proc_id

  end subroutine latlon_halo_init_nest

  subroutine latlon_halo_clear(this)

    class(latlon_halo_type), intent(inout) :: this

    integer i, j, k
    integer ierr

    do k = 1, 2
      do j = 1, 2
        do i = 1, 2
          if (this%send_type_3d(i,j,k) /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%send_type_3d(i,j,k), ierr)
          if (this%recv_type_3d(i,j,k) /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%recv_type_3d(i,j,k), ierr)
        end do
      end do
    end do

    do j = 1, 2
      do i = 1, 2
        if (this%send_type_2d(i,j) /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%send_type_2d(i,j), ierr)
        if (this%recv_type_2d(i,j) /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%recv_type_2d(i,j), ierr)
      end do
    end do

    this%initialized = .false.

  end subroutine latlon_halo_clear

  subroutine latlon_halo_final(this)

    type(latlon_halo_type), intent(inout) :: this

    call this%clear()

  end subroutine latlon_halo_final

end module latlon_halo_mod
