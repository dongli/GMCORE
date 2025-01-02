! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module interp_mod

  use const_mod
  use namelist_mod
  use latlon_mesh_mod, mesh_type => latlon_mesh_type
  use latlon_field_types_mod

  implicit none

  private

  !                              / lev_lat
  !               o-------------o------------o lev_vtx
  !              /|            /            /|
  !             / |                        / |
  !            /  |        |              /  |
  !           o   |        o lev        -o- lev_lon
  !          /    |        |    /       /    |
  !         /     o vtx        o lat   /     o vtx
  !        /      |           /       /      |
  !       o-------+-----o------------o       |
  !       |       |                  |       |
  !      lon -o-  |        o cell    |  -o- lon
  !       |       |                  |       |
  !       |       o------------------+-------o
  !       |      /       /           |      /
  !       o vtx /       o lat        o vtx /
  !       |    /       /   |         |    /
  !       |   o            o lev     |   o
  !       |  /             |         |  /
  !       | /                        | /
  !       |/                         |/
  !       o-------------o------------o
  !

  public interp_init
  public interp_final
  public interp_run
  public average_run

  interface interp_run
    module procedure interp_run_3d
  end interface interp_run

  interface average_run
    module procedure average_run_3d
  end interface average_run

contains

  subroutine interp_init()

    call interp_final()

  end subroutine interp_init

  subroutine interp_final()

  end subroutine interp_final

  subroutine interp_run_3d(x, y, extrap)

    type(latlon_field3d_type), intent(in   ) :: x
    type(latlon_field3d_type), intent(inout) :: y
    logical, intent(in), optional :: extrap

    logical extrap_opt
    real(r8) a, b, c, x1, x2, x3
    integer i, j, k

    extrap_opt = .true.; if (present(extrap)) extrap_opt = extrap

    select case (trim(x%loc) // '>' // trim(y%loc))
    ! --------------------------------------------------------------------------
    case ('cell>lon')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%full_jds_no_pole, x%mesh%full_jde_no_pole
          do i = x%mesh%half_ids, x%mesh%half_ide
            y%d(i,j,k) = (x%mesh%area_lon_west(j) * x%d(i  ,j,k) + &
                          x%mesh%area_lon_east(j) * x%d(i+1,j,k)   &
                         ) / x%mesh%area_lon(j)
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('cell>lat')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%half_jds, x%mesh%half_jde
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = (x%mesh%area_lat_north(j) * x%d(i,j+1,k) + &
                          x%mesh%area_lat_south(j) * x%d(i,j  ,k)   &
                         ) / x%mesh%area_lat(j)
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('cell>lev')
      if (x%mesh%full_nlev == 1) return
      ! -------
      !
      ! ===o=== k-1
      !
      ! ---?--- k
      !
      ! ===o=== k
      !
      ! -------
      do k = x%mesh%half_kds + 1, x%mesh%half_kde - 1
        a = x%mesh%full_dlev(k-1) / (2 * x%mesh%half_dlev(k))
        b = x%mesh%full_dlev(k  ) / (2 * x%mesh%half_dlev(k))
        do j = x%mesh%full_jds, x%mesh%full_jde
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = a * x%d(i,j,k-1) + b * x%d(i,j,k)
          end do
        end do
      end do
      if (extrap_opt) then
        k = x%mesh%half_kds
        x1 = x%mesh%full_lev(k  ) - x%mesh%half_lev(k)
        x2 = x%mesh%full_lev(k+1) - x%mesh%half_lev(k)
        a =  x2 / (x2 - x1)
        b = -x1 / (x2 - x1)
        do j = x%mesh%full_jds, x%mesh%full_jde
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = a * x%d(i,j,k) + b * x%d(i,j,k+1)
          end do
        end do
        k = x%mesh%half_kde
        x1 = x%mesh%half_lev(k) - x%mesh%full_lev(k-1)
        x2 = x%mesh%half_lev(k) - x%mesh%full_lev(k-2)
        a =  x2 / (x2 - x1)
        b = -x1 / (x2 - x1)
        do j = x%mesh%full_jds, x%mesh%full_jde
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = a * x%d(i,j,k-1) + b * x%d(i,j,k-2)
          end do
        end do
      else
        k = x%mesh%half_kds
        y%d(:,:,k) = x%d(:,:,k)
        k = x%mesh%half_kde
        y%d(:,:,k) = x%d(:,:,k-1)
      end if
    ! --------------------------------------------------------------------------
    case ('cell>vtx')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%half_jds, x%mesh%half_jde
          do i = x%mesh%half_ids, x%mesh%half_ide
            y%d(i,j,k) = (                                                   &
              (x%d(i,j  ,k) + x%d(i+1,j  ,k)) * x%mesh%area_subcell(2,j  ) + &
              (x%d(i,j+1,k) + x%d(i+1,j+1,k)) * x%mesh%area_subcell(1,j+1)   &
            ) / x%mesh%area_vtx(j)
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('lon>cell')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%full_jds_no_pole, x%mesh%full_jde_no_pole
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = (x%mesh%area_lon_east(j) * x%d(i-1,j,k) + &
                          x%mesh%area_lon_west(j) * x%d(i  ,j,k)   &
                         ) / x%mesh%area_lon(j)
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('lon>lat')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%half_jds, x%mesh%half_jde
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = x%mesh%tg_wgt_lat(1,j) * (x%d(i-1,j  ,k) + x%d(i,j  ,k)) + &
                         x%mesh%tg_wgt_lat(2,j) * (x%d(i-1,j+1,k) + x%d(i,j+1,k))
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('lat>cell')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%full_jds_no_pole, x%mesh%full_jde_no_pole
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = (x%mesh%area_lat_south(j  ) * x%d(i,j  ,k) + &
                          x%mesh%area_lat_north(j-1) * x%d(i,j-1,k)   &
                         ) / (x%mesh%area_lat_south(j) + x%mesh%area_lat_north(j-1))
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('lat>lon')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%full_jds_no_pole, x%mesh%full_jde_no_pole
          do i = x%mesh%half_ids, x%mesh%half_ide
            y%d(i,j,k) = x%mesh%tg_wgt_lon(1,j) * (x%d(i,j-1,k) + x%d(i+1,j-1,k)) + &
                         x%mesh%tg_wgt_lon(2,j) * (x%d(i,j  ,k) + x%d(i+1,j  ,k))
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('lev>cell')
      ! =======
      !
      ! ---o--- k
      !
      ! ===?=== k
      !
      ! ---o--- k+1
      !
      ! =======
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%full_jds, x%mesh%full_jde
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = 0.5_r8 * (x%d(i,j,k) + x%d(i,j,k+1))
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('lon>lev_lon')
      ! -------
      !
      ! ===o=== k-1
      !
      ! ---?--- k
      !
      ! ===o=== k
      !
      ! ----o--
      do k = x%mesh%half_kds + 1, x%mesh%half_kde - 1
        a = x%mesh%full_dlev(k-1) / (x%mesh%full_dlev(k-1) + x%mesh%full_dlev(k))
        b = x%mesh%full_dlev(k  ) / (x%mesh%full_dlev(k-1) + x%mesh%full_dlev(k))
        do j = x%mesh%full_jds, x%mesh%full_jde
          do i = x%mesh%half_ids, x%mesh%half_ide
            y%d(i,j,k) = a * x%d(i,j,k) + b * x%d(i,j,k-1)
          end do
        end do
      end do
      k = x%mesh%half_kds
      ! ---?--- 1
      !
      ! ===o=== 1   x1
      !
      ! -------
      !
      ! ===o=== 2   x2
      !
      ! -------
      !
      ! ===o=== 3   x3
      x1 = x%mesh%full_lev(k  ) - x%mesh%half_lev(k)
      x2 = x%mesh%full_lev(k+1) - x%mesh%half_lev(k)
      x3 = x%mesh%full_lev(k+2) - x%mesh%half_lev(k)
      a =  x2 * x3 / (x1**2 - x1 * x2 - x1 * x3 + x2 * x3)
      b =  x1 * x3 / (x2**2 - x2 * x1 - x2 * x3 + x1 * x3)
      c =  x1 * x2 / (x3**2 - x3 * x1 - x3 * x2 + x1 * x2)
      do j = x%mesh%full_jds, x%mesh%full_jde
        do i = x%mesh%half_ids, x%mesh%half_ide
          y%d(i,j,k) = a * x%d(i,j,k) + b * x%d(i,j,k+1) + c * x%d(i,j,k+2)
        end do
      end do
      k = x%mesh%half_kde
      ! ===o=== NLEV - 2  x3
      !
      ! -------
      !
      ! ===o=== NLEV - 1  x2
      !
      ! -------
      !
      ! ===o=== NLEV      x1
      !
      ! ---?--- NLEV + 1
      x1 = x%mesh%half_lev(k) - x%mesh%full_lev(k-1)
      x2 = x%mesh%half_lev(k) - x%mesh%full_lev(k-2)
      x3 = x%mesh%half_lev(k) - x%mesh%full_lev(k-3)
      a =  x2 * x3 / (x1**2 - x1 * x2 - x1 * x3 + x2 * x3)
      b =  x1 * x3 / (x2**2 - x2 * x1 - x2 * x3 + x1 * x3)
      c =  x1 * x2 / (x3**2 - x3 * x1 - x3 * x2 + x1 * x2)
      do j = x%mesh%full_jds, x%mesh%full_jde
        do i = x%mesh%half_ids, x%mesh%half_ide
          y%d(i,j,k) = a * x%d(i,j,k-1) + b * x%d(i,j,k-2) + c * x%d(i,j,k-3)
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('lat>lev_lat')
      ! -------
      !
      ! ===o=== k-1
      !
      ! ---?--- k
      !
      ! ===o=== k
      !
      ! -------
      do k = x%mesh%half_kds + 1, x%mesh%half_kde - 1
        a = x%mesh%full_dlev(k-1) / (x%mesh%full_dlev(k-1) + x%mesh%full_dlev(k))
        b = x%mesh%full_dlev(k  ) / (x%mesh%full_dlev(k-1) + x%mesh%full_dlev(k))
        do j = x%mesh%half_jds, x%mesh%half_jde
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = a * x%d(i,j,k) + b * x%d(i,j,k-1)
          end do
        end do
      end do
      k = x%mesh%half_kds
      ! ---?--- 1
      !
      ! ===o=== 1   x1
      !
      ! -------
      !
      ! ===o=== 2   x2
      !
      ! -------
      !
      ! ===o=== 3   x3
      x1 = x%mesh%full_lev(k  ) - x%mesh%half_lev(k)
      x2 = x%mesh%full_lev(k+1) - x%mesh%half_lev(k)
      x3 = x%mesh%full_lev(k+2) - x%mesh%half_lev(k)
      a =  x2 * x3 / (x1**2 - x1 * x2 - x1 * x3 + x2 * x3)
      b =  x1 * x3 / (x2**2 - x2 * x1 - x2 * x3 + x1 * x3)
      c =  x1 * x2 / (x3**2 - x3 * x1 - x3 * x2 + x1 * x2)
      do j = x%mesh%half_jds, x%mesh%half_jde
        do i = x%mesh%full_ids, x%mesh%full_ide
          y%d(i,j,k) = a * x%d(i,j,k) + b * x%d(i,j,k+1) + c * x%d(i,j,k+2)
        end do
      end do
      k = x%mesh%half_kde
      ! ===o=== NLEV - 2  x3
      !
      ! -------
      !
      ! ===o=== NLEV - 1  x2
      !
      ! -------
      !
      ! ===o=== NLEV      x1
      !
      ! ---?--- NLEV + 1
      x1 = x%mesh%half_lev(k) - x%mesh%full_lev(k-1)
      x2 = x%mesh%half_lev(k) - x%mesh%full_lev(k-2)
      x3 = x%mesh%half_lev(k) - x%mesh%full_lev(k-3)
      a =  x2 * x3 / (x1**2 - x1 * x2 - x1 * x3 + x2 * x3)
      b =  x1 * x3 / (x2**2 - x2 * x1 - x2 * x3 + x1 * x3)
      c =  x1 * x2 / (x3**2 - x3 * x1 - x3 * x2 + x1 * x2)
      do j = x%mesh%half_jds, x%mesh%half_jde
        do i = x%mesh%full_ids, x%mesh%full_ide
          y%d(i,j,k) = a * x%d(i,j,k-1) + b * x%d(i,j,k-2) + c * x%d(i,j,k-3)
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('lev>lev_lon')
      do k = x%mesh%half_kds, x%mesh%half_kde
        do j = x%mesh%full_jds_no_pole, x%mesh%full_jde_no_pole
          do i = x%mesh%half_ids, x%mesh%half_ide
            y%d(i,j,k) = (x%mesh%area_lon_west(j) * x%d(i  ,j,k) + &
                          x%mesh%area_lon_east(j) * x%d(i+1,j,k)   &
                         ) / x%mesh%area_lon(j)
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('lev>lev_lat')
      do k = x%mesh%half_kds, x%mesh%half_kde
        do j = x%mesh%half_jds, x%mesh%half_jde
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = (x%mesh%area_lat_north(j) * x%d(i,j+1,k) + &
                          x%mesh%area_lat_south(j) * x%d(i,j  ,k)   &
                         ) / x%mesh%area_lat(j)
          end do
        end do
      end do
    end select

  end subroutine interp_run_3d

  subroutine average_run_3d(x, y)

    type(latlon_field3d_type), intent(in   ) :: x
    type(latlon_field3d_type), intent(inout) :: y

    integer i, j, k


    select case (trim(x%loc) // '>' // trim(y%loc))
    ! --------------------------------------------------------------------------
    case ('cell>lon')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%full_jds_no_pole, x%mesh%full_jde_no_pole
          do i = x%mesh%half_ids, x%mesh%half_ide
            y%d(i,j,k) = (x%d(i,j,k) + x%d(i+1,j,k)) * 0.5_r8
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('cell>lat')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%half_jds, x%mesh%half_jde
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = (x%d(i,j,k) + x%d(i,j+1,k)) * 0.5_r8
          end do
        end do
      end do
    case ('lev>cell')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%full_jds, x%mesh%full_jde
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = (x%d(i,j,k) + x%d(i,j,k+1)) * 0.5_r8
          end do
        end do
      end do
    end select

  end subroutine average_run_3d

end module interp_mod
