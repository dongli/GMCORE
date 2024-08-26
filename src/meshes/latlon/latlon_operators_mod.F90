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
!   This module implements basic differential operators on the lat-lon meshes.
!
! History:
!
!   20240304: Initial creation.
!
! Authors:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
! ==============================================================================

module latlon_operators_mod

  use const_mod
  use latlon_mesh_mod
  use latlon_field_types_mod
  use latlon_parallel_mod

  implicit none

  private

  public divx_operator
  public divy_operator
  public div_operator
  public curl_operator

contains

  subroutine divx_operator(fx, divx)

    type(latlon_field3d_type), intent(in   ) :: fx
    type(latlon_field3d_type), intent(inout) :: divx

    integer ks, ke, i, j, k

    associate (mesh => divx%mesh)
    ks = merge(mesh%full_kds, mesh%half_kds, divx%loc(1:3) /= 'lev')
    ke = merge(mesh%full_kde, mesh%half_kde, divx%loc(1:3) /= 'lev')
    do k = ks, ke
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          divx%d(i,j,k) = (fx%d(i,j,k) - fx%d(i-1,j,k)) * mesh%le_lon(j) / mesh%area_cell(j)
        end do
      end do
    end do
    if (mesh%has_south_pole()) divx%d(:,mesh%full_jds,:) = 0
    if (mesh%has_north_pole()) divx%d(:,mesh%full_jde,:) = 0
    end associate

  end subroutine divx_operator

  subroutine divy_operator(fy, divy)

    type(latlon_field3d_type), intent(in   ) :: fy
    type(latlon_field3d_type), intent(inout) :: divy

    real(r8) work(divy%mesh%full_ids:divy%mesh%full_ide,divy%nlev)
    real(r8) pole(divy%nlev)
    integer ks, ke, i, j, k

    associate (mesh => divy%mesh)
    ks = merge(mesh%full_kds, mesh%half_kds, divy%loc(1:3) /= 'lev')
    ke = merge(mesh%full_kde, mesh%half_kde, divy%loc(1:3) /= 'lev')
    do k = ks, ke
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          divy%d(i,j,k) = (                    &
            fy%d(i,j  ,k) * mesh%le_lat(j  ) - &
            fy%d(i,j-1,k) * mesh%le_lat(j-1)   &
          ) / mesh%area_cell(j)
        end do
      end do
    end do
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = ks, ke
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = fy%d(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%area_pole_cap
      do k = ks, ke
        do i = mesh%full_ids, mesh%full_ide
          divy%d(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = ks, ke
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = fy%d(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%area_pole_cap
      do k = ks, ke
        do i = mesh%full_ids, mesh%full_ide
          divy%d(i,j,k) = -pole(k)
        end do
      end do
    end if
    end associate

  end subroutine divy_operator

  subroutine div_operator(fx, fy, div, with_halo)

    type(latlon_field3d_type), intent(in   ) :: fx
    type(latlon_field3d_type), intent(in   ) :: fy
    type(latlon_field3d_type), intent(inout) :: div
    logical, intent(in), optional :: with_halo

    real(r8) work(div%mesh%full_ids:div%mesh%full_ide,div%nlev)
    real(r8) pole(div%nlev)
    integer i, j, k, is, ie, js, je, ks, ke

    associate (mesh => div%mesh)
    ks = merge(mesh%full_kds, mesh%half_kds, div%loc(1:3) /= 'lev')
    ke = merge(mesh%full_kde, mesh%half_kde, div%loc(1:3) /= 'lev')
    is = mesh%full_ids
    ie = mesh%full_ide; if (present(with_halo)) ie = merge(ie + 1, ie, with_halo)
    js = mesh%full_jds_no_pole
    je = mesh%full_jde_no_pole; if (present(with_halo)) je = merge(je + 1, je, with_halo)
    do k = ks, ke
      do j = js, je
        do i = is, ie
          div%d(i,j,k) = ((                    &
            fx%d(i,j,k) - fx%d(i-1,j,k)        &
          ) * mesh%le_lon(j) + (               &
            fy%d(i,j  ,k) * mesh%le_lat(j  ) - &
            fy%d(i,j-1,k) * mesh%le_lat(j-1)   &
          )) / mesh%area_cell(j)
        end do
      end do
    end do
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = ks, ke
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = fy%d(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%area_pole_cap
      do k = ks, ke
        do i = mesh%full_ids, mesh%full_ide
          div%d(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = ks, ke
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = fy%d(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%area_pole_cap
      do k = ks, ke
        do i = mesh%full_ids, mesh%full_ide
          div%d(i,j,k) = -pole(k)
        end do
      end do
    end if
    end associate

  end subroutine div_operator

  subroutine curl_operator(fx, fy, curl, with_halo)

    type(latlon_field3d_type), intent(in   ) :: fx
    type(latlon_field3d_type), intent(in   ) :: fy
    type(latlon_field3d_type), intent(inout) :: curl
    logical, intent(in), optional :: with_halo

    integer i, j, k, is, ie, js, je, ks, ke

    associate (mesh => curl%mesh)
    ks = merge(mesh%full_kds, mesh%half_kds, curl%loc(1:3) /= 'lev')
    ke = merge(mesh%full_kde, mesh%half_kde, curl%loc(1:3) /= 'lev')
    is = mesh%half_ids
    ie = mesh%half_ide; if (present(with_halo)) ie = merge(ie + 1, ie, with_halo)
    js = mesh%half_jds
    je = mesh%half_jde; if (present(with_halo)) je = merge(je + 1, je, with_halo)
    do k = ks, ke
      do j = js, je
        do i = is, ie
          curl%d(i,j,k) = (                                                     &
            fx%d(i  ,j,k) * mesh%de_lon(j) - fx%d(i,j+1,k) * mesh%de_lon(j+1) + &
            fy%d(i+1,j,k) * mesh%de_lat(j) - fy%d(i,j  ,k) * mesh%de_lat(j  )   &
          ) / mesh%area_vtx(j)
        end do
      end do
    end do
    end associate

  end subroutine curl_operator

end module latlon_operators_mod
