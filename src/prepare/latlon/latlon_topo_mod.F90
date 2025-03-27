! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module latlon_topo_mod

  use const_mod
  use namelist_mod
  use topo_reader_mod
  use block_mod
  use latlon_parallel_mod
  use latlon_field_types_mod
  use math_mod
  use filter_mod

  implicit none

  private

  public latlon_topo_regrid
  public latlon_topo_smooth

contains

  subroutine fill_grid(lon1, lon2, lat1, lat2, gzs, std, lnd, cnt)

    real(r8), intent(in) :: lon1
    real(r8), intent(in) :: lon2
    real(r8), intent(in) :: lat1
    real(r8), intent(in) :: lat2
    real(r8), intent(inout) :: gzs
    real(r8), intent(inout) :: std
    real(r8), intent(inout) :: lnd
    integer , intent(inout) :: cnt

    integer is, ie, js, je, i, j

    do is = 1, size(topo_lon)
      if (lon1 <= topo_lon(is)) exit
    end do
    do ie = size(topo_lon), 1, -1
      if (lon2 >= topo_lon(ie)) exit
    end do
    do js = 1, size(topo_lat)
      if (lat1 <= topo_lat(js)) exit
    end do
    do je = size(topo_lat), 1, -1
      if (lat2 >= topo_lat(je)) exit
    end do

    select case (planet)
    case ('earth')
      do j = js, je
        do i = is, ie
          if (topo_lon(i) < lon1 .or. topo_lon(i) > lon2 .or. topo_lat(j) < lat1 .or. topo_lat(j) > lat2) then
            stop 999
          end if
          if (allocated(topo_mask)) then
            if (topo_mask(i,j) == 1) then
              gzs = gzs + topo_gzs(i,j)
              std = std + topo_gzs(i,j)**2
            end if
          else
            if (topo_gzs(i,j) > 0) then
              gzs = gzs + topo_gzs(i,j)
              std = std + topo_gzs(i,j)**2
            end if
          end if
        end do
      end do
      lnd = lnd + count(topo_gzs(is:ie,js:je) > 0)
    case ('mars')
      gzs = gzs + sum(topo_gzs(is:ie,js:je))
      std = std + sum(topo_gzs(is:ie,js:je)**2)
      lnd = lnd + (ie - is + 1) * (je - js + 1)
    end select
    cnt = cnt + (ie - is + 1) * (je - js + 1)

  end subroutine fill_grid

  subroutine latlon_topo_regrid(block)

    type(block_type), intent(inout) :: block

    real(r8) min_lon, max_lon, min_lat, max_lat
    real(r8) lon1, lon2, lat1, lat2, pole_gzs, pole_std, pole_lnd
    integer i, j, pole_n
    integer n(block%mesh%full_ids:block%mesh%full_ide)

    associate (mesh  => block%mesh           , &
               lnd   => block%static%landmask, &
               gzs   => block%static%gzs     , &
               std   => block%static%zs_std  , &
               dzsdx => block%static%dzsdx   , &
               dzsdy => block%static%dzsdy   )
    lnd%d = 0; gzs%d = 0; std%d = 0
    do j = mesh%full_jds, mesh%full_jde
      lat1 = mesh%half_lat_deg(j-1); lat1 = merge(lat1, -90.0_r8, lat1 /= inf)
      lat2 = mesh%half_lat_deg(j  ); lat2 = merge(lat2,  90.0_r8, lat2 /= inf)
      n = 0
      do i = mesh%full_ids, mesh%full_ide
        lon1 = mesh%half_lon_deg(i-1)
        lon2 = mesh%half_lon_deg(i  )
        call fill_grid(lon1, lon2, lat1, lat2, gzs%d(i,j), std%d(i,j), lnd%d(i,j), n(i))
        if (.not. mesh%is_pole(j)) then
          gzs%d(i,j) = gzs%d(i,j) / n(i)
          std%d(i,j) = (std%d(i,j) - 2 * gzs%d(i,j)**2 * n(i) + gzs%d(i,j)**2) / n(i) / g
          lnd%d(i,j) = lnd%d(i,j) / n(i)
        end if
      end do
      if (mesh%is_pole(j)) then
        call zonal_sum(proc%zonal_circle, gzs%d(mesh%full_ids:mesh%full_ide,j), pole_gzs)
        call zonal_sum(proc%zonal_circle, std%d(mesh%full_ids:mesh%full_ide,j), pole_std)
        call zonal_sum(proc%zonal_circle, lnd%d(mesh%full_ids:mesh%full_ide,j), pole_lnd)
        call zonal_sum(proc%zonal_circle, n, pole_n)
        gzs%d(mesh%full_ids:mesh%full_ide,j) = pole_gzs / pole_n
        std%d(mesh%full_ids:mesh%full_ide,j) = (pole_std - 2 * pole_gzs**2 * pole_n + pole_gzs**2) / pole_n / g
        lnd%d(mesh%full_ids:mesh%full_ide,j) = pole_lnd / pole_n
      end if
    end do
    call fill_halo(gzs)
    call fill_halo(std)
    call fill_halo(lnd)
    call calc_zs_slope(gzs, dzsdx, dzsdy)
    end associate

  end subroutine latlon_topo_regrid

  subroutine latlon_topo_smooth(block)

    type(block_type), intent(inout) :: block

    real(r8) tmp

    if (.not. use_topo_smooth) return
    if (proc%is_root()) call log_notice('Filter topography.')

    associate (lnd    => block%static%landmask, &
               gzs    => block%static%gzs     , &
               dzsdx  => block%static%dzsdx   , &
               dzsdy  => block%static%dzsdy   )
    tmp = global_max(proc%comm_model, maxval(gzs%d / g))
    if (proc%is_root()) call log_notice('Maximum zs is ' // to_str(tmp, 'F8.1') // ' m.')
    tmp = global_max(proc%comm_model, max(dzsdx%absmax(), dzsdy%absmax()))
    if (proc%is_root()) call log_notice('Maximum topography slope angle before filter is ' // to_str(atan(tmp) * deg, 'F7.2') // ' deg.')
    if (use_zs_polar_filter) call zs_polar_filter(block, gzs, dzsdx, dzsdy)
    call zs_grad_filter(block, lnd, gzs, dzsdx, dzsdy)
    tmp = global_max(proc%comm_model, maxval(gzs%d / g))
    if (proc%is_root()) call log_notice('Maximum zs is ' // to_str(tmp, 'F8.1') // ' m.')
    tmp = global_max(proc%comm_model, max(dzsdx%absmax(), dzsdy%absmax()))
    if (proc%is_root()) call log_notice('Maximum topography slope angle after filter is ' // to_str(atan(tmp) * deg, 'F7.2') // ' deg.')
    end associate

  end subroutine latlon_topo_smooth

  subroutine calc_zs_slope(gzs, dzsdx, dzsdy)

    type(latlon_field2d_type), intent(in) :: gzs
    type(latlon_field2d_type), intent(inout) :: dzsdx
    type(latlon_field2d_type), intent(inout) :: dzsdy

    integer i, j

    associate (mesh => gzs%mesh)
    do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
      do i = mesh%half_ids, mesh%half_ide
        dzsdx%d(i,j) = (gzs%d(i+1,j) - gzs%d(i,j)) / g / mesh%de_lon(j)
      end do
    end do
    do j = mesh%half_jds, mesh%half_jde
      do i = mesh%full_ids, mesh%full_ide
        dzsdy%d(i,j) = (gzs%d(i,j+1) - gzs%d(i,j)) / g / mesh%de_lat(j)
      end do
    end do
    call fill_halo(dzsdx)
    call fill_halo(dzsdy)
    end associate

  end subroutine calc_zs_slope

  subroutine zs_polar_filter(block, gzs, dzsdx, dzsdy)

    type(block_type), intent(in) :: block
    type(latlon_field2d_type), intent(inout) :: gzs
    type(latlon_field2d_type), intent(inout) :: dzsdx
    type(latlon_field2d_type), intent(inout) :: dzsdy

    type(latlon_field2d_type) gzs_f
    real(r8) lat0, wgt
    integer j, cyc

    associate (mesh => block%filter_mesh, halo => block%filter_halo)
    call gzs_f%init('', '', '', 'cell', mesh, halo)
    lat0 = abs(global_mesh%full_lat_deg(2))
    do cyc = 1, 3
      call filter_run(block%big_filter, gzs, gzs_f)
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        wgt = exp_two_values(1.0_r8, 0.0_r8, lat0, 60.0_r8, abs(mesh%full_lat_deg(j)))
        gzs%d(:,j) = wgt * gzs_f%d(:,j) + (1 - wgt) * gzs%d(:,j)
      end do
      call fill_halo(gzs)
    end do
    call calc_zs_slope(gzs, dzsdx, dzsdy)
    call gzs_f%clear()
    end associate

  end subroutine zs_polar_filter

  subroutine zs_grad_filter(block, lnd, gzs, dzsdx, dzsdy)

    type(block_type), intent(inout) :: block
    type(latlon_field2d_type), intent(in) :: lnd
    type(latlon_field2d_type), intent(inout) :: gzs
    type(latlon_field2d_type), intent(inout) :: dzsdx
    type(latlon_field2d_type), intent(inout) :: dzsdy

    type(latlon_field2d_type) g1, g2, fx, fy
    real(r8) c0, c1, max_slope
    real(r8) work(gzs%mesh%full_ids:gzs%mesh%full_ide), pole
    integer i, j, k, cyc

    ! Do Laplacian damping on target grid with terrain slope larger than a threshold.

    associate (mesh  => gzs%mesh          , &
               halo  => gzs%halo          , &
               order => topo_smooth_order , &
               dfx   => block%dtend%dmgsdt, &
               f     => gzs               )
    c0 = (-1)**(order / 2) * topo_smooth_coef
    call g1%init('g1', '', '', 'cell', mesh, halo)
    call g2%init('g2', '', '', 'cell', mesh, halo)
    call fx%init('fx', '', '', 'lon' , mesh, halo)
    call fy%init('fy', '', '', 'lat' , mesh, halo)
    do cyc = 1, topo_smooth_cycles
      call g1%copy(f, with_halo=.true.)
      do k = 1, (order - 2) / 2
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids - 1, mesh%half_ide
            fx%d(i,j) = (g1%d(i+1,j) - g1%d(i,j)) / mesh%de_lon(j)
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            fy%d(i,j) = (g1%d(i,j+1) - g1%d(i,j)) / mesh%de_lat(j)
          end do
        end do
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide
            g2%d(i,j) = (                                  &
              (fx%d(i,j) - fx%d(i-1,j)) * mesh%le_lon(j) + &
               fy%d(i,j  ) * mesh%le_lat(j  ) -            &
               fy%d(i,j-1) * mesh%le_lat(j-1)              &
            ) / mesh%area_cell(j)
          end do
        end do
        if (mesh%has_south_pole()) then
          j = mesh%full_jds
          do i = mesh%full_ids, mesh%full_ide
            work(i) = fy%d(i,j)
          end do
          call zonal_sum(proc%zonal_circle, work, pole)
          pole = pole * mesh%le_lat(j) / global_mesh%area_pole_cap
          do i = mesh%full_ids, mesh%full_ide
            g2%d(i,j) = pole
          end do
        end if
        if (mesh%has_north_pole()) then
          j = mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            work(i) = fy%d(i,j-1)
          end do
          call zonal_sum(proc%zonal_circle, work, pole)
          pole = -pole * mesh%le_lat(j-1) / global_mesh%area_pole_cap
          do i = mesh%full_ids, mesh%full_ide
            g2%d(i,j) = pole
          end do
        end if
        call fill_halo(g2)
        g1%d = g2%d
      end do
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        c1 = exp_two_values(1.0e2_r8, 1.0_r8, 90.0_r8, 45.0_r8, abs(mesh%full_lat_deg(j)))
        do i = mesh%half_ids - 1, mesh%half_ide
          fx%d(i,j) = (g1%d(i+1,j) - g1%d(i,j)) / mesh%de_lon(j)
          fx%d(i,j) = c1 * abs(abs(dzsdx%d(i,j)) - topo_max_slope) / topo_max_slope * fx%d(i,j)
          fx%d(i,j) = max(0.0_r8, min(lnd%d(i,j), lnd%d(i+1,j))) * fx%d(i,j)
        end do
      end do
      do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          fy%d(i,j) = (g1%d(i,j+1) - g1%d(i,j)) / mesh%de_lat(j)
          fy%d(i,j) = abs(abs(dzsdy%d(i,j)) - topo_max_slope) / topo_max_slope * fy%d(i,j)
          fy%d(i,j) = max(0.0_r8, min(lnd%d(i,j), lnd%d(i,j+1))) * fy%d(i,j)
        end do
      end do
      if (order > 2) then
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids - 1, mesh%half_ide
            fx%d(i,j) = fx%d(i,j) * max(0.0_r8, sign(1.0_r8, -fx%d(i,j) * (f%d(i+1,j) - f%d(i,j))))
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            fy%d(i,j) = fy%d(i,j) * max(0.0_r8, sign(1.0_r8, -fy%d(i,j) * (f%d(i,j+1) - f%d(i,j))))
          end do
        end do
      end if
      ! Filter the zonal tendency.
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          dfx%d(i,j) = (fx%d(i,j) - fx%d(i-1,j)) * mesh%le_lon(j) / mesh%area_cell(j)
        end do
      end do
      call filter_run(block%small_filter, dfx)
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          f%d(i,j) = f%d(i,j) - c0 * (dfx%d(i,j) + ( &
              fy%d(i,j  ) * mesh%le_lat(j  ) -       &
              fy%d(i,j-1) * mesh%le_lat(j-1)         &
            ) / mesh%area_cell(j))
        end do
      end do
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do i = mesh%full_ids, mesh%full_ide
          work(i) = fy%d(i,j)
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = c0 * pole * mesh%le_lat(j) / global_mesh%area_pole_cap
        do i = mesh%full_ids, mesh%full_ide
          f%d(i,j) = f%d(i,j) - pole
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          work(i) = fy%d(i,j-1)
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = -c0 * pole * mesh%le_lat(j-1) / global_mesh%area_pole_cap
        do i = mesh%full_ids, mesh%full_ide
          f%d(i,j) = f%d(i,j) - pole
        end do
      end if
      call fill_halo(f)
      call calc_zs_slope(f, dzsdx, dzsdy)
    end do
    call g1%clear()
    call g2%clear()
    call fx%clear()
    call fy%clear()
    end associate

  end subroutine zs_grad_filter

end module latlon_topo_mod
