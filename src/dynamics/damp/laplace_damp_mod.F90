! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module laplace_damp_mod

  ! ∂q       n+1    2n
  ! -- = (-1)    K ∇ q
  ! ∂t

  use const_mod
  use namelist_mod
  use latlon_mesh_mod, only: global_mesh
  use latlon_field_types_mod
  use latlon_parallel_mod
  use process_mod, only: proc
  use block_mod
  use filter_mod

  implicit none

  private

  public laplace_damp_init
  public laplace_damp_final
  public laplace_damp_run
  public laplace_damp

  interface laplace_damp
    module procedure laplace_damp_2d
    module procedure laplace_damp_3d
  end interface laplace_damp

contains

  subroutine laplace_damp_init()

  end subroutine laplace_damp_init

  subroutine laplace_damp_final()

  end subroutine laplace_damp_final

  subroutine laplace_damp_run(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    call laplace_damp(block, dstate%u_lon, laplace_damp_order, laplace_damp_coef, block%aux%dudt_damp)
    call fill_halo(dstate%u_lon, async=.true.)
    call laplace_damp(block, dstate%v_lat, laplace_damp_order, laplace_damp_coef, block%aux%dvdt_damp)
    call fill_halo(dstate%v_lat, async=.true.)
    if (nonhydrostatic) then
      call laplace_damp(block, dstate%w_lev, laplace_damp_order, laplace_damp_coef, block%aux%dwdt_damp)
      call fill_halo(dstate%w_lev, async=.true.)
    end if

  end subroutine laplace_damp_run

  subroutine laplace_damp_2d(block, f, order, coef)

    type(block_type), intent(inout) :: block
    type(latlon_field2d_type), intent(inout) :: f
    integer, intent(in) :: order
    real(r8), intent(in) :: coef

    real(r8) work(f%mesh%full_ids:f%mesh%full_ide), pole
    real(r8) c0
    integer i, j, l

    call wait_halo(f)

    c0 = (-1)**(order / 2) * coef

    select case (f%loc)
    case ('cell')
      associate (mesh => block%mesh       , &
                 g1   => block%aux%g1_2d  , &
                 g2   => block%aux%g2_2d  , &
                 fx   => block%aux%fx_2d  , &
                 fy   => block%aux%fy_2d  , &
                 dxdt => block%aux%dxdt_2d)
      call g1%copy(f, with_halo=.true.)
      do l = 1, (order - 2) / 2
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
        do i = mesh%half_ids - 1, mesh%half_ide
          fx%d(i,j) = (g1%d(i+1,j) - g1%d(i,j)) / mesh%de_lon(j)
        end do
      end do
      do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          fy%d(i,j) = (g1%d(i,j+1) - g1%d(i,j)) / mesh%de_lat(j)
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
          dxdt%d(i,j) = (fx%d(i,j) - fx%d(i-1,j)) * mesh%le_lon(j) / mesh%area_cell(j)
        end do
      end do
      call filter_run(block%big_filter, dxdt)
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          f%d(i,j) = f%d(i,j) - c0 * (dxdt%d(i,j) + ( &
              fy%d(i,j  ) * mesh%le_lat(j  ) - &
              fy%d(i,j-1) * mesh%le_lat(j-1)   &
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
      end associate
    end select

    call fill_halo(f)

  end subroutine laplace_damp_2d

  subroutine laplace_damp_3d(block, f, order, coef, dfx)

    type(block_type), intent(inout) :: block
    type(latlon_field3d_type), intent(inout) :: f
    integer, intent(in) :: order
    real(r8), intent(in) :: coef
    type(latlon_field3d_type), intent(inout) :: dfx

    real(r8) work(f%mesh%full_ids:f%mesh%full_ide,f%mesh%full_nlev), pole(f%mesh%full_nlev)
    real(r8) c0
    integer i, j, k, l

    call wait_halo(f)

    c0 = (-1)**(order / 2) * coef

    select case (f%loc)
    case ('cell')
      associate (mesh => block%mesh     , &
                 g1   => block%aux%g1_3d, &
                 g2   => block%aux%g2_3d, &
                 fx   => block%aux%fx_3d, &
                 fy   => block%aux%fy_3d)
      call g1%copy(f, with_halo=.true.)
      do l = 1, (order - 2) / 2
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids - 1, mesh%half_ide
              fx%d(i,j,k) = (g1%d(i+1,j,k) - g1%d(i,j,k)) / mesh%de_lon(j)
            end do
          end do
          do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              fy%d(i,j,k) = (g1%d(i,j+1,k) - g1%d(i,j,k)) / mesh%de_lat(j)
            end do
          end do
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%full_ids, mesh%full_ide
              g2%d(i,j,k) = (                                      &
                (fx%d(i,j  ,k) - fx%d(i-1,j,k)) * mesh%le_lon(j) + &
                 fy%d(i,j  ,k) * mesh%le_lat(j  ) -                &
                 fy%d(i,j-1,k) * mesh%le_lat(j-1)                  &
              ) / mesh%area_cell(j)
            end do
          end do
        end do
        if (mesh%has_south_pole()) then
          j = mesh%full_jds
          do k = mesh%full_kds, mesh%full_kde
            do i = mesh%full_ids, mesh%full_ide
              work(i,k) = fy%d(i,j,k)
            end do
          end do
          call zonal_sum(proc%zonal_circle, work, pole)
          pole = pole * mesh%le_lat(j) / global_mesh%area_pole_cap
          do k = mesh%full_kds, mesh%full_kde
            do i = mesh%full_ids, mesh%full_ide
              g2%d(i,j,k) = pole(k)
            end do
          end do
        end if
        if (mesh%has_north_pole()) then
          j = mesh%full_jde
          do k = mesh%full_kds, mesh%full_kde
            do i = mesh%full_ids, mesh%full_ide
              work(i,k) = fy%d(i,j-1,k)
            end do
          end do
          call zonal_sum(proc%zonal_circle, work, pole)
          pole = -pole * mesh%le_lat(j-1) / global_mesh%area_pole_cap
          do k = mesh%full_kds, mesh%full_kde
            do i = mesh%full_ids, mesh%full_ide
              g2%d(i,j,k) = pole(k)
            end do
          end do
        end if
        call fill_halo(g2)
        g1%d = g2%d
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids - 1, mesh%half_ide
            fx%d(i,j,k) = (g1%d(i+1,j,k) - g1%d(i,j,k)) / mesh%de_lon(j)
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            fy%d(i,j,k) = (g1%d(i,j+1,k) - g1%d(i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do
      if (order > 2) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids - 1, mesh%half_ide
              fx%d(i,j,k) = fx%d(i,j,k) * max(0.0_r8, sign(1.0_r8, -fx%d(i,j,k) * (f%d(i+1,j,k) - f%d(i,j,k))))
            end do
          end do
          do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              fy%d(i,j,k) = fy%d(i,j,k) * max(0.0_r8, sign(1.0_r8, -fy%d(i,j,k) * (f%d(i,j+1,k) - f%d(i,j,k))))
            end do
          end do
        end do
      end if
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide
            dfx%d(i,j,k) = (fx%d(i,j,k) - fx%d(i-1,j,k)) * mesh%le_lon(j) / mesh%area_cell(j)
          end do
        end do
      end do
      call filter_run(block%big_filter, dfx)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide
            f%d(i,j,k) = f%d(i,j,k) - c0 * (dfx%d(i,j,k) + ( &
                fy%d(i,j  ,k) * mesh%le_lat(j  ) -           &
                fy%d(i,j-1,k) * mesh%le_lat(j-1)             &
              ) / mesh%area_cell(j))
          end do
        end do
      end do
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = fy%d(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = c0 * pole * mesh%le_lat(j) / global_mesh%area_pole_cap
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            f%d(i,j,k) = f%d(i,j,k) - pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = fy%d(i,j-1,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = -c0 * pole * mesh%le_lat(j-1) / global_mesh%area_pole_cap
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            f%d(i,j,k) = f%d(i,j,k) - pole(k)
          end do
        end do
      end if
      end associate
    case ('lon')
      associate (mesh => block%mesh         , &
                 g1   => block%aux%g1_3d_lon, &
                 g2   => block%aux%g2_3d_lon, &
                 fx   => block%aux%fx_3d_lon, &
                 fy   => block%aux%fy_3d_lon)
      call g1%copy(f, with_halo=.true.)
      do l = 1, (order - 2) / 2
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%full_ids, mesh%full_ide + 1
              fx%d(i,j,k) = (g1%d(i,j,k) - g1%d(i-1,j,k)) / mesh%de_lon(j)
            end do
          end do
          do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
            do i = mesh%half_ids, mesh%half_ide
              fy%d(i,j,k) = (g1%d(i,j+1,k) - g1%d(i,j,k)) / mesh%de_lat(j)
            end do
          end do
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              g2%d(i,j,k) = (                                    &
                (fx%d(i+1,j,k) - fx%d(i,j,k)) * mesh%le_lon(j) + &
                 fy%d(i,j  ,k) * mesh%le_lat(j  ) -              &
                 fy%d(i,j-1,k) * mesh%le_lat(j-1)                &
              ) / (2 * mesh%area_lon(j))
            end do
          end do
        end do
        call fill_halo(g2)
        g1%d = g2%d
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide + 1
            fx%d(i,j,k) = (g1%d(i,j,k) - g1%d(i-1,j,k)) / mesh%de_lon(j)
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%half_ids, mesh%half_ide
            fy%d(i,j,k) = (g1%d(i,j+1,k) - g1%d(i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do
      if (order > 2) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%full_ids, mesh%full_ide + 1
              fx%d(i,j,k) = fx%d(i,j,k) * max(0.0_r8, sign(1.0_r8, -fx%d(i,j,k) * (f%d(i,j,k) - f%d(i-1,j,k))))
            end do
          end do
          do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
            do i = mesh%half_ids, mesh%half_ide
              fy%d(i,j,k) = fy%d(i,j,k) * max(0.0_r8, sign(1.0_r8, -fy%d(i,j,k) * (f%d(i,j+1,k) - f%d(i,j,k))))
            end do
          end do
        end do
      end if
      ! Filter the zonal tendency.
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            dfx%d(i,j,k) = (fx%d(i+1,j,k) - fx%d(i,j,k)) * mesh%le_lon(j) / (2 * mesh%area_lon(j))
          end do
        end do
      end do
      call filter_run(block%big_filter, dfx)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            f%d(i,j,k) = f%d(i,j,k) - c0 * (dfx%d(i,j,k) + ( &
              fy%d(i,j  ,k) * mesh%le_lat(j  ) -             &
              fy%d(i,j-1,k) * mesh%le_lat(j-1)               &
            ) / (2 * mesh%area_lon(j)))
          end do
        end do
      end do
      end associate
    case ('lat')
      associate (mesh => block%mesh         , &
                 g1   => block%aux%g1_3d_lat, &
                 g2   => block%aux%g2_3d_lat, &
                 fx   => block%aux%fx_3d_lat, &
                 fy   => block%aux%fy_3d_lat)
      call g1%copy(f, with_halo=.true.)
      do l = 1, (order - 2) / 2
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%half_ids - 1, mesh%half_ide
              fx%d(i,j,k) = (g1%d(i+1,j,k) - g1%d(i,j,k)) / mesh%le_lat(j)
            end do
          end do
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
            do i = mesh%full_ids, mesh%full_ide
              fy%d(i,j,k) = (g1%d(i,j,k) - g1%d(i,j-1,k)) / mesh%le_lon(j)
            end do
          end do
          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              g2%d(i,j,k) = (                                    &
                (fx%d(i,j,k) - fx%d(i-1,j,k)) * mesh%de_lat(j) + &
                 fy%d(i,j+1,k) * mesh%de_lon(j+1) -              &
                 fy%d(i,j  ,k) * mesh%de_lon(j  )                &
              ) / (2 * mesh%area_lat(j))
            end do
          end do
        end do
        call fill_halo(g2)
        g1%d = g2%d
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%half_ids - 1, mesh%half_ide
            fx%d(i,j,k) = (g1%d(i+1,j,k) - g1%d(i,j,k)) / mesh%le_lat(j)
          end do
        end do
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_ids, mesh%full_ide
            fy%d(i,j,k) = (g1%d(i,j,k) - g1%d(i,j-1,k)) / mesh%le_lon(j)
          end do
        end do
      end do
      if (order > 2) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%half_ids - 1, mesh%half_ide
              fx%d(i,j,k) = fx%d(i,j,k) * max(0.0_r8, sign(1.0_r8, -fx%d(i,j,k) * (f%d(i+1,j,k) - f%d(i,j,k))))
            end do
          end do
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
            do i = mesh%full_ids, mesh%full_ide
              fy%d(i,j,k) = fy%d(i,j,k) * max(0.0_r8, sign(1.0_r8, -fy%d(i,j,k) * (f%d(i,j,k) - f%d(i,j-1,k))))
            end do
          end do
        end do
      end if
      ! Filter the zonal tendency.
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            dfx%d(i,j,k) = (fx%d(i,j,k) - fx%d(i-1,j,k)) * mesh%de_lat(j) / (2 * mesh%area_lat(j))
          end do
        end do
      end do
      call filter_run(block%big_filter, dfx)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            f%d(i,j,k) = f%d(i,j,k) - c0 * (dfx%d(i,j,k) + ( &
              fy%d(i,j+1,k) * mesh%de_lon(j+1) -             &
              fy%d(i,j  ,k) * mesh%de_lon(j  )               &
            ) / (2 * mesh%area_lat(j)))
          end do
        end do
      end do
      end associate
    case ('lev')
      associate (mesh => block%mesh         , &
                 g1   => block%aux%g1_3d_lev, &
                 g2   => block%aux%g2_3d_lev, &
                 fx   => block%aux%fx_3d_lev, &
                 fy   => block%aux%fy_3d_lev)
      call g1%copy(f, with_halo=.true.)
      do l = 1, (order - 2) / 2
        do k = mesh%half_kds + 1, mesh%half_kde - 1
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids - 1, mesh%half_ide
              fx%d(i,j,k) = (g1%d(i+1,j,k) - g1%d(i,j,k)) / mesh%de_lon(j)
            end do
          end do
          do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              fy%d(i,j,k) = (g1%d(i,j+1,k) - g1%d(i,j,k)) / mesh%de_lat(j)
            end do
          end do
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%full_ids, mesh%full_ide
              g2%d(i,j,k) = (                                    &
                (fx%d(i,j,k) - fx%d(i-1,j,k)) * mesh%le_lon(j) + &
                 fy%d(i,j  ,k) * mesh%le_lat(j  ) -              &
                 fy%d(i,j-1,k) * mesh%le_lat(j-1)                &
              ) / mesh%area_cell(j)
            end do
          end do
        end do
        if (mesh%has_south_pole()) then
          j = mesh%full_jds
          do k = mesh%half_kds + 1, mesh%half_kde - 1
            do i = mesh%full_ids, mesh%full_ide
              work(i,k) = fy%d(i,j,k)
            end do
          end do
          call zonal_sum(proc%zonal_circle, work, pole)
          pole = pole * mesh%le_lat(j) / global_mesh%area_pole_cap
          do k = mesh%half_kds + 1, mesh%half_kde - 1
            do i = mesh%full_ids, mesh%full_ide
              g2%d(i,j,k) = pole(k)
            end do
          end do
        end if
        if (mesh%has_north_pole()) then
          j = mesh%full_jde
          do k = mesh%half_kds + 1, mesh%half_kde - 1
            do i = mesh%full_ids, mesh%full_ide
              work(i,k) = fy%d(i,j-1,k)
            end do
          end do
          call zonal_sum(proc%zonal_circle, work, pole)
          pole = -pole * mesh%le_lat(j-1) / global_mesh%area_pole_cap
          do k = mesh%half_kds + 1, mesh%half_kde - 1
            do i = mesh%full_ids, mesh%full_ide
              g2%d(i,j,k) = pole(k)
            end do
          end do
        end if
        call fill_halo(g2)
        g1%d = g2%d
      end do
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids - 1, mesh%half_ide
            fx%d(i,j,k) = (g1%d(i+1,j,k) - g1%d(i,j,k)) / mesh%de_lon(j)
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            fy%d(i,j,k) = (g1%d(i,j+1,k) - g1%d(i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do
      if (order > 2) then
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids - 1, mesh%half_ide
            fx%d(i,j,k) = fx%d(i,j,k) * max(0.0_r8, sign(1.0_r8, -fx%d(i,j,k) * (f%d(i+1,j,k) - f%d(i,j,k))))
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            fy%d(i,j,k) = fy%d(i,j,k) * max(0.0_r8, sign(1.0_r8, -fy%d(i,j,k) * (f%d(i,j+1,k) - f%d(i,j,k))))
          end do
        end do
      end if
      ! Filter the zonal tendency.
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide
            dfx%d(i,j,k) = (fx%d(i,j,k) - fx%d(i-1,j,k)) * mesh%le_lon(j) / mesh%area_cell(j)
          end do
        end do
      end do
      call filter_run(block%big_filter, dfx)
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide
            f%d(i,j,k) = f%d(i,j,k) - c0 * (dfx%d(i,j,k) + ( &
              fy%d(i,j  ,k) * mesh%le_lat(j  ) -             &
              fy%d(i,j-1,k) * mesh%le_lat(j-1)               &
            ) / mesh%area_cell(j))
          end do
        end do
      end do
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do k = mesh%half_kds + 1, mesh%half_kde - 1
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = fy%d(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = c0 * pole * mesh%le_lat(j) / global_mesh%area_pole_cap
        do k = mesh%half_kds + 1, mesh%half_kde - 1
          do i = mesh%full_ids, mesh%full_ide
            f%d(i,j,k) = f%d(i,j,k) - pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do k = mesh%half_kds + 1, mesh%half_kde - 1
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = fy%d(i,j-1,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = -c0 * pole * mesh%le_lat(j-1) / global_mesh%area_pole_cap
        do k = mesh%half_kds + 1, mesh%half_kde - 1
          do i = mesh%full_ids, mesh%full_ide
            f%d(i,j,k) = f%d(i,j,k) - pole(k)
          end do
        end do
      end if
      end associate
    end select

  end subroutine laplace_damp_3d

end module laplace_damp_mod
