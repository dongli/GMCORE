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
!   This module implements FFSL (Flux-Form Semi-Lagrangian) advection scheme on
!   the lat-lon grid.
!
! Authors:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
! ==============================================================================

module ffsl_mod

  use const_mod
  use namelist_mod
  use latlon_mesh_mod, only: global_mesh
  use latlon_field_types_mod
  use latlon_parallel_mod
  use process_mod, only: proc
  use adv_batch_mod
  use ppm_mod
  use limiter_mod
  use perf_mod

  implicit none

  private

  public ffsl_init
  public ffsl_calc_mass_hflx
  public ffsl_calc_tracer_hflx
  public ffsl_calc_tracer_vflx

  interface
    subroutine hflx_interface(batch, u, v, mx, my, mfx, mfy, dt)
      import  adv_batch_type, latlon_field3d_type, r8
      type(adv_batch_type     ), intent(inout) :: batch
      type(latlon_field3d_type), intent(in   ) :: u
      type(latlon_field3d_type), intent(in   ) :: v
      type(latlon_field3d_type), intent(in   ) :: mx
      type(latlon_field3d_type), intent(in   ) :: my
      type(latlon_field3d_type), intent(inout) :: mfx
      type(latlon_field3d_type), intent(inout) :: mfy
      real(r8), intent(in) :: dt
    end subroutine hflx_interface
    subroutine vflx_interface(batch, w, m, mfz, dt)
      import adv_batch_type, latlon_field3d_type, r8
      type(adv_batch_type     ), intent(inout) :: batch
      type(latlon_field3d_type), intent(in   ) :: w
      type(latlon_field3d_type), intent(in   ) :: m
      type(latlon_field3d_type), intent(inout) :: mfz
      real(r8), intent(in) :: dt
    end subroutine vflx_interface
  end interface

  procedure(hflx_interface), pointer :: hflx => null()
  procedure(vflx_interface), pointer :: vflx => null()

contains

  subroutine ffsl_init()

    select case (ffsl_flux_type)
    case ('van_leer')
      hflx => hflx_van_leer
      vflx => vflx_van_leer
    case ('ppm')
      hflx => hflx_ppm
      vflx => vflx_ppm
    case default
      if (proc%is_root()) call log_error('Invalid ffsl_flux_type ' // trim(ffsl_flux_type) // '!')
    end select

    call limiter_init()

  end subroutine ffsl_init

  subroutine ffsl_calc_mass_hflx(batch, m, mfx, mfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: m
    type(latlon_field3d_type), intent(inout) :: mfx
    type(latlon_field3d_type), intent(inout) :: mfy
    real(r8), intent(in), optional :: dt

    integer i, j, k
    real(r8) work(m%mesh%full_ids:m%mesh%full_ide,m%mesh%full_nlev)
    real(r8) pole(m%mesh%full_nlev)
    real(r8) dt_opt

    call perf_start('ffsl_calc_mass_hflx')

    dt_opt = batch%dt; if (present(dt)) dt_opt = dt

    associate (mesh => m%mesh    , &
               u    => batch%u   , & ! in
               v    => batch%v   , & ! in
               divx => batch%divx, & ! in
               divy => batch%divy, & ! in
               mx   => batch%qx  , & ! work array
               my   => batch%qy  )   ! work array
    ! Run inner advective operators.
    call hflx(batch, u, v, m, m, mfx, mfy, dt_opt)
    call fill_halo(mfx, south_halo=.false., north_halo=.false., east_halo=.false.)
    call fill_halo(mfy, west_halo=.false., east_halo=.false., north_halo=.false.)
    ! Calculate intermediate tracer density due to advective operators.
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          ! Subtract divergence terms from flux to form advective operators.
          mx%d(i,j,k) = m%d(i,j,k) - 0.5_r8 * (      &
            (                                        &
              mfx%d(i,j,k) - mfx%d(i-1,j,k)          &
            ) * mesh%le_lon(j) / mesh%area_cell(j) - &
            divx%d(i,j,k) * m%d(i,j,k)               &
          ) * dt_opt
          my%d(i,j,k) = m%d(i,j,k) - 0.5_r8 * (     &
            (                                       &
              mfy%d(i,j  ,k) * mesh%le_lat(j  ) -   &
              mfy%d(i,j-1,k) * mesh%le_lat(j-1)     &
            ) / mesh%area_cell(j) -                 &
            divy%d(i,j,k) * m%d(i,j,k)              &
          ) * dt_opt
        end do
      end do
    end do
    ! Handle the Pole boundary conditions.
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = mfy%d(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%area_pole_cap
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          mx%d(i,j,k) = m%d(i,j,k)
          my%d(i,j,k) = m%d(i,j,k) - 0.5_r8 * (pole(k) - divy%d(i,j,k) * m%d(i,j,k)) * dt_opt
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = mfy%d(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%area_pole_cap
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          mx%d(i,j,k) = m%d(i,j,k)
          my%d(i,j,k) = m%d(i,j,k) + 0.5_r8 * (pole(k) - divy%d(i,j,k) * m%d(i,j,k)) * dt_opt
        end do
      end do
    end if
    call fill_halo(mx, west_halo=.false., east_halo=.false.)
    call fill_halo(my, south_halo=.false., north_halo=.false.)
    ! Run outer flux form operators.
    call hflx(batch, u, v, my, mx, mfx, mfy, dt_opt)
    end associate

    call perf_stop('ffsl_calc_mass_hflx')

  end subroutine ffsl_calc_mass_hflx

  subroutine ffsl_calc_tracer_hflx(batch, q, qmfx, qmfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: q
    type(latlon_field3d_type), intent(inout) :: qmfx
    type(latlon_field3d_type), intent(inout) :: qmfy
    real(r8), intent(in), optional :: dt

    integer ks, ke, i, j, k
    real(r8) work(q%mesh%full_ids:q%mesh%full_ide,q%mesh%half_nlev)
    real(r8) pole(q%mesh%half_nlev)
    real(r8) dt_opt

    call perf_start('ffsl_calc_tracer_hflx')

    dt_opt = batch%dt; if (present(dt)) dt_opt = dt

    associate (mesh  => q%mesh     , &
               m     => batch%m    , & ! in
               u     => batch%u    , & ! in
               v     => batch%v    , & ! in
               mfx   => batch%mfx  , & ! in
               mfy   => batch%mfy  , & ! in
               divx  => batch%divx , & ! in
               divy  => batch%divy , & ! in
               qx    => batch%qx   , & ! work array
               qy    => batch%qy   , & ! work array
               qmfx0 => batch%qmfx0, & ! out
               qmfy0 => batch%qmfy0)   ! out
    ! Run inner advective operators.
    call hflx(batch, mfx, mfy, q, q, qmfx0, qmfy0, dt_opt)
    call fill_halo(qmfx0, south_halo=.false., north_halo=.false., east_halo=.false.)
    call fill_halo(qmfy0, west_halo=.false., east_halo=.false., north_halo=.false.)
    select case (batch%loc)
    case ('cell', 'lev')
      ks = merge(mesh%full_kds, mesh%half_kds, batch%loc == 'cell')
      ke = merge(mesh%full_kde, mesh%half_kde, batch%loc == 'cell')
      ! Calculate intermediate tracer density due to advective operators.
      do k = ks, ke
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide
            qx%d(i,j,k) = (q%d(i,j,k) - dt_opt * (     &
              (                                        &
                qmfx0%d(i,j,k) - qmfx0%d(i-1,j,k)      &
              ) * mesh%le_lon(j) / mesh%area_cell(j)   &
            ) / m%d(i,j,k)) / (1 - dt_opt * divx%d(i,j,k))
            qy%d(i,j,k) = (q%d(i,j,k) - dt_opt * (     &
              (                                        &
                qmfy0%d(i,j  ,k) * mesh%le_lat(j  ) -  &
                qmfy0%d(i,j-1,k) * mesh%le_lat(j-1)    &
              ) / mesh%area_cell(j)                    &
            ) / m%d(i,j,k)) / (1 - dt_opt * divy%d(i,j,k))
          end do
        end do
      end do
      ! Handle the Pole boundary conditions.
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = qmfy0%d(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work(:,ks:ke), pole(ks:ke))
        pole(ks:ke) = dt_opt * pole(ks:ke) * mesh%le_lat(j) / global_mesh%area_pole_cap
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            qx%d(i,j,k) = q%d(i,j,k)
            qy%d(i,j,k) = (q%d(i,j,k) - pole(k) / m%d(i,j,k)) / (1 - dt_opt * divy%d(i,j,k))
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = qmfy0%d(i,j-1,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work(:,ks:ke), pole(ks:ke))
        pole(ks:ke) = dt_opt * pole(ks:ke) * mesh%le_lat(j-1) / global_mesh%area_pole_cap
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            qx%d(i,j,k) = q%d(i,j,k)
            qy%d(i,j,k) = (q%d(i,j,k) + pole(k) / m%d(i,j,k)) / (1 - dt_opt * divy%d(i,j,k))
          end do
        end do
      end if
    case ('vtx')
    end select
    call fill_halo(qx, west_halo=.false., east_halo=.false.)
    call fill_halo(qy, south_halo=.false., north_halo=.false.)
    ! Run outer flux form operators.
    call hflx(batch, mfx, mfy, qy, qx, qmfx, qmfy, dt_opt)
    ! Do SWIFT splitting.
    do k = ks, ke
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          qmfx%d(i,j,k) = 0.5_r8 * (qmfx0%d(i,j,k) + qmfx%d(i,j,k))
        end do
      end do
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          qmfy%d(i,j,k) = 0.5_r8 * (qmfy0%d(i,j,k) + qmfy%d(i,j,k))
        end do
      end do
    end do
    end associate

    call perf_stop('ffsl_calc_tracer_hflx')

  end subroutine ffsl_calc_tracer_hflx

  subroutine ffsl_calc_tracer_vflx(batch, q, qmfz, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: q
    type(latlon_field3d_type), intent(inout) :: qmfz
    real(r8), intent(in), optional :: dt

    real(r8) dt_opt

    call perf_start('ffsl_calc_tracer_vflx')

    dt_opt = batch%dt; if (present(dt)) dt_opt = dt

    call vflx(batch, batch%we, q, qmfz, dt_opt)

    call perf_stop('ffsl_calc_tracer_vflx')

  end subroutine ffsl_calc_tracer_vflx

  subroutine hflx_van_leer(batch, u, v, qx, qy, qmfx, qmfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: u
    type(latlon_field3d_type), intent(in   ) :: v
    type(latlon_field3d_type), intent(in   ) :: qx
    type(latlon_field3d_type), intent(in   ) :: qy
    type(latlon_field3d_type), intent(inout) :: qmfx
    type(latlon_field3d_type), intent(inout) :: qmfy
    real(r8), intent(in) :: dt

    integer ks, ke, i, j, k, iu, ju, ci
    real(r8) cf, dq

    associate (mesh => u%mesh    , &
               m    => batch%m   , & ! in
               cflx => batch%cflx, & ! in
               cfly => batch%cfly)   ! in
    select case (batch%loc)
    case ('cell', 'lev')
      ks = merge(mesh%full_kds, mesh%half_kds, batch%loc == 'cell')
      ke = merge(mesh%full_kde, mesh%half_kde, batch%loc == 'cell')
      do k = ks, ke
        ! Along x-axis
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            ci = int(cflx%d(i,j,k))
            cf = cflx%d(i,j,k) - ci
            if (abs(cflx%d(i,j,k)) < 1.0e-16_r8) then
              qmfx%d(i,j,k) = 0
            else if (cflx%d(i,j,k) > 0) then
              iu = i - ci
              dq = slope(qx%d(iu-1,j,k), qx%d(iu,j,k), qx%d(iu+1,j,k))
              qmfx%d(i,j,k) = (cf * (qx%d(iu,j,k) + dq * 0.5_r8 * (1 - cf))) * m%d(iu,j,k) * mesh%de_lon(j) / dt &
                            + sum(m%d(iu+1:i,j,k) * qx%d(iu+1:i,j,k)) * mesh%de_lon(j) / dt
            else
              iu = i - ci + 1
              dq = slope(qx%d(iu-1,j,k), qx%d(iu,j,k), qx%d(iu+1,j,k))
              qmfx%d(i,j,k) = (cf * (qx%d(iu,j,k) - dq * 0.5_r8 * (1 + cf))) * m%d(iu,j,k) * mesh%de_lon(j) / dt &
                            - sum(m%d(i+1:iu-1,j,k) * qx%d(i+1:iu-1,j,k)) * mesh%de_lon(j) / dt
            end if
          end do
        end do
        ! Along y-axis
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            cf = cfly%d(i,j,k)
            ju = merge(j, j + 1, cf > 0)
            dq = slope(qy%d(i,ju-1,k), qy%d(i,ju,k), qy%d(i,ju+1,k))
            qmfy%d(i,j,k) = v%d(i,j,k) * (qy%d(i,ju,k) + dq * 0.5_r8 * (sign(1.0_r8, cf) - cf))
          end do
        end do
      end do
    case ('vtx')
    end select
    end associate

  end subroutine hflx_van_leer

  subroutine vflx_van_leer(batch, w, q, qmfz, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: w
    type(latlon_field3d_type), intent(in   ) :: q
    type(latlon_field3d_type), intent(inout) :: qmfz
    real(r8), intent(in) :: dt

    integer i, j, k, ku, ci
    real(r8) cf, dq

    associate (mesh => q%mesh    , &
               m    => batch%m   , & ! in
               cflz => batch%cflz)   ! in
    select case (batch%loc)
    case ('cell')
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            ci = int(cflz%d(i,j,k))
            cf = cflz%d(i,j,k) - ci
            if (abs(cflz%d(i,j,k)) < 1.0e-16_r8) then
              qmfz%d(i,j,k) = 0
            else if (cflz%d(i,j,k) > 0) then
              ku = k - ci - 1
              dq = slope(q%d(i,j,ku-1), q%d(i,j,ku), q%d(i,j,ku+1))
              qmfz%d(i,j,k) = (cf * (q%d(i,j,ku) + dq * 0.5_r8 * (1 - cf))) * w%d(i,j,k) / cflz%d(i,j,k) &
                            + sum(m%d(i,j,ku+1:k-1) * q%d(i,j,ku+1:k-1)) / dt 
            else
              ku = k - ci
              dq = slope(q%d(i,j,ku-1), q%d(i,j,ku), q%d(i,j,ku+1))
              qmfz%d(i,j,k) = (cf * (q%d(i,j,ku) - dq * 0.5_r8 * (1 + cf))) * w%d(i,j,k) / cflz%d(i,j,k) &
                            - sum(m%d(i,j,k:ku-1) * q%d(i,j,k:ku-1)) / dt
            end if
          end do
        end do
      end do
    case ('lev')
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            ci = int(cflz%d(i,j,k))
            cf = cflz%d(i,j,k) - ci
            if (abs(cflz%d(i,j,k)) < 1.0e-16_r8) then
              qmfz%d(i,j,k) = 0
            else if (cflz%d(i,j,k) > 0) then
              ku = k - ci
              dq = slope(q%d(i,j,ku-1), q%d(i,j,ku), q%d(i,j,ku+1))
              qmfz%d(i,j,k) = (cf * (q%d(i,j,ku) + dq * 0.5_r8 * (1 - cf))) * w%d(i,j,k) / cflz%d(i,j,k) &
                            + sum(m%d(i,j,ku+1:k) * q%d(i,j,ku+1:k)) / dt
            else
              ku = k - ci + 1
              dq = slope(q%d(i,j,ku-1), q%d(i,j,ku), q%d(i,j,ku+1))
              qmfz%d(i,j,k) = (cf * (q%d(i,j,ku) - dq * 0.5_r8 * (1 + cf))) * w%d(i,j,k) / cflz%d(i,j,k) &
                            - sum(m%d(i,j,k+1:ku-1) * q%d(i,j,k+1:ku-1)) / dt
            end if
          end do
        end do
      end do
    end select
    end associate

  end subroutine vflx_van_leer

  subroutine hflx_ppm(batch, u, v, qx, qy, qmfx, qmfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: u
    type(latlon_field3d_type), intent(in   ) :: v
    type(latlon_field3d_type), intent(in   ) :: qx
    type(latlon_field3d_type), intent(in   ) :: qy
    type(latlon_field3d_type), intent(inout) :: qmfx
    type(latlon_field3d_type), intent(inout) :: qmfy
    real(r8), intent(in) :: dt

    integer ks, ke, i, j, k, iu, ju, ci
    real(r8) cf, s1, s2, ds1, ds2, ds3, ql, dq, q6

    associate (mesh => u%mesh    , &
               m    => batch%m   , & ! in
               cflx => batch%cflx, & ! in
               cfly => batch%cfly)   ! in
    select case (batch%loc)
    case ('cell', 'lev')
      ks = merge(mesh%full_kds, mesh%half_kds, batch%loc == 'cell')
      ke = merge(mesh%full_kde, mesh%half_kde, batch%loc == 'cell')
      do k = ks, ke
        ! Along x-axis
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            ci = int(cflx%d(i,j,k))
            cf = cflx%d(i,j,k) - ci
            if (abs(cflx%d(i,j,k)) < 1.0e-16_r8) then
              qmfx%d(i,j,k) = 0
            else if (cflx%d(i,j,k) > 0) then
              iu = i - ci
              call ppm(qx%d(iu-2,j,k), qx%d(iu-1,j,k), qx%d(iu,j,k), qx%d(iu+1,j,k), qx%d(iu+2,j,k), ql, dq, q6)
              s1 = 1 - cf
              s2 = 1
              ds1 = s2    - s1
              ds2 = s2**2 - s1**2
              ds3 = s2**3 - s1**3
              qmfx%d(i,j,k) = ((ql * ds1 + 0.5_r8 * dq * ds2 + q6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) * m%d(iu,j,k) &
                            + sum(m%d(iu+1:i,j,k) * qx%d(iu+1:i,j,k))) * mesh%de_lon(j) / dt
            else
              iu = i - ci + 1
              call ppm(qx%d(iu-2,j,k), qx%d(iu-1,j,k), qx%d(iu,j,k), qx%d(iu+1,j,k), qx%d(iu+2,j,k), ql, dq, q6)
              s1 = 0
              s2 = -cf
              ds1 = s2    - s1
              ds2 = s2**2 - s1**2
              ds3 = s2**3 - s1**3
              qmfx%d(i,j,k) = -((ql * ds1 + 0.5_r8 * dq * ds2 + q6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) * m%d(iu,j,k) &
                            + sum(m%d(i+1:iu-1,j,k) * qx%d(i+1:iu-1,j,k))) * mesh%de_lon(j) / dt
            end if
          end do
        end do
        ! Along y-axis
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            if (abs(cfly%d(i,j,k)) < 1.0e-16_r8) then
              qmfy%d(i,j,k) = 0
            else if (cfly%d(i,j,k) > 0) then
              ju = j
              call ppm(qy%d(i,ju-2,k), qy%d(i,ju-1,k), qy%d(i,ju,k), qy%d(i,ju+1,k), qy%d(i,ju+2,k), ql, dq, q6)
              s1 = 1 - cfly%d(i,j,k)
              s2 = 1
              ds1 = s2    - s1
              ds2 = s2**2 - s1**2
              ds3 = s2**3 - s1**3
              qmfy%d(i,j,k) =  v%d(i,j,k) * (ql * ds1 + 0.5_r8 * dq * ds2 + q6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cfly%d(i,j,k)
            else if (cfly%d(i,j,k) < 0) then
              ju = j + 1
              call ppm(qy%d(i,ju-2,k), qy%d(i,ju-1,k), qy%d(i,ju,k), qy%d(i,ju+1,k), qy%d(i,ju+2,k), ql, dq, q6)
              s1 = 0
              s2 = -cfly%d(i,j,k)
              ds1 = s2    - s1
              ds2 = s2**2 - s1**2
              ds3 = s2**3 - s1**3
              qmfy%d(i,j,k) = -v%d(i,j,k) * (ql * ds1 + 0.5_r8 * dq * ds2 + q6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cfly%d(i,j,k)
            end if
          end do
        end do
      end do
    case ('vtx')
    end select
    end associate

  end subroutine hflx_ppm

  subroutine vflx_ppm(batch, w, q, qmfz, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: w
    type(latlon_field3d_type), intent(in   ) :: q
    type(latlon_field3d_type), intent(inout) :: qmfz
    real(r8), intent(in) :: dt

    integer i, j, k, ku, ci
    real(r8) cf, s1, s2, ds1, ds2, ds3, ql, dq, q6

    associate (mesh => q%mesh    , &
               m    => batch%m   , & ! in
               cflz => batch%cflz)   ! in
    select case (batch%loc)
    case ('cell')
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            ci = int(cflz%d(i,j,k))
            cf = cflz%d(i,j,k) - ci
            if (abs(cflz%d(i,j,k)) < 1.0e-16_r8) then
              qmfz%d(i,j,k) = 0
            else if (cflz%d(i,j,k) > 0) then
              ku = k - ci - 1
              ku = min(max(ku, mesh%full_kms + 2), mesh%full_kme - 2)
              call ppm(q%d(i,j,ku-2), q%d(i,j,ku-1), q%d(i,j,ku), q%d(i,j,ku+1), q%d(i,j,ku+2), ql, dq, q6)
              s1 = 1 - cf
              s2 = 1
              ds1 = s2    - s1
              ds2 = s2**2 - s1**2
              ds3 = s2**3 - s1**3
              qmfz%d(i,j,k) = (ql * ds1 + 0.5_r8 * dq * ds2 + q6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) * w%d(i,j,k) / cflz%d(i,j,k) &
                            + sum(m%d(i,j,ku+1:k-1) * q%d(i,j,ku+1:k-1)) / dt
                              
            else
              ku = k - ci
              ku = min(max(ku, mesh%full_kms + 2), mesh%full_kme - 2)
              call ppm(q%d(i,j,ku-2), q%d(i,j,ku-1), q%d(i,j,ku), q%d(i,j,ku+1), q%d(i,j,ku+2), ql, dq, q6)
              s1 = 0
              s2 = -cf
              ds1 = s2    - s1
              ds2 = s2**2 - s1**2
              ds3 = s2**3 - s1**3
              qmfz%d(i,j,k) = -(ql * ds1 + 0.5_r8 * dq * ds2 + q6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) * w%d(i,j,k) / cflz%d(i,j,k) &
                            - sum(m%d(i,j,k:ku-1) * q%d(i,j,k:ku-1)) / dt
                              
            end if
          end do
        end do
      end do
    case ('lev')
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            ci = int(cflz%d(i,j,k))
            cf = cflz%d(i,j,k) - ci
            if (abs(cflz%d(i,j,k)) < 1.0e-16_r8) then
              qmfz%d(i,j,k) = 0
            else if (cflz%d(i,j,k) > 0) then
              ku = k - ci
              call ppm(q%d(i,j,ku-2), q%d(i,j,ku-1), q%d(i,j,ku), q%d(i,j,ku+1), q%d(i,j,ku+2), ql, dq, q6)
              s1 = 1 - cf
              s2 = 1
              ds1 = s2    - s1
              ds2 = s2**2 - s1**2
              ds3 = s2**3 - s1**3
              qmfz%d(i,j,k) = (ql * ds1 + 0.5_r8 * dq * ds2 + q6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) * w%d(i,j,k) / cflz%d(i,j,k) &
                            + sum(m%d(i,j,ku+1:k) * q%d(i,j,ku+1:k)) / dt
            else
              ku = k - ci + 1
              call ppm(q%d(i,j,ku-2), q%d(i,j,ku-1), q%d(i,j,ku), q%d(i,j,ku+1), q%d(i,j,ku+2), ql, dq, q6)
              s1 = 0
              s2 = -cf
              ds1 = s2    - s1
              ds2 = s2**2 - s1**2
              ds3 = s2**3 - s1**3
              qmfz%d(i,j,k) = -(ql * ds1 + 0.5_r8 * dq * ds2 + q6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) * w%d(i,j,k) / cflz%d(i,j,k) &
                            - sum(m%d(i,j,k+1:ku-1) * q%d(i,j,k+1:ku-1)) / dt
            end if
          end do
        end do
      end do
    end select
    end associate

  end subroutine vflx_ppm

end module ffsl_mod
