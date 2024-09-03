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
!   In August 2024, SWIFT splitting replaced the orginal Lin and Rood (1996)
!   splitting.
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
  use latlon_operators_mod
  use latlon_parallel_mod
  use process_mod, only: proc
  use adv_batch_mod
  use ppm_mod
  use limiter_mod
  use perf_mod

  implicit none

  private

  public ffsl_init
  public ffsl_calc_mass_hflx_lin96
  public ffsl_calc_mass_hflx_swift
  public ffsl_calc_mass_vflx
  public ffsl_calc_tracer_hflx_lin96
  public ffsl_calc_tracer_hflx_swift
  public ffsl_calc_tracer_vflx

  interface
    subroutine hflx_mass_interface(batch, u, u_frac, v, mx, my, mfx, mfy, dt)
      import  adv_batch_type, latlon_field3d_type, r8
      type(adv_batch_type     ), intent(inout) :: batch
      type(latlon_field3d_type), intent(in   ) :: u
      type(latlon_field3d_type), intent(in   ) :: u_frac
      type(latlon_field3d_type), intent(in   ) :: v
      type(latlon_field3d_type), intent(in   ) :: mx
      type(latlon_field3d_type), intent(in   ) :: my
      type(latlon_field3d_type), intent(inout) :: mfx
      type(latlon_field3d_type), intent(inout) :: mfy
      real(r8), intent(in) :: dt
    end subroutine hflx_mass_interface
    subroutine hflx_tracer_interface(batch, mx, my, cflx, cfly, mfx, mfx_frac, mfy, qx, qy, qmfx, qmfy, dt)
      import  adv_batch_type, latlon_field3d_type, r8
      type(adv_batch_type     ), intent(inout) :: batch
      type(latlon_field3d_type), intent(in   ) :: mx
      type(latlon_field3d_type), intent(in   ) :: my
      type(latlon_field3d_type), intent(in   ) :: cflx
      type(latlon_field3d_type), intent(in   ) :: cfly
      type(latlon_field3d_type), intent(in   ) :: mfx
      type(latlon_field3d_type), intent(in   ) :: mfx_frac
      type(latlon_field3d_type), intent(in   ) :: mfy
      type(latlon_field3d_type), intent(in   ) :: qx
      type(latlon_field3d_type), intent(in   ) :: qy
      type(latlon_field3d_type), intent(inout) :: qmfx
      type(latlon_field3d_type), intent(inout) :: qmfy
      real(r8), intent(in) :: dt
    end subroutine hflx_tracer_interface
    subroutine vflx_mass_interface(batch, w, w_frac, m, mfz, dt)
      import adv_batch_type, latlon_field3d_type, r8
      type(adv_batch_type     ), intent(inout) :: batch
      type(latlon_field3d_type), intent(in   ) :: w
      type(latlon_field3d_type), intent(in   ) :: w_frac
      type(latlon_field3d_type), intent(in   ) :: m
      type(latlon_field3d_type), intent(inout) :: mfz
      real(r8), intent(in) :: dt
    end subroutine vflx_mass_interface
    subroutine vflx_tracer_interface(batch, w, w_frac, m, mfz, dt)
      import adv_batch_type, latlon_field3d_type, r8
      type(adv_batch_type     ), intent(inout) :: batch
      type(latlon_field3d_type), intent(in   ) :: w
      type(latlon_field3d_type), intent(in   ) :: w_frac
      type(latlon_field3d_type), intent(in   ) :: m
      type(latlon_field3d_type), intent(inout) :: mfz
      real(r8), intent(in) :: dt
    end subroutine vflx_tracer_interface
  end interface

  procedure(hflx_mass_interface), pointer :: hflx_mass => null()
  procedure(vflx_mass_interface), pointer :: vflx_mass => null()
  procedure(hflx_tracer_interface), pointer :: hflx_tracer => null()
  procedure(vflx_tracer_interface), pointer :: vflx_tracer => null()

contains

  subroutine ffsl_init()

    select case (ffsl_flux_type)
    case ('ppm')
      hflx_mass => hflx_ppm_mass
      vflx_mass => vflx_ppm_mass
      hflx_tracer => hflx_ppm_tracer
      vflx_tracer => vflx_ppm_tracer
    case default
      if (proc%is_root()) call log_error('Invalid ffsl_flux_type ' // trim(ffsl_flux_type) // '!')
    end select

    call limiter_init()

  end subroutine ffsl_init

  subroutine ffsl_calc_mass_hflx_lin96(batch, m, mfx, mfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: m
    type(latlon_field3d_type), intent(inout) :: mfx
    type(latlon_field3d_type), intent(inout) :: mfy
    real(r8), intent(in), optional :: dt

    integer i, j, k
    real(r8) dt_opt, half_dt

    call perf_start('ffsl_calc_mass_hflx_lin96')

    dt_opt = batch%dt; if (present(dt)) dt_opt = dt
    half_dt = 0.5_r8 * dt_opt

    associate (mesh   => m%mesh      , &
               u      => batch%u     , & ! in
               u_frac => batch%u_frac, & ! in
               v      => batch%v     , & ! in
               divx   => batch%divx  , & ! in
               divy   => batch%divy  , & ! in
               mx     => batch%qx    , & ! work array
               my     => batch%qy    , & ! work array
               dmxdt  => batch%qx    , & ! borrowed array
               dmydt  => batch%qy    )   ! borrowed array
    ! Run inner advective operators.
    ! FIXME: Avoid duplicated calculation.
    call batch%calc_cflxy_mass(half_dt)
    call hflx_mass(batch, u, u_frac, v, m, m, mfx, mfy, half_dt)
    call fill_halo(mfx, south_halo=.false., north_halo=.false., east_halo =.false.)
    call fill_halo(mfy, west_halo =.false., east_halo =.false., north_halo=.false.)
    ! Calculate intermediate tracer density due to advective operators.
    call divx_operator(mfx, dmxdt)
    call divy_operator(mfy, dmydt)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          mx%d(i,j,k) = (m%d(i,j,k) - half_dt * dmxdt%d(i,j,k)) / (1 - half_dt * divx%d(i,j,k))
          my%d(i,j,k) = (m%d(i,j,k) - half_dt * dmydt%d(i,j,k)) / (1 - half_dt * divy%d(i,j,k))
        end do
      end do
    end do
    call fill_halo(mx, west_halo =.false., east_halo =.false.)
    call fill_halo(my, south_halo=.false., north_halo=.false.)
    ! Run outer flux form operators.
    ! FIXME: Avoid duplicated calculation.
    call batch%calc_cflxy_mass(dt_opt)
    call hflx_mass(batch, u, u_frac, v, my, mx, mfx, mfy, dt_opt)
    call fill_halo(mfx, south_halo=.false., north_halo=.false., east_halo =.false.)
    call fill_halo(mfy, west_halo =.false., east_halo =.false., north_halo=.false.)
    end associate

    call perf_stop('ffsl_calc_mass_hflx_lin96')

  end subroutine ffsl_calc_mass_hflx_lin96

  subroutine ffsl_calc_tracer_hflx_lin96(batch, q, qmfx, qmfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: q
    type(latlon_field3d_type), intent(inout) :: qmfx
    type(latlon_field3d_type), intent(inout) :: qmfy
    real(r8), intent(in), optional :: dt

    integer ks, ke, i, j, k
    real(r8) dt_opt, half_dt

    call perf_start('ffsl_calc_tracer_hflx_lin96')

    dt_opt = batch%dt; if (present(dt)) dt_opt = dt
    half_dt = 0.5_r8 * dt_opt

    associate (mesh     => q%mesh        , &
               m        => batch%m       , & ! in
               u        => batch%u       , & ! in
               v        => batch%v       , & ! in
               mfx      => batch%mfx     , & ! in
               mfx_frac => batch%mfx_frac, & ! in
               mfy      => batch%mfy     , & ! in
               cflx     => batch%cflx    , & ! in
               cfly     => batch%cfly    , & ! in
               divx     => batch%divx    , & ! in
               divy     => batch%divy    , & ! in
               qx       => batch%qx      , & ! work array
               qy       => batch%qy      , & ! work array
               dqmxdt   => batch%qx      , & ! borrowed array
               dqmydt   => batch%qy      )   ! borrowed array
    ! Run inner advective operators.
    ! FIXME: Avoid duplicated calculation.
    call batch%calc_cflxy_tracer(m, m, mfx, mfy, cflx, cfly, mfx_frac, half_dt)
    call hflx_tracer(batch, m, m, cflx, cfly, mfx, mfx_frac, mfy, q, q, qmfx, qmfy, half_dt)
    call fill_halo(qmfx, south_halo=.false., north_halo=.false., east_halo =.false.)
    call fill_halo(qmfy, west_halo =.false., east_halo =.false., north_halo=.false.)
    select case (batch%loc)
    case ('cell', 'lev')
      ks = merge(mesh%full_kds, mesh%half_kds, batch%loc == 'cell')
      ke = merge(mesh%full_kde, mesh%half_kde, batch%loc == 'cell')
      ! Calculate intermediate tracer density due to advective operators.
      call divx_operator(qmfx, dqmxdt)
      call divy_operator(qmfy, dqmydt)
      do k = ks, ke
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            qx%d(i,j,k) = (q%d(i,j,k) - half_dt * dqmxdt%d(i,j,k) / m%d(i,j,k)) / (1 - half_dt * divx%d(i,j,k))
            qy%d(i,j,k) = (q%d(i,j,k) - half_dt * dqmydt%d(i,j,k) / m%d(i,j,k)) / (1 - half_dt * divy%d(i,j,k))
          end do
        end do
      end do
    case ('vtx')
    end select
    call fill_halo(qx, west_halo =.false., east_halo =.false.)
    call fill_halo(qy, south_halo=.false., north_halo=.false.)
    ! Run outer flux form operators.
    ! FIXME: Avoid duplicated calculation.
    call batch%calc_cflxy_tracer(m, m, mfx, mfy, cflx, cfly, mfx_frac, dt_opt)
    call hflx_tracer(batch, m, m, cflx, cfly, mfx, mfx_frac, mfy, qy, qx, qmfx, qmfy, dt_opt)
    call fill_halo(qmfx, south_halo=.false., north_halo=.false., east_halo =.false.)
    call fill_halo(qmfy, west_halo =.false., east_halo =.false., north_halo=.false.)
    end associate

    call perf_stop('ffsl_calc_tracer_hflx_lin96')

  end subroutine ffsl_calc_tracer_hflx_lin96

  subroutine ffsl_calc_mass_hflx_swift(batch, m, mfx, mfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: m
    type(latlon_field3d_type), intent(inout) :: mfx
    type(latlon_field3d_type), intent(inout) :: mfy
    real(r8), intent(in), optional :: dt

    integer i, j, k
    real(r8) dt_opt

    call perf_start('ffsl_calc_mass_hflx_swift')

    dt_opt = batch%dt; if (present(dt)) dt_opt = dt

    associate (mesh     => m%mesh        , &
               u        => batch%u       , & ! in
               u_frac   => batch%u_frac  , & ! in
               v        => batch%v       , & ! in
               cflx     => batch%cflx    , & ! inout
               cfly     => batch%cfly    , & ! inout
               divx     => batch%divx    , & ! in
               divy     => batch%divy    , & ! in
               mx       => batch%qx      , & ! work array
               my       => batch%qy      , & ! work array
               sx       => batch%qx      , & ! borrowed array
               sy       => batch%qy      , & ! borrowed array
               mfx0     => batch%qmfx0   , & ! work array
               mfy0     => batch%qmfy0   , & ! work array
               mfx_frac => batch%mfx_frac, & ! out
               dmxdt    => batch%qx      , & ! borrowed array
               dmydt    => batch%qy      )   ! borrowed array
    ! Run inner advective operators.
    call hflx_mass(batch, u, u_frac, v, m, m, mfx0, mfy0, dt_opt)
    call fill_halo(mfx0, south_halo=.false., north_halo=.false., east_halo =.false.)
    call fill_halo(mfy0, west_halo =.false., east_halo =.false., north_halo=.false.)
    ! Update cflx and cfly for outer step to use sy and sx as density.
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          sx%d(i,j,k) = 1 - dt_opt * divx%d(i,j,k) ! (33a)
          sy%d(i,j,k) = 1 - dt_opt * divy%d(i,j,k) ! (33b)
        end do
      end do
    end do
    call fill_halo(sx, west_halo =.false., east_halo =.false.)
    call fill_halo(sy, south_halo=.false., north_halo=.false.)
    ! NOTE: Swap sy and sx.
    call batch%calc_cflxy_tracer(sy, sx, u, v, cflx, cfly, u_frac, dt_opt)
    ! Calculate intermediate tracer density due to advective operators.
    call divx_operator(mfx0, dmxdt)
    call divy_operator(mfy0, dmydt)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          mx%d(i,j,k) = (m%d(i,j,k) - dt_opt * dmxdt%d(i,j,k)) / (1 - dt_opt * divx%d(i,j,k)) ! (34a)
          my%d(i,j,k) = (m%d(i,j,k) - dt_opt * dmydt%d(i,j,k)) / (1 - dt_opt * divy%d(i,j,k)) ! (34b)
        end do
      end do
    end do
    call fill_halo(mx, west_halo =.false., east_halo =.false.)
    call fill_halo(my, south_halo=.false., north_halo=.false.)
    ! Run outer flux form operators with updated cflx, cfly and u_frac.
    call hflx_mass(batch, u, u_frac, v, my, mx, mfx, mfy, dt_opt)
    ! Do SWIFT splitting.
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          mfx%d(i,j,k) = 0.5_r8 * (mfx0%d(i,j,k) + mfx%d(i,j,k)) ! (37a)
        end do
      end do
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          mfy%d(i,j,k) = 0.5_r8 * (mfy0%d(i,j,k) + mfy%d(i,j,k)) ! (37b)
        end do
      end do
    end do
    call fill_halo(mfx, south_halo=.false., north_halo=.false., east_halo =.false.)
    call fill_halo(mfy, west_halo =.false., east_halo =.false., north_halo=.false.)
    end associate

    call perf_stop('ffsl_calc_mass_hflx_swift')

  end subroutine ffsl_calc_mass_hflx_swift

  subroutine swift_hflx_connector(batch, m, mfx, mfy, mx, my, cflx, cfly, mfx_frac, dt)

    type(adv_batch_type), intent(inout) :: batch
    type(latlon_field3d_type), intent(in) :: m
    type(latlon_field3d_type), intent(in) :: mfx
    type(latlon_field3d_type), intent(in) :: mfy
    type(latlon_field3d_type), intent(inout) :: mx
    type(latlon_field3d_type), intent(inout) :: my
    type(latlon_field3d_type), intent(inout) :: cflx
    type(latlon_field3d_type), intent(inout) :: cfly
    type(latlon_field3d_type), intent(inout) :: mfx_frac
    real(r8), intent(in), optional :: dt

    integer ks, ke, i, j, k
    real(r8) dt_opt

    call perf_start('swift_hflx_connector')

    dt_opt = batch%dt; if (present(dt)) dt_opt = dt

    associate (mesh  => batch%mesh, &
               dmxdt => batch%qx  , & ! borrowed array
               dmydt => batch%qy  )   ! borrowed array
    ks = merge(mesh%full_kds, mesh%half_kds, batch%loc(1:3) /= 'lev')
    ke = merge(mesh%full_kde, mesh%half_kde, batch%loc(1:3) /= 'lev')
    ! Calculate new mx and my from final mfx and mfy.
    call divx_operator(mfx, dmxdt)
    call divy_operator(mfy, dmydt)
    do k = ks, ke
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          mx%d(i,j,k) = m%d(i,j,k) - dt_opt * dmxdt%d(i,j,k) ! (37c)
          my%d(i,j,k) = m%d(i,j,k) - dt_opt * dmydt%d(i,j,k) ! (37d)
        end do
      end do
    end do
    call fill_halo(mx, west_halo =.false., east_halo =.false.)
    call fill_halo(my, south_halo=.false., north_halo=.false.)
    ! Update cflx, cfly and mfx_frac for tracer advection.
    ! NOTE: Swap mx and my.
    call batch%calc_cflxy_tracer(my, mx, mfx, mfy, cflx, cfly, mfx_frac, dt_opt)
    end associate

    call perf_stop('swift_hflx_connector')

  end subroutine swift_hflx_connector

  subroutine ffsl_calc_tracer_hflx_swift(batch, q, qmfx, qmfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: q
    type(latlon_field3d_type), intent(inout) :: qmfx
    type(latlon_field3d_type), intent(inout) :: qmfy
    real(r8), intent(in), optional :: dt

    integer ks, ke, i, j, k
    real(r8) dt_opt

    call perf_start('ffsl_calc_tracer_hflx_swift')

    dt_opt = batch%dt; if (present(dt)) dt_opt = dt

    associate (mesh        => q%mesh           , &
               m           => batch%m          , & ! in
               mx          => batch%bg%qx      , & ! borrowed array
               my          => batch%bg%qy      , & ! borrowed array
               u           => batch%u          , & ! in
               v           => batch%v          , & ! in
               mfx         => batch%mfx        , & ! in
               mfx_frac    => batch%mfx_frac   , & ! in
               mfx_my_frac => batch%bg%mfx_frac, & ! borrowed array
               mfy         => batch%mfy        , & ! in
               cflx        => batch%cflx       , & ! in
               cflx_my     => batch%bg%cflx    , & ! borrowed array
               cfly        => batch%cfly       , & ! in
               cfly_mx     => batch%bg%cfly    , & ! borrowed array
               divx        => batch%divx       , & ! in
               divy        => batch%divy       , & ! in
               qx          => batch%qx         , & ! work array
               qy          => batch%qy         , & ! work array
               qmfx0       => batch%qmfx0      , & ! out
               qmfy0       => batch%qmfy0      , & ! out
               dqmxdt      => batch%qx         , & ! borrowed array
               dqmydt      => batch%qy         )   ! borrowed array
    ! FIXME: This should only run once for each batch.
    call swift_hflx_connector(batch, m, mfx, mfy, mx, my, cflx_my, cfly_mx, mfx_my_frac, dt_opt)
    ! Run inner advective operators.
    call hflx_tracer(batch, m, m, cflx, cfly, mfx, mfx_frac, mfy, q, q, qmfx0, qmfy0, dt_opt)
    call fill_halo(qmfx0, south_halo=.false., north_halo=.false., east_halo =.false.)
    call fill_halo(qmfy0, west_halo =.false., east_halo =.false., north_halo=.false.)
    select case (batch%loc)
    case ('cell', 'lev')
      ks = merge(mesh%full_kds, mesh%half_kds, batch%loc == 'cell')
      ke = merge(mesh%full_kde, mesh%half_kde, batch%loc == 'cell')
      ! Calculate intermediate tracer density due to advective operators.
      call divx_operator(qmfx0, dqmxdt)
      call divy_operator(qmfy0, dqmydt)
      do k = ks, ke
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            qx%d(i,j,k) = (q%d(i,j,k) * m%d(i,j,k) - dt_opt * dqmxdt%d(i,j,k)) / mx%d(i,j,k) ! (38a, 39a)
            qy%d(i,j,k) = (q%d(i,j,k) * m%d(i,j,k) - dt_opt * dqmydt%d(i,j,k)) / my%d(i,j,k) ! (38b, 39b)
          end do
        end do
      end do
    case ('vtx')
    end select
    call fill_halo(qx, west_halo =.false., east_halo =.false.)
    call fill_halo(qy, south_halo=.false., north_halo=.false.)
    ! Run outer flux form operators.
    call hflx_tracer(batch, my, mx, cflx_my, cfly_mx, mfx, mfx_my_frac, mfy, qy, qx, qmfx, qmfy, dt_opt)
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
    call fill_halo(qmfx, south_halo=.false., north_halo=.false., east_halo =.false.)
    call fill_halo(qmfy, west_halo =.false., east_halo =.false., north_halo=.false.)
    end associate

    call perf_stop('ffsl_calc_tracer_hflx_swift')

  end subroutine ffsl_calc_tracer_hflx_swift

  subroutine ffsl_calc_mass_vflx(batch, m, mfz, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: m
    type(latlon_field3d_type), intent(inout) :: mfz
    real(r8), intent(in), optional :: dt

    real(r8) dt_opt

    call perf_start('ffsl_calc_mass_vflx')

    dt_opt = batch%dt; if (present(dt)) dt_opt = dt

    call vflx_mass(batch, batch%w, batch%w_frac, m, mfz, dt_opt)

    call perf_stop('ffsl_calc_mass_vflx')

  end subroutine ffsl_calc_mass_vflx

  subroutine ffsl_calc_tracer_vflx(batch, q, qmfz, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: q
    type(latlon_field3d_type), intent(inout) :: qmfz
    real(r8), intent(in), optional :: dt

    real(r8) dt_opt

    call perf_start('ffsl_calc_tracer_vflx')

    dt_opt = batch%dt; if (present(dt)) dt_opt = dt

    call vflx_tracer(batch, batch%mfz, batch%mfz_frac, q, qmfz, dt_opt)

    call perf_stop('ffsl_calc_tracer_vflx')

  end subroutine ffsl_calc_tracer_vflx

  subroutine hflx_ppm_mass(batch, u, u_frac, v, mx, my, mfx, mfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: u
    type(latlon_field3d_type), intent(in   ) :: u_frac
    type(latlon_field3d_type), intent(in   ) :: v
    type(latlon_field3d_type), intent(in   ) :: mx
    type(latlon_field3d_type), intent(in   ) :: my
    type(latlon_field3d_type), intent(inout) :: mfx
    type(latlon_field3d_type), intent(inout) :: mfy
    real(r8), intent(in) :: dt

    integer ks, ke, i, j, k, iu, ju, ci
    real(r8) cf, mr

    associate (mesh => u%mesh    , &
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
              mfx%d(i,j,k) = 0
            else if (cflx%d(i,j,k) > 0) then
              iu = i - ci
              mr = ppm3(cf, mx%d(iu-2,j,k), mx%d(iu-1,j,k), mx%d(iu,j,k), mx%d(iu+1,j,k), mx%d(iu+2,j,k))
              mfx%d(i,j,k) = u_frac%d(i,j,k) * mr + sum(mx%d(iu+1:i,j,k)) * mesh%de_lon(j) / dt
            else
              iu = i - ci + 1
              mr = ppm3(cf, mx%d(iu-2,j,k), mx%d(iu-1,j,k), mx%d(iu,j,k), mx%d(iu+1,j,k), mx%d(iu+2,j,k))
              mfx%d(i,j,k) = u_frac%d(i,j,k) * mr - sum(mx%d(i+1:iu-1,j,k)) * mesh%de_lon(j) / dt
            end if
          end do
        end do
        ! Along y-axis
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            if (abs(cfly%d(i,j,k)) < 1.0e-16_r8) then
              mfy%d(i,j,k) = 0
            else if (cfly%d(i,j,k) > 0) then
              ju = j
              mr = ppm3(cfly%d(i,j,k), my%d(i,ju-2,k), my%d(i,ju-1,k), my%d(i,ju,k), my%d(i,ju+1,k), my%d(i,ju+2,k))
              mfy%d(i,j,k) = v%d(i,j,k) * mr
            else if (cfly%d(i,j,k) < 0) then
              ju = j + 1
              mr = ppm3(cfly%d(i,j,k), my%d(i,ju-2,k), my%d(i,ju-1,k), my%d(i,ju,k), my%d(i,ju+1,k), my%d(i,ju+2,k))
              mfy%d(i,j,k) = v%d(i,j,k) * mr
            end if
          end do
        end do
      end do
    case ('vtx')
    end select
    end associate

  end subroutine hflx_ppm_mass

  subroutine hflx_ppm_tracer(batch, mx, my, cflx, cfly, mfx, mfx_frac, mfy, qx, qy, qmfx, qmfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: mx
    type(latlon_field3d_type), intent(in   ) :: my
    type(latlon_field3d_type), intent(in   ) :: cflx
    type(latlon_field3d_type), intent(in   ) :: cfly
    type(latlon_field3d_type), intent(in   ) :: mfx
    type(latlon_field3d_type), intent(in   ) :: mfx_frac
    type(latlon_field3d_type), intent(in   ) :: mfy
    type(latlon_field3d_type), intent(in   ) :: qx
    type(latlon_field3d_type), intent(in   ) :: qy
    type(latlon_field3d_type), intent(inout) :: qmfx
    type(latlon_field3d_type), intent(inout) :: qmfy
    real(r8), intent(in) :: dt

    integer ks, ke, i, j, k, iu, ju, ci
    real(r8) cf, qr

    associate (mesh => mfx%mesh)
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
              qr = ppm3(cf, qx%d(iu-2,j,k), qx%d(iu-1,j,k), qx%d(iu,j,k), qx%d(iu+1,j,k), qx%d(iu+2,j,k))
              qmfx%d(i,j,k) = mfx_frac%d(i,j,k) * qr + sum(mx%d(iu+1:i,j,k) * qx%d(iu+1:i,j,k)) * mesh%de_lon(j) / dt
            else
              iu = i - ci + 1
              qr = ppm3(cf, qx%d(iu-2,j,k), qx%d(iu-1,j,k), qx%d(iu,j,k), qx%d(iu+1,j,k), qx%d(iu+2,j,k))
              qmfx%d(i,j,k) = mfx_frac%d(i,j,k) * qr - sum(mx%d(i+1:iu-1,j,k) * qx%d(i+1:iu-1,j,k)) * mesh%de_lon(j) / dt
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
              qr = ppm3(cfly%d(i,j,k), qy%d(i,ju-2,k), qy%d(i,ju-1,k), qy%d(i,ju,k), qy%d(i,ju+1,k), qy%d(i,ju+2,k))
              qmfy%d(i,j,k) = mfy%d(i,j,k) * qr
            else if (cfly%d(i,j,k) < 0) then
              ju = j + 1
              qr = ppm3(cfly%d(i,j,k), qy%d(i,ju-2,k), qy%d(i,ju-1,k), qy%d(i,ju,k), qy%d(i,ju+1,k), qy%d(i,ju+2,k))
              qmfy%d(i,j,k) = mfy%d(i,j,k) * qr
            end if
          end do
        end do
      end do
    case ('vtx')
    end select
    end associate

  end subroutine hflx_ppm_tracer

  subroutine vflx_ppm_mass(batch, w, w_frac, m, mfz, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: w
    type(latlon_field3d_type), intent(in   ) :: w_frac
    type(latlon_field3d_type), intent(in   ) :: m
    type(latlon_field3d_type), intent(inout) :: mfz
    real(r8), intent(in) :: dt

    integer i, j, k, ku, ci
    real(r8) cf, mr

    associate (mesh => m%mesh    , &
               cflz => batch%cflz)   ! in
    select case (batch%loc)
    case ('cell')
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            ci = int(cflz%d(i,j,k))
            cf = cflz%d(i,j,k) - ci
            if (abs(cflz%d(i,j,k)) < 1.0e-16_r8) then
              mfz%d(i,j,k) = 0
            else if (cflz%d(i,j,k) > 0) then
              ku = k - ci - 1
              mr = ppm3(cf, m%d(i,j,ku-2), m%d(i,j,ku-1), m%d(i,j,ku), m%d(i,j,ku+1), m%d(i,j,ku+2))
              mfz%d(i,j,k) = w_frac%d(i,j,k) * mr + sum(m%d(i,j,ku+1:k-1)) / dt
            else
              ku = k - ci
              mr = ppm3(cf, m%d(i,j,ku-2), m%d(i,j,ku-1), m%d(i,j,ku), m%d(i,j,ku+1), m%d(i,j,ku+2))
              mfz%d(i,j,k) = w_frac%d(i,j,k) * mr - sum(m%d(i,j,k:ku-1)) / dt
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
              mfz%d(i,j,k) = 0
            else if (cflz%d(i,j,k) > 0) then
              ku = k - ci
              mr = ppm3(cf, m%d(i,j,ku-2), m%d(i,j,ku-1), m%d(i,j,ku), m%d(i,j,ku+1), m%d(i,j,ku+2))
              mfz%d(i,j,k) = w_frac%d(i,j,k) * mr + sum(m%d(i,j,ku+1:k)) / dt
            else
              ku = k - ci + 1
              mr = ppm3(cf, m%d(i,j,ku-2), m%d(i,j,ku-1), m%d(i,j,ku), m%d(i,j,ku+1), m%d(i,j,ku+2))
              mfz%d(i,j,k) = w_frac%d(i,j,k) * mr - sum(m%d(i,j,k+1:ku-1)) / dt
            end if
          end do
        end do
      end do
    end select
    end associate

  end subroutine vflx_ppm_mass

  subroutine vflx_ppm_tracer(batch, mfz, mfz_frac, q, qmfz, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: mfz
    type(latlon_field3d_type), intent(in   ) :: mfz_frac
    type(latlon_field3d_type), intent(in   ) :: q
    type(latlon_field3d_type), intent(inout) :: qmfz
    real(r8), intent(in) :: dt

    integer i, j, k, ku, ci
    real(r8) cf, qr

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
              qr = ppm3(cf, q%d(i,j,ku-2), q%d(i,j,ku-1), q%d(i,j,ku), q%d(i,j,ku+1), q%d(i,j,ku+2))
              qmfz%d(i,j,k) = mfz_frac%d(i,j,k) * qr + sum(m%d(i,j,ku+1:k-1) * q%d(i,j,ku+1:k-1)) / dt
            else
              ku = k - ci
              qr = ppm3(cf, q%d(i,j,ku-2), q%d(i,j,ku-1), q%d(i,j,ku), q%d(i,j,ku+1), q%d(i,j,ku+2))
              qmfz%d(i,j,k) = mfz_frac%d(i,j,k) * qr - sum(m%d(i,j,k:ku-1) * q%d(i,j,k:ku-1)) / dt
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
              qr = ppm3(cf, q%d(i,j,ku-2), q%d(i,j,ku-1), q%d(i,j,ku), q%d(i,j,ku+1), q%d(i,j,ku+2))
              qmfz%d(i,j,k) = mfz_frac%d(i,j,k) * qr + sum(m%d(i,j,ku+1:k) * q%d(i,j,ku+1:k)) / dt
            else
              ku = k - ci + 1
              qr = ppm3(cf, q%d(i,j,ku-2), q%d(i,j,ku-1), q%d(i,j,ku), q%d(i,j,ku+1), q%d(i,j,ku+2))
              qmfz%d(i,j,k) = mfz_frac%d(i,j,k) * qr - sum(m%d(i,j,k+1:ku-1) * q%d(i,j,k+1:ku-1)) / dt
            end if
          end do
        end do
      end do
    end select
    end associate

  end subroutine vflx_ppm_tracer

end module ffsl_mod
