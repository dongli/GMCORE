! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module adv_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use math_mod
  use time_mod
  use block_mod
  use latlon_field_types_mod
  use latlon_operators_mod
  use latlon_parallel_mod
  use process_mod, only: proc
  use tracer_mod
  use interp_mod
  use adv_batch_mod
  use ffsl_mod
  use upwind_mod
  use weno_mod
  use physics_mod
  use perf_mod

  implicit none

  private

  public adv_init
  public adv_prepare
  public adv_run_mass
  public adv_run_tracers
  public adv_final
  public adv_fill_vhalo
  public adv_set_wind
  public adv_accum_wind
  public adv_calc_tracer_hflx
  public adv_calc_tracer_vflx
  public adv_batch_type

  public upwind1
  public upwind3
  public upwind5
  public weno3
  public weno5

contains

  subroutine adv_init()

    integer iblk, ibat, itra, n, idx(1000)

    call adv_final()

    call ffsl_init()

    ! Initialize advection batches.
    do iblk = 1, size(blocks)
      if (advection) then
        call blocks(iblk)%adv_batch_bg%init(                  &
          blocks(iblk)%big_filter                           , &
          blocks(iblk)%filter_mesh, blocks(iblk)%filter_halo, &
          blocks(iblk)%mesh, blocks(iblk)%halo              , &
          bg_adv_scheme, 'cell', 'bg', dt_dyn, dynamic=.true., slave=.false.)
      end if
      if (baroclinic) then
        call blocks(iblk)%adv_batch_pt%init(                  &
          blocks(iblk)%big_filter                           , &
          blocks(iblk)%filter_mesh, blocks(iblk)%filter_halo, &
          blocks(iblk)%mesh, blocks(iblk)%halo              , &
          pt_adv_scheme, 'cell', 'pt', dt_dyn, dynamic=.true., slave=.true.)
      end if
      if (nonhydrostatic) then
        call blocks(iblk)%adv_batch_nh%init(                  &
          blocks(iblk)%big_filter                           , &
          blocks(iblk)%filter_mesh, blocks(iblk)%filter_halo, &
          blocks(iblk)%mesh, blocks(iblk)%halo              , &
          nh_adv_scheme, 'lev', 'nh', dt_dyn, dynamic=.true., slave=.true.)
      end if
    end do

    if (nbatches == 0) then
      if (proc%is_root()) call log_warning('No advection batches have been defined yet!')
      return
    end if

    do iblk = 1, size(blocks)
      allocate(blocks(iblk)%adv_batches(nbatches))
      do ibat = 1, nbatches
        n = 0
        do itra = 1, ntracers
          if (batch_names(ibat) == tracer_batches(itra)) then
            n = n + 1
            idx(n) = itra
          end if
        end do
        call blocks(iblk)%adv_batches(ibat)%init(             &
          blocks(iblk)%big_filter                           , &
          blocks(iblk)%filter_mesh, blocks(iblk)%filter_halo, &
          blocks(iblk)%mesh, blocks(iblk)%halo              , &
          'ffsl', 'cell', batch_names(ibat), batch_dts(ibat), &
          dynamic=.false., slave=.true., idx=idx(1:n)       , &
          bg=blocks(iblk)%adv_batch_bg)
      end do
    end do

    if (proc%is_root()) then
      call log_notice('There are ' // to_str(size(blocks(1)%adv_batches)) // ' advection batches.')
      do ibat = 1, size(blocks(1)%adv_batches)
        write(*, *) '- ', trim(blocks(1)%adv_batches(ibat)%name), int(blocks(1)%adv_batches(ibat)%dt), blocks(1)%adv_batches(ibat)%ntracers
      end do
    end if

  end subroutine adv_init

  subroutine adv_prepare(itime)

    integer, intent(in) :: itime

    integer iblk, m

    call perf_start('adv_prepare')

    if (.not. restart) then
      do iblk = 1, size(blocks)
        associate (dmg => blocks(iblk)%dstate(itime)%dmg)
        if (allocated(blocks(iblk)%adv_batches)) then
          do m = 1, size(blocks(iblk)%adv_batches)
            call blocks(iblk)%adv_batches(m)%copy_m_old(dmg)
          end do
        end if
        end associate
      end do
      call adv_accum_wind(itime)
    end if

    call perf_stop('adv_prepare')

  end subroutine adv_prepare

  subroutine adv_calc_mass_hflx(batch, m, mfx, mfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: m
    type(latlon_field3d_type), intent(inout) :: mfx
    type(latlon_field3d_type), intent(inout) :: mfy
    real(r8), intent(in), optional :: dt

    select case (batch%scheme_h)
    case ('ffsl')
      call ffsl_calc_mass_hflx_swift1(batch, m, mfx, mfy, dt)
    end select

  end subroutine adv_calc_mass_hflx

  subroutine adv_calc_mass_vflx(batch, m, mfz, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: m
    type(latlon_field3d_type), intent(inout) :: mfz
    real(r8), intent(in), optional :: dt

    select case (batch%scheme_v)
    case ('ffsl')
      call ffsl_calc_mass_vflx(batch, m, mfz, dt)
    end select

  end subroutine adv_calc_mass_vflx

  subroutine adv_calc_tracer_hflx(batch, q, qmfx, qmfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: q
    type(latlon_field3d_type), intent(inout) :: qmfx
    type(latlon_field3d_type), intent(inout) :: qmfy
    real(r8), intent(in), optional :: dt

    select case (batch%scheme_h)
    case ('upwind')
      call upwind_calc_tracer_hflx(batch, q, qmfx, qmfy, dt)
    case ('ffsl')
      call ffsl_calc_tracer_hflx_swift1(batch, q, qmfx, qmfy, dt)
    end select

  end subroutine adv_calc_tracer_hflx

  subroutine adv_calc_tracer_vflx(batch, q, qmfz, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: q
    type(latlon_field3d_type), intent(inout) :: qmfz
    real(r8), intent(in), optional :: dt

    select case (batch%scheme_v)
    case ('upwind')
      call upwind_calc_tracer_vflx(batch, q, qmfz, dt)
    case ('ffsl')
      call ffsl_calc_tracer_vflx(batch, q, qmfz, dt)
    end select

  end subroutine adv_calc_tracer_vflx

  subroutine adv_run_mass(old, new, dt)

    integer, intent(in) :: old
    integer, intent(in) :: new
    real(r8), intent(in), optional :: dt

    real(r8) dt_opt
    integer iblk, i, j, k

    if (.not. blocks(1)%adv_batch_bg%initialized) return

    dt_opt = blocks(1)%adv_batch_bg%dt; if (present(dt)) dt_opt = dt

    call adv_set_wind(old)

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)                 , &
                 batch => blocks(iblk)%adv_batch_bg    , &
                 mesh  => blocks(iblk)%mesh            , &
                 m_old => blocks(iblk)%dstate(old)%dmg , & ! in
                 m_new => blocks(iblk)%dstate(new)%dmg , & ! out
                 mfx   => blocks(iblk)%adv_batch_bg%mfx, & ! work array
                 mfy   => blocks(iblk)%adv_batch_bg%mfy, & ! work array
                 mfz   => blocks(iblk)%adv_batch_bg%mfz, & ! work array
                 dmdt  => blocks(iblk)%aux%pv          )   ! borrowed array
      call adv_calc_mass_hflx(batch, m_old, mfx, mfy, dt_opt)
      call div_operator(mfx, mfy, dmdt)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            m_new%d(i,j,k) = m_old%d(i,j,k) - dt_opt * dmdt%d(i,j,k)
          end do
        end do
      end do
      call adv_fill_vhalo(m_new)
      call adv_calc_mass_vflx(batch, m_new, mfz, dt_opt)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            m_new%d(i,j,k) = m_new%d(i,j,k) - dt_opt * (mfz%d(i,j,k+1) - mfz%d(i,j,k))
          end do
        end do
      end do
      call fill_halo(m_new)
      end associate
    end do

  end subroutine adv_run_mass

  subroutine adv_run_tracers(itime)

    integer, intent(in) :: itime

    integer iblk, i, j, k, l, m, idx
    type(latlon_field3d_type) q_old, q_new

    if (.not. allocated(blocks(1)%adv_batches)) return

    call adv_accum_wind(itime)

    do iblk = 1, size(blocks)
      associate (block     => blocks(iblk)                  , &
                 mesh      => blocks(iblk)%filter_mesh      , &
                 m_new     => blocks(iblk)%dstate(itime)%dmg, & ! in
                 dqdt      => blocks(iblk)%aux%pv           )   ! borrowed array
      do m = 1, size(block%adv_batches)
        if (time_is_alerted(block%adv_batches(m)%name)) then
          if (m == 1 .and. pdc_type == 2) call physics_update_dynamics(block, itime, dt_adv)
          associate (batch => block%adv_batches(m))
          call q_old%link(batch%dmf) ! Borrow array.
          do l = 1, block%adv_batches(m)%ntracers
            idx = batch%idx(l)
            call q_new%link(tracers(iblk)%q, idx)
            q_old%d = q_new%d
            associate (m_old => batch%m   , & ! in
                       qmfx  => batch%qmfx, & ! work array
                       qmfy  => batch%qmfy, & ! work array
                       qmfz  => batch%qmfz)   ! work array
            ! Calculate horizontal tracer mass flux.
            call adv_calc_tracer_hflx(batch, q_old, qmfx, qmfy)
            call div_operator(qmfx, qmfy, dqdt)
            ! Update tracer mixing ratio due to horizontal advection.
            do k = mesh%full_kds, mesh%full_kde
              do j = mesh%full_jds, mesh%full_jde
                do i = mesh%full_ids, mesh%full_ide
                  q_new%d(i,j,k) = (m_old%d(i,j,k) * q_old%d(i,j,k) - dt_adv * dqdt%d(i,j,k)) / m_new%d(i,j,k)
                end do
              end do
            end do
            ! Calculate vertical tracer mass flux.
            call adv_fill_vhalo(q_new)
            call adv_calc_tracer_vflx(block%adv_batches(m), q_new, qmfz)
            do k = mesh%full_kds, mesh%full_kde
              do j = mesh%full_jds, mesh%full_jde
                do i = mesh%full_ids, mesh%full_ide
                  q_new%d(i,j,k) = q_new%d(i,j,k) - dt_adv * (qmfz%d(i,j,k+1) - qmfz%d(i,j,k)) / m_new%d(i,j,k)
                end do
              end do
            end do
            call fill_halo(q_new)
            end associate
          end do
          end associate
          call block%adv_batches(m)%copy_m_old(m_new)
        end if
      end do
      call tracer_calc_qm(block)
      end associate
    end do

  end subroutine adv_run_tracers

  subroutine adv_set_wind(itime)

    integer, intent(in) :: itime

    integer iblk

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk))
      if (block%adv_batch_bg%initialized) then
        call block%adv_batch_bg%set_wind( &
          u=block%dstate(itime)%u_lon   , & ! in
          v=block%dstate(itime)%v_lat   , & ! in
          w=block%dstate(itime)%we_lev  , & ! in
          mfx=block%aux%mfx_lon         , & ! out
          mfy=block%aux%mfy_lat         )   ! out
      end if
      end associate
    end do

  end subroutine adv_set_wind

  subroutine adv_accum_wind(itime)

    integer, intent(in) :: itime

    integer iblk, l

    call perf_start('adv_accum_wind')

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)                  , &
        m     => blocks(iblk)%dstate(itime)%dmg, & ! in
        mx    => blocks(iblk)%aux%dmg_lon      , & ! out
        my    => blocks(iblk)%aux%dmg_lat      , & ! out
        mfx   => blocks(iblk)%aux%mfx_lon      , & ! in
        mfy   => blocks(iblk)%aux%mfy_lat      )   ! in
      if (advection) then
        call interp_run(m, mx)
        call fill_halo(mx)
        call interp_run(m, my)
        call fill_halo(my)
      end if
      if (allocated(block%adv_batches)) then
        do l = 1, size(block%adv_batches)
          select case (block%adv_batches(l)%loc)
          case ('cell')
            call block%adv_batches(l)%accum_wind(mx, my, mfx, mfy)
          end select
        end do
      end if
      end associate
    end do

    call perf_stop('adv_accum_wind')

  end subroutine adv_accum_wind

  subroutine adv_final()

  end subroutine adv_final

end module adv_mod
