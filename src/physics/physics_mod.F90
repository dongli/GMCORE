! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module physics_mod

  use fiona
  use const_mod
  use namelist_mod
  use time_mod
  use formula_mod
  use block_mod
  use tracer_mod
  use physics_types_mod
  use dp_coupling_mod
  use latlon_parallel_mod
  use simple_physics_driver_mod, simple_objects => objects
#ifdef HAS_CAM
  use cam_physics_driver_mod, cam_objects => objects
#endif
  use mars_nasa_physics_driver_mod
  use perf_mod

  implicit none

  private

  public physics_init_stage1
  public physics_init_stage2
  public physics_init_stage3
  public physics_run
  public physics_update
  public physics_update_dynamics
  public physics_final
  public physics_add_output
  public physics_output

  type(physics_mesh_type), allocatable :: mesh(:)

contains

  subroutine physics_init_stage1(namelist_path)

    character(*), intent(in) :: namelist_path

    select case (physics_suite)
    case ('simple_physics')
      call simple_physics_init_stage1(namelist_path, dt_adv, dt_phys)
    end select

  end subroutine physics_init_stage1

  subroutine physics_init_stage2(namelist_path)

    character(*), intent(in) :: namelist_path

    integer nblk, iblk, icol, i, j
    integer , allocatable :: ncol(:)
    real(r8), allocatable :: lon (:,:)
    real(r8), allocatable :: lat (:,:)
    real(r8), allocatable :: area(:,:)

    if (physics_suite /= 'N/A') then
      call time_add_alert('phys', seconds=dt_phys)
    end if

    nblk = size(blocks)
    allocate(ncol(nblk))
    ncol = blocks(:)%mesh%full_nlon * blocks(:)%mesh%full_nlat
    allocate(lon (maxval(ncol),nblk))
    allocate(lat (maxval(ncol),nblk))
    allocate(area(maxval(ncol),nblk))
    do iblk = 1, nblk
      ! This is specific to lat-lon mesh. Maybe it should be put into mesh type.
      associate (mesh => blocks(iblk)%mesh)
      icol = 1
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          lon (icol,iblk) = mesh%full_lon (i)
          lat (icol,iblk) = mesh%full_lat (j)
          area(icol,iblk) = mesh%area_cell(j)
          icol = icol + 1
        end do
      end do
      ncol(iblk) = icol - 1
      end associate
    end do
    ! Create mesh objects.
    allocate(mesh(nblk))
    do iblk = 1, nblk
      associate (dmesh => blocks(iblk)%mesh)
      call mesh(iblk)%init(ncol(iblk), nlev, lon(:,iblk), lat(:,iblk), &
                           blocks(iblk)%mesh%full_lev, blocks(iblk)%mesh%full_dlev, &
                           area(:,iblk), ptop=ptop)
      ! For output dimensions
      mesh(iblk)%cell_start_2d = [dmesh%full_ids ,dmesh%full_jds ,1]
      mesh(iblk)%cell_count_2d = [dmesh%full_nlon,dmesh%full_nlat,1]
      mesh(iblk)%cell_start_3d = [dmesh%full_ids ,dmesh%full_jds ,dmesh%full_kds ,1]
      mesh(iblk)%cell_count_3d = [dmesh%full_nlon,dmesh%full_nlat,dmesh%full_nlev,1]
      mesh(iblk)%lev_start = [dmesh%full_ids ,dmesh%full_jds ,dmesh%half_kds ,1]
      mesh(iblk)%lev_count = [dmesh%full_nlon,dmesh%full_nlat,dmesh%half_nlev,1]
      end associate
    end do
    deallocate(ncol, lon, lat, area)

    select case (physics_suite)
    case ('simple_physics')
      call simple_physics_init_stage2(mesh)
    case ('cam')
#ifdef HAS_CAM
      call cam_physics_init_stage2(namelist_path, mesh, dt_adv, dt_phys)
#else
      if (proc%is_root()) call log_error('CAM physics is not compiled!')
#endif
    case ('mars_nasa')
      call mars_nasa_init_stage2(namelist_path, mesh, dt_adv, dt_phys, &
        min_lon, max_lon, min_lat, max_lat, input_ngroups, gmcore_root)
    end select

  end subroutine physics_init_stage2

  subroutine physics_init_stage3()

    integer iblk

    do iblk = 1, size(blocks)
      call dp_coupling_d2p(blocks(iblk), 1)
    end do

    select case (physics_suite)
    case ('mars_nasa')
      call mars_nasa_init_stage3()
    end select

  end subroutine physics_init_stage3

  subroutine physics_run(block, itime, dt)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime
    real(r8), intent(in) :: dt

    if (.not. time_is_alerted('phys')) return

    call perf_start('physics_run')

    if (proc%is_root()) call log_notice('Run ' // to_upper(trim(physics_suite)) // ' physics.')

    call dp_coupling_d2p(block, itime)

    select case (physics_suite)
    case ('simple_physics')
      call simple_physics_run()
#ifdef HAS_CAM
    case ('cam')
      call cam_physics_run_stage1()
      call cam_physics_sfc_flux()
      call cam_physics_run_stage2()
#endif
    case ('mars_nasa')
      call mars_nasa_run(curr_time)
    end select

    call dp_coupling_p2d(block, itime)

    select case (physics_suite)
    case ('simple_physics')
      call physics_set_tracers(block, simple_objects(block%id)%state)
#ifdef HAS_CAM
    case ('cam')
      call physics_set_tracers(block, cam_objects(block%id)%state)
#endif
    end select

    call perf_stop('physics_run')

  contains

    subroutine physics_set_tracers(block, pstate)

      type(block_type), intent(in) :: block
      class(physics_state_type), intent(in) :: pstate

      integer i, j, k, m, icol

      associate (mesh => block%mesh, q => tracers(block%id)%q, qm => tracers(block%id)%qm)
      qm%d = 0
      do m = 1, ntracers
        if (is_water_tracer(m)) then
          do k = mesh%full_kds, mesh%full_kde
            icol = 1
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%full_ids, mesh%full_ide
                qm%d(i,j,k) = pstate%q(icol,k,m)
                icol = icol + 1
              end do
            end do
          end do
        end if
      end do
      do m = 1, ntracers
        do k = mesh%full_kds, mesh%full_kde
          icol = 1
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              if (physics_use_wet_tracers(m)) then
                q%d(i,j,k,m) = pstate%q(icol,k,m) / (1 - qm%d(i,j,k))
              else
                q%d(i,j,k,m) = pstate%q(icol,k,m)
              end if
              icol = icol + 1
            end do
          end do
        end do
        call tracer_fill_negative_values(block, itime, q%d(:,:,:,m), m, __FILE__, __LINE__)
        call fill_halo(q, m)
      end do
      call tracer_calc_qm(block)
      end associate

    end subroutine physics_set_tracers

  end subroutine physics_run

  subroutine physics_update(block, itime, dt)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime
    real(r8), intent(in) :: dt

    class(physics_tend_type), pointer :: tend
    integer i, j, k, m

    call perf_start('physics_update')

    associate (mesh  => block%mesh               , &
               dudt  => block%aux%dudt_phys      , & ! in
               dvdt  => block%aux%dvdt_phys      , & ! in
               dptdt => block%aux%dptdt_phys     , & ! in
               dqdt  => block%aux%dqdt_phys      , & ! in
               dmg   => block%dstate(itime)%dmg  , & ! in
               u_lon => block%dstate(itime)%u_lon, & ! inout
               v_lat => block%dstate(itime)%v_lat, & ! inout
               pt    => block%dstate(itime)%pt   , & ! inout
               q     => tracers(block%id)%q      )   ! inout
    ! Update dynamics.
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          u_lon%d(i,j,k) = u_lon%d(i,j,k) + dt * 0.5_r8 * (dudt%d(i,j,k) + dudt%d(i+1,j,k))
        end do
      end do
    end do
    call fill_halo(u_lon)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          v_lat%d(i,j,k) = v_lat%d(i,j,k) + dt * 0.5_r8 * (dvdt%d(i,j,k) + dvdt%d(i,j+1,k))
        end do
      end do
    end do
    call fill_halo(v_lat)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          pt%d(i,j,k) = pt%d(i,j,k) + dt * dptdt%d(i,j,k) / dmg%d(i,j,k)
        end do
      end do
    end do
    call fill_halo(pt)
    ! Update tracers.
    do m = 1, ntracers
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            q%d(i,j,k,m) = q%d(i,j,k,m) + dt * dqdt%d(i,j,k,m) / dmg%d(i,j,k)
          end do
        end do
      end do
      call tracer_fill_negative_values(block, itime, q%d(:,:,:,m), m, __FILE__, __LINE__)
      call fill_halo(q, m)
    end do
    call tracer_calc_qm(block)
    end associate

    call perf_stop('physics_update')

  end subroutine physics_update

  subroutine physics_update_dynamics(block, itime, dt)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime
    real(r8), intent(in) :: dt

    integer i, j, k

    call perf_start('physics_update_dynamics')

    associate (mesh  => block%mesh               , &
               dudt  => block%aux%dudt_phys      , & ! in
               dvdt  => block%aux%dvdt_phys      , & ! in
               dptdt => block%aux%dptdt_phys     , & ! in
               dmg   => block%dstate(itime)%dmg  , & ! in
               u_lon => block%dstate(itime)%u_lon, & ! inout
               v_lat => block%dstate(itime)%v_lat, & ! inout
               pt    => block%dstate(itime)%pt   )   ! inout
    ! Update dynamics.
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          u_lon%d(i,j,k) = u_lon%d(i,j,k) + dt * 0.5_r8 * (dudt%d(i,j,k) + dudt%d(i+1,j,k))
        end do
      end do
    end do
    call fill_halo(u_lon)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          v_lat%d(i,j,k) = v_lat%d(i,j,k) + dt * 0.5_r8 * (dvdt%d(i,j,k) + dvdt%d(i,j+1,k))
        end do
      end do
    end do
    call fill_halo(v_lat)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          pt%d(i,j,k) = pt%d(i,j,k) + dt * dptdt%d(i,j,k) / dmg%d(i,j,k)
        end do
      end do
    end do
    call fill_halo(pt)
    end associate

    call perf_stop('physics_update_dynamics')

  end subroutine physics_update_dynamics

  subroutine physics_final()

    if (allocated(mesh)) deallocate(mesh)

    select case (physics_suite)
    case ('simple_physics')
      call simple_physics_final()
#ifdef HAS_CAM
    case ('cam')
      call cam_physics_final()
#endif
    case ('mars_nasa')
      call mars_nasa_final()
    end select

    if (allocated(physics_use_wet_tracers)) deallocate(physics_use_wet_tracers)

  end subroutine physics_final

  subroutine physics_add_output(tag)

    character(*), intent(in) :: tag

    if (physics_suite == 'N/A') return

    select case (physics_suite)
    case ('simple_physics')
      call simple_physics_add_output(tag, output_h0_dtype)
#ifdef HAS_CAM
    case ('cam')
      call cam_physics_add_output(tag, output_h0_dtype)
#endif
    case ('mars_nasa')
      call mars_nasa_add_output(tag, output_h0_dtype)
    end select

  end subroutine physics_add_output

  subroutine physics_output(tag, iblk)

    character(*), intent(in) :: tag
    integer, intent(in) :: iblk

    if (physics_suite == 'N/A') return

    select case (physics_suite)
    case ('simple_physics')
      call simple_physics_output(tag, iblk)
#ifdef HAS_CAM
    case ('cam')
      call cam_physics_output(tag, iblk)
#endif
    case ('mars_nasa')
      call mars_nasa_output(tag, iblk)
    end select

  end subroutine physics_output

end module physics_mod
