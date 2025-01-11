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
  use gomars_v1_driver_mod
  use gomars_v2_driver_mod
  use perf_mod

  implicit none

  private

  public physics_init_stage1
  public physics_init_stage2
  public physics_init_stage3
  public physics_run
  public physics_update_after_dynamics
  public physics_update_after_advection
  public physics_update_after_physics
  public physics_update_after_rk_substep
  public physics_final
  public physics_add_output
  public physics_output

  type(physics_mesh_type), allocatable :: mesh(:)

contains

  subroutine physics_init_stage1(namelist_path)

    character(*), intent(in) :: namelist_path

    select case (physics_suite)
    case ('simple_physics:v6', 'simple_physics:kessler')
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
      call time_add_alert('phys', seconds=dt_phys/time_scale)
    else
      return
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
                           blocks(iblk)%mesh%full_lev , &
                           blocks(iblk)%mesh%half_lev , &
                           blocks(iblk)%mesh%full_dlev, &
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
    case ('simple_physics:v6', 'simple_physics:kessler')
      call simple_physics_init_stage2(mesh)
    case ('cam')
#ifdef HAS_CAM
      call cam_physics_init_stage2(namelist_path, mesh, dt_adv, dt_phys)
#else
      if (proc%is_root()) call log_error('CAM physics is not compiled!')
#endif
    case ('gomars_v1')
      call gomars_v1_init_stage2(namelist_path, mesh, dt_adv, dt_phys, input_ngroups, gmcore_root)
    case ('gomars_v2')
      call gomars_v2_init_stage2(namelist_path, mesh, dt_adv, dt_phys, input_ngroups, gmcore_root)
    end select

  end subroutine physics_init_stage2

  subroutine physics_init_stage3()

    integer iblk

    do iblk = 1, size(blocks)
      call dp_coupling_d2p(blocks(iblk), 1)
    end do

    select case (physics_suite)
    case ('gomars_v1')
      call gomars_v1_init_stage3()
    case ('gomars_v2')
      call gomars_v2_init_stage3()
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
    case ('simple_physics:v6', 'simple_physics:kessler')
      call simple_physics_run()
#ifdef HAS_CAM
    case ('cam')
      call cam_physics_run_stage1()
      call cam_physics_sfc_flux()
      call cam_physics_run_stage2()
#endif
    case ('gomars_v1')
      call gomars_v1_run(curr_time)
    case ('gomars_v2')
      call gomars_v2_run(curr_time)
    end select

    call dp_coupling_p2d(block, itime)

    call perf_stop('physics_run')

  end subroutine physics_run

  subroutine physics_update_uv(block, dstate, dt)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    integer i, j, k

    associate (mesh      => block%mesh         , &
               dudt_phys => block%aux%dudt_phys, & ! in
               dvdt_phys => block%aux%dvdt_phys, & ! in
               dudt_damp => block%aux%dudt_damp, & ! in
               dvdt_damp => block%aux%dvdt_damp, & ! in
               u_lon     => dstate%u_lon       , & ! inout
               v_lat     => dstate%v_lat       )   ! inout
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          u_lon%d(i,j,k) = u_lon%d(i,j,k) + dt * 0.5_r8 * (dudt_phys%d(i,j,k) + dudt_phys%d(i+1,j,k))
        end do
      end do
    end do
    call fill_halo(u_lon)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          v_lat%d(i,j,k) = v_lat%d(i,j,k) + dt * 0.5_r8 * (dvdt_phys%d(i,j,k) + dvdt_phys%d(i,j+1,k))
        end do
      end do
    end do
    call fill_halo(v_lat)
    end associate

  end subroutine physics_update_uv

  subroutine physics_update_mgs(block, dstate, dt)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    integer i, j

    associate (mesh  => block%mesh          , &
               dpsdt => block%aux%dpsdt_phys, & ! in
               mgs   => dstate%mgs          )   ! inout
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        mgs%d(i,j) = mgs%d(i,j) + dt * dpsdt%d(i,j)
      end do
    end do
    call fill_halo(mgs)
    end associate

  end subroutine physics_update_mgs

  subroutine physics_update_pt(block, dstate, dt)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    integer i, j, k

    associate (mesh       => block%mesh          , &
               dptdt_phys => block%aux%dptdt_phys, & ! in
               dmg        => dstate%dmg          , & ! in
               pt         => dstate%pt           )   ! inout
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          pt%d(i,j,k) = pt%d(i,j,k) + dt * dptdt_phys%d(i,j,k) / dmg%d(i,j,k)
        end do
      end do
    end do
    call fill_halo(pt)
    end associate

  end subroutine physics_update_pt

  subroutine physics_update_q(block, dstate, dt)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(in) :: dstate
    real(r8), intent(in) :: dt

    integer i, j, k, m

    associate (mesh      => block%mesh         , &
               dqdt_phys => block%aux%dqdt_phys, & ! in
               dmg       => dstate%dmg         , & ! in
               q         => tracers(block%id)%q)   ! inout
    do m = 1, ntracers
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            q%d(i,j,k,m) = q%d(i,j,k,m) + dt * dqdt_phys%d(i,j,k,m) / dmg%d(i,j,k)
          end do
        end do
      end do
      call fill_halo(q, m)
    end do
    call tracer_calc_qm(block)
    end associate

  end subroutine physics_update_q

  subroutine physics_update_after_dynamics(block, itime, dt)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime
    real(r8), intent(in) :: dt

    call perf_start('physics_update_dynamics')

    select case (pdc_type)
    case (1)
      call physics_update_uv (block, block%dstate(itime), dt)
      call physics_update_mgs(block, block%dstate(itime), dt)
      call physics_update_pt (block, block%dstate(itime), dt)
      call physics_update_q  (block, block%dstate(itime), dt)
    case (14)
      call physics_update_uv (block, block%dstate(itime), dt)
      call physics_update_mgs(block, block%dstate(itime), dt)
      call physics_update_q  (block, block%dstate(itime), dt)
    end select

    call perf_stop('physics_update_after_dynamics')

  end subroutine physics_update_after_dynamics

  subroutine physics_update_after_advection(block, itime, dt)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime
    real(r8), intent(in) :: dt

    call perf_start('physics_update_after_advection')

    select case (pdc_type)
    case (2)
      call physics_update_uv (block, block%dstate(itime), dt)
      call physics_update_mgs(block, block%dstate(itime), dt)
      call physics_update_pt (block, block%dstate(itime), dt)
      call physics_update_q  (block, block%dstate(itime), dt)
    case (24)
      call physics_update_uv (block, block%dstate(itime), dt)
      call physics_update_mgs(block, block%dstate(itime), dt)
      call physics_update_q  (block, block%dstate(itime), dt)
    case (5)
      call physics_update_q  (block, block%dstate(itime), dt)
    end select

    call perf_stop('physics_update_after_advection')

  end subroutine physics_update_after_advection

  subroutine physics_update_after_physics(block, itime, dt)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime
    real(r8), intent(in) :: dt

    call perf_start('physics_update_after_physics')

    select case (pdc_type)
    case (3)
      call physics_update_uv (block, block%dstate(itime), dt)
      call physics_update_mgs(block, block%dstate(itime), dt)
      call physics_update_pt (block, block%dstate(itime), dt)
      call physics_update_q  (block, block%dstate(itime), dt)
    case (34)
      call physics_update_uv (block, block%dstate(itime), dt)
      call physics_update_mgs(block, block%dstate(itime), dt)
      call physics_update_q  (block, block%dstate(itime), dt)
    end select

    call perf_stop('physics_update_after_physics')

  end subroutine physics_update_after_physics

  subroutine physics_update_after_rk_substep(block, dstate, dt, pass)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt
    integer, intent(in) :: pass

    call perf_start('physics_update_after_rk_substep')

    select case (pdc_type)
    case (4)
      select case (pass)
      case (backward_pass)
        call physics_update_uv (block, dstate, dt)
      case (forward_pass)
        call physics_update_mgs(block, dstate, dt)
        call physics_update_pt (block, dstate, dt)
        call physics_update_q  (block, dstate, dt)
      end select
    case (14, 24, 34)
      select case (pass)
      case (forward_pass)
        call physics_update_pt (block, dstate, dt)
      end select
    case (5)
      select case (pass)
      case (backward_pass)
        call physics_update_uv (block, dstate, dt)
      case (forward_pass)
        call physics_update_pt (block, dstate, dt)
      end select
    end select

    call perf_stop('physics_update_after_rk_substep')

  end subroutine physics_update_after_rk_substep

  subroutine physics_final()

    if (allocated(mesh)) deallocate(mesh)

    select case (physics_suite)
    case ('simple_physics:v6', 'simple_physics:kessler')
      call simple_physics_final()
#ifdef HAS_CAM
    case ('cam')
      call cam_physics_final()
#endif
    case ('gomars_v1')
      call gomars_v1_final()
    case ('gomars_v2')
      call gomars_v2_final()
    end select

    if (allocated(physics_use_wet_tracers)) deallocate(physics_use_wet_tracers)

  end subroutine physics_final

  subroutine physics_add_output(tag)

    character(*), intent(in) :: tag

    if (physics_suite == 'N/A') return

    select case (physics_suite)
    case ('simple_physics:v6', 'simple_physics:kessler')
      call simple_physics_add_output(tag, output_h0_dtype)
#ifdef HAS_CAM
    case ('cam')
      call cam_physics_add_output(tag, output_h0_dtype)
#endif
    case ('gomars_v2')
      call gomars_v2_add_output(tag, output_h0_dtype)
    end select

  end subroutine physics_add_output

  subroutine physics_output(tag, iblk)

    character(*), intent(in) :: tag
    integer, intent(in) :: iblk

    if (physics_suite == 'N/A') return

    select case (physics_suite)
    case ('simple_physics:v6', 'simple_physics:kessler')
      call simple_physics_output(tag, iblk)
#ifdef HAS_CAM
    case ('cam')
      call cam_physics_output(tag, iblk)
#endif
    case ('gomars_v2')
      call gomars_v2_output(tag, iblk)
    end select

  end subroutine physics_output

end module physics_mod
