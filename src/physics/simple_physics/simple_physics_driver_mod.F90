! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module simple_physics_driver_mod

  use formula_mod
  use namelist_mod
  use tracer_mod
  use simple_physics_types_mod
  use simple_physics_objects_mod
  use simple_physics_output_mod

  implicit none

  private

  public simple_physics_init_stage1
  public simple_physics_init_stage2
  public simple_physics_final
  public simple_physics_run
  public simple_physics_add_output
  public simple_physics_output
  public objects

  real(r8) dt

contains

  subroutine simple_physics_init_stage1(namelist_path, dt_adv, dt_phys)

    character(*), intent(in) :: namelist_path
    real(r8), intent(in) :: dt_adv
    real(r8), intent(in) :: dt_phys

    integer ignore

    call simple_physics_final()

    select case (physics_suite)
    case ('simple_physics:v6')
      if (idx_qv == 0) call tracer_add('moist', dt_adv, 'qv', 'Water vapor', 'kg kg-1')
      allocate(physics_use_wet_tracers(ntracers))
      physics_use_wet_tracers = .true.
    case ('simple_physics:kessler')
      if (idx_qv == 0) call tracer_add('moist', dt_adv, 'qv', 'Water vapor', 'kg kg-1')
      if (idx_qc == 0) call tracer_add('moist', dt_adv, 'qc', 'Cloud water', 'kg kg-1')
      if (idx_qr == 0) call tracer_add('moist', dt_adv, 'qr', 'Rain water' , 'kg kg-1')
      allocate(physics_use_wet_tracers(ntracers))
      physics_use_wet_tracers = .false.
    end select

    dt = dt_phys

  end subroutine simple_physics_init_stage1

  subroutine simple_physics_init_stage2(mesh)

    type(physics_mesh_type), intent(in), target :: mesh(:)

    call simple_physics_objects_init(mesh)

  end subroutine simple_physics_init_stage2

  subroutine simple_physics_final()

    if (allocated(physics_use_wet_tracers)) deallocate(physics_use_wet_tracers)

    call simple_physics_objects_final()

  end subroutine simple_physics_final

  subroutine simple_physics_run()

    integer iblk, icol, k, m

    do iblk = 1, size(objects)
      associate (mesh  => objects(iblk)%mesh , &
                 state => objects(iblk)%state, &
                 tend  => objects(iblk)%tend )
      select case (physics_suite)
      case ('simple_physics:v6')
        call simple_physics(     &
          mesh%ncol            , &
          mesh%nlev            , &
          dt                   , &
          mesh%lat             , &
          state%t              , &
          state%qv             , &
          state%u              , &
          state%v              , &
          state%p              , &
          state%p_lev          , &
          state%dp             , &
          state%ps             , &
          tend%dudt            , &
          tend%dvdt            , &
          tend%dtdt            , &
          tend%dqdt(:,:,idx_qv), &
          state%precl          , &
          0                    , & ! test
          .true.               , & ! RJ2012_precip
          .false.                & ! TC_PBL_mod
        )
        tend%updated_u         = .true.
        tend%updated_v         = .true.
        tend%updated_t         = .true.
        tend%updated_q(idx_qv) = .true.
      case ('simple_physics:kessler')
        do icol = 1, mesh%ncol
          call kessler(          &
            mesh%nlev          , &
            state%pt   (icol,:), &
            state%qv   (icol,:), &
            state%qc   (icol,:), &
            state%qr   (icol,:), &
            state%rhod (icol,:), &
            state%pk   (icol,:), &
            dt                 , &
            state%z    (icol,:), &
            state%precl(icol  )  &
          )
        end do
        do m = 1, ntracers
          do k = 1, mesh%nlev
            do icol = 1, mesh%ncol
              tend%dqdt(icol,k,m) = (state%q(icol,k,m) - state%q_old(icol,k,m)) / dt
            end do
          end do
        end do
        tend%updated_q = .true.
        do k = 1, mesh%nlev
          do icol = 1, mesh%ncol
            tend%dptdt(icol,k) = (1 + rv_o_rd * state%qv_old(icol,k)) * (state%pt(icol,k) - state%pt_old(icol,k)) / dt + &
                                 rv_o_rd * state%pt_old(icol,k) * tend%dqvdt(icol,k)
          end do
        end do
        tend%updated_pt = .true.
      end select
      end associate
    end do

  end subroutine simple_physics_run

end module simple_physics_driver_mod
