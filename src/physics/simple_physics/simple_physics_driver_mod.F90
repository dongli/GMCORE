! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module simple_physics_driver_mod

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
  public simple_physics_p2d
  public simple_physics_add_output
  public simple_physics_output
  public objects

  real(r8) dt

contains

  subroutine simple_physics_init_stage1(namelist_path, dt_adv, dt_phys)

    character(*), intent(in) :: namelist_path
    real(r8), intent(in) :: dt_adv
    real(r8), intent(in) :: dt_phys

    call simple_physics_final()

    call tracer_add('moist', dt_adv, 'qv', 'Water vapor', 'kg kg-1')

    dt = dt_phys

  end subroutine simple_physics_init_stage1

  subroutine simple_physics_init_stage2(mesh)

    type(physics_mesh_type), intent(in), target :: mesh(:)

    call simple_physics_objects_init(mesh)

  end subroutine simple_physics_init_stage2

  subroutine simple_physics_final()

    call simple_physics_objects_final()

  end subroutine simple_physics_final

  subroutine simple_physics_run()

    integer iblk

    do iblk = 1, size(objects)
      associate (mesh  => objects(iblk)%mesh , &
                 state => objects(iblk)%state, &
                 tend  => objects(iblk)%tend )
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
      end associate
    end do

  end subroutine simple_physics_run

  subroutine simple_physics_p2d()

    integer iblk, icol, k

    do iblk = 1, size(objects)
      associate (mesh  => objects(iblk)%mesh , &
                 state => objects(iblk)%state, &
                 tend  => objects(iblk)%tend )
      do k = 1, mesh%nlev
        do icol = 1, mesh%ncol
          tend%dptdt(icol,k) = (p0 / state%p(icol,k))**rd_o_cpd * ( &
            (1 + rv_o_rd * state%qv(icol,k)) * tend%dtdt(icol,k) +  &
            rv_o_rd * state%t(icol,k) * tend%dqvdt(icol,k))
          ! Convert to dry mixing ratio tendency.
          tend%dqvdt(icol,k) = tend%dqvdt(icol,k) / (1 - state%qv(icol,k))**2
        end do
      end do
      end associate
    end do

  end subroutine simple_physics_p2d

end module simple_physics_driver_mod
