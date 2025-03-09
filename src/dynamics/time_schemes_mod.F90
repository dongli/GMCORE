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
!   This module implements time integration schemes, including PC2 and WRFRK3,
!   and with forward-backward technique.
!
! Authors:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
!   - Jianghao Li
! ==============================================================================

module time_schemes_mod

  use flogger
  use const_mod
  use namelist_mod
  use block_mod
  use operators_mod
  use latlon_parallel_mod
  use process_mod, only: proc
  use physics_mod
  use filter_mod
  use math_mod
  use perf_mod

  implicit none

  private

  public time_scheme_init
  public time_scheme_final
  public time_integrator
  public update_state
  public space_operators_interface

  interface
    subroutine space_operators_interface(block, old_dstate, star_dstate, new_dstate, dtend, dt, pass, substep)
      import block_type, dstate_type, dtend_type, r8
      type(block_type ), intent(inout) :: block
      type(dstate_type), intent(in   ) :: old_dstate
      type(dstate_type), intent(inout) :: star_dstate
      type(dstate_type), intent(inout) :: new_dstate
      type(dtend_type ), intent(inout) :: dtend
      real(r8), intent(in) :: dt
      integer, intent(in) :: pass
      integer, intent(in) :: substep
    end subroutine space_operators_interface

    subroutine step_interface(space_operators, block, old_dstate, star_dstate, new_dstate, dtend, dt, substep)
      import space_operators_interface, block_type, dstate_type, dtend_type, r8
      procedure(space_operators_interface) space_operators
      type(block_type ), intent(inout) :: block
      type(dstate_type), intent(in   ) :: old_dstate
      type(dstate_type), intent(inout) :: star_dstate
      type(dstate_type), intent(inout) :: new_dstate
      type(dtend_type ), intent(inout) :: dtend
      real(r8), intent(in) :: dt
      integer, intent(in) :: substep
    end subroutine step_interface

    subroutine time_integrator_interface(space_operators, block, old, new, dt)
      import block_type, dtend_type, dstate_type, space_operators_interface, r8
      procedure(space_operators_interface) space_operators
      type(block_type), intent(inout) :: block
      integer, intent(in) :: old
      integer, intent(in) :: new
      real(r8), intent(in) :: dt
    end subroutine time_integrator_interface
  end interface

  procedure(step_interface), pointer :: step
  procedure(time_integrator_interface), pointer :: time_integrator

contains

  subroutine time_scheme_init()

    call time_scheme_final()

    select case (time_scheme)
    case ('euler')
      time_integrator => euler
      total_substeps = 1
    case ('rk2')
      time_integrator => rk2
      total_substeps = 2
    case ('pc2')
      time_integrator => pc2
      total_substeps = 3
    case ('wrfrk3')
      time_integrator => wrfrk3
      total_substeps = 3
    case default
      time_integrator => pc2
      total_substeps = 3
    end select

    step => step_forward_backward

  end subroutine time_scheme_init

  subroutine time_scheme_final()

  end subroutine time_scheme_final

  subroutine step_all(space_operators, block, old_dstate, star_dstate, new_dstate, dtend, dt, substep)

    procedure(space_operators_interface) space_operators
    type(block_type ), intent(inout) :: block
    type(dstate_type), intent(in   ) :: old_dstate
    type(dstate_type), intent(inout) :: star_dstate
    type(dstate_type), intent(inout) :: new_dstate
    type(dtend_type ), intent(inout) :: dtend
    real(r8), intent(in) :: dt
    integer, intent(in) :: substep

    call space_operators(block, old_dstate, star_dstate, new_dstate, dtend, dt, all_pass, substep)
    call update_state(block, dtend, old_dstate, new_dstate, dt, all_pass, substep)

  end subroutine step_all

  subroutine step_forward_backward(space_operators, block, old_dstate, star_dstate, new_dstate, dtend, dt, substep)

    procedure(space_operators_interface) space_operators
    type(block_type ), intent(inout) :: block
    type(dstate_type), intent(in   ) :: old_dstate
    type(dstate_type), intent(inout) :: star_dstate
    type(dstate_type), intent(inout) :: new_dstate
    type(dtend_type ), intent(inout) :: dtend
    real(r8), intent(in) :: dt
    integer, intent(in) :: substep

    call space_operators(block, old_dstate, star_dstate, new_dstate, dtend, dt, forward_pass, substep)
    call update_state(block, dtend, old_dstate, new_dstate, dt, forward_pass, substep)
    call space_operators(block, old_dstate, star_dstate, new_dstate, dtend, dt, backward_pass, substep)
    call update_state(block, dtend, old_dstate, new_dstate, dt, backward_pass, substep)

  end subroutine step_forward_backward

  subroutine update_state(block, dtend, old_dstate, new_dstate, dt, pass, substep)

    type(block_type ), intent(inout) :: block
    type(dtend_type ), intent(inout) :: dtend
    type(dstate_type), intent(in   ) :: old_dstate
    type(dstate_type), intent(inout) :: new_dstate
    real(r8), intent(in) :: dt
    integer, intent(in) :: pass
    integer, intent(in) :: substep

    integer i, j, k
    real(r8) c, tmp(block%mesh%full_ids-1:block%mesh%full_ide+1)

    associate (mesh   => block%mesh  , &
               dmgsdt => dtend%dmgsdt, &
               dgzdt  => dtend%dgzdt , &
               dptdt  => dtend%dptdt , &
               dudt   => dtend%dudt  , &
               dvdt   => dtend%dvdt  )
    if (baroclinic) then
      if (dtend%update_mgs) then
        ! ----------------------------------------------------------------------
        call filter_run(block%big_filter, dmgsdt)
        ! ----------------------------------------------------------------------
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            new_dstate%mgs%d(i,j) = old_dstate%mgs%d(i,j) + dt * dmgsdt%d(i,j)
          end do
        end do
        call fill_halo(new_dstate%mgs)
        call calc_mg (block, new_dstate)
        call calc_dmg(block, new_dstate)
        call calc_ph (block, new_dstate)
      end if

      if (dtend%update_pt) then
        if (.not. dtend%update_mgs .and. proc%is_root()) call log_error('Mass is not updated or copied!')
        ! ----------------------------------------------------------------------
        call filter_run(block%big_filter, dptdt)
        ! ----------------------------------------------------------------------
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              new_dstate%pt%d(i,j,k) = (old_dstate%pt%d(i,j,k) * old_dstate%dmg%d(i,j,k) + dt * dptdt%d(i,j,k)) / new_dstate%dmg%d(i,j,k)
            end do
          end do
        end do
        call fill_halo(new_dstate%pt, async=.true.)
      end if
    else
      if (dtend%update_gz) then
        ! ----------------------------------------------------------------------
        call filter_run(block%big_filter, dgzdt)
        ! ----------------------------------------------------------------------
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              new_dstate%gz%d(i,j,k) = old_dstate%gz%d(i,j,k) + dt * dgzdt%d(i,j,k)
            end do
          end do
        end do
        call fill_halo(new_dstate%gz)
        call calc_dmg(block, new_dstate)
      end if
    end if

    if (dtend%update_u .and. dtend%update_v) then
      ! ------------------------------------------------------------------------
      call filter_run(block%big_filter, dudt)
      call filter_run(block%big_filter, dvdt)
      ! ------------------------------------------------------------------------
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            new_dstate%u_lon%d(i,j,k) = old_dstate%u_lon%d(i,j,k) + dt * dudt%d(i,j,k)
          end do
        end do
      end do
      call fill_halo(new_dstate%u_lon, async=.true.)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            new_dstate%v_lat%d(i,j,k) = old_dstate%v_lat%d(i,j,k) + dt * dvdt%d(i,j,k)
          end do
        end do
      end do
      ! ------------------------------------------------------------------------
      ! This nudging of polar v helps to keep the flow neat around the poles.
      ! NOTE: DO NOT REMOVE IT!
      c = 0.8_r8
      if (mesh%has_south_pole()) then
        j = mesh%half_jds
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            new_dstate%v_lat%d(i,j,k) = c * new_dstate%v_lat%d(i,j,k) + (1 - c) * new_dstate%v_lat%d(i,j+1,k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%half_jde
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            new_dstate%v_lat%d(i,j,k) = c * new_dstate%v_lat%d(i,j,k) + (1 - c) * new_dstate%v_lat%d(i,j-1,k)
          end do
        end do
      end if
      ! ------------------------------------------------------------------------
      call fill_halo(new_dstate%v_lat, south_halo=.false., north_halo=.false.)
      do j = mesh%half_jds, mesh%half_jde
        c = exp_two_values(1.0_r8, 0.0_r8, global_mesh%half_lat_deg(global_mesh%half_nlat), 80.0_r8, abs(mesh%half_lat_deg(j)))
        do k = mesh%full_kds, mesh%full_kde
          tmp = new_dstate%v_lat%d(mesh%full_ids-1:mesh%full_ide+1,j,k)
          do i = mesh%full_ids, mesh%full_ide
            new_dstate%v_lat%d(i,j,k) = (1 - 0.5_r8 * c) * new_dstate%v_lat%d(i,j,k) + 0.25_r8 * c * (tmp(i-1) + tmp(i+1))
          end do
        end do
      end do
      ! ------------------------------------------------------------------------
      call fill_halo(new_dstate%v_lat, async=.true.)
    end if
    end associate

    call physics_update_after_rk_substep(block, new_dstate, dt, pass)

  end subroutine update_state

  subroutine rk2(space_operators, block, old, new, dt)

    procedure(space_operators_interface) space_operators
    type(block_type ), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(r8), intent(in) :: dt

    associate (dstate => block%dstate, dtend => block%dtend)
    call step(space_operators, block, dstate(old), dstate(old), dstate(new), dtend, dt / 2.0_r8, 1)
    call step(space_operators, block, dstate(old), dstate(new), dstate(new), dtend, dt         , 2)
    end associate

  end subroutine rk2

  subroutine pc2(space_operators, block, old, new, dt)

    procedure(space_operators_interface) space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(r8), intent(in) :: dt

    associate (dstate => block%dstate, dtend => block%dtend)
    call step(space_operators, block, dstate(old), dstate(old), dstate(new), dtend, dt / 2.0_r8, 1)
    call step(space_operators, block, dstate(old), dstate(new), dstate(3  ), dtend, dt / 2.0_r8, 2)
    call step(space_operators, block, dstate(old), dstate(3  ), dstate(new), dtend, dt         , 3)
    end associate

  end subroutine pc2

  subroutine wrfrk3(space_operators, block, old, new, dt)

    procedure(space_operators_interface) space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(r8), intent(in) :: dt

    associate (dstate => block%dstate, dtend => block%dtend)
    call step(space_operators, block, dstate(old), dstate(old), dstate(new), dtend, dt / 3.0_r8, 1)
    call step(space_operators, block, dstate(old), dstate(new), dstate(3  ), dtend, dt / 2.0_r8, 2)
    call step(space_operators, block, dstate(old), dstate(3  ), dstate(new), dtend, dt         , 3)
    end associate

  end subroutine wrfrk3

  subroutine euler(space_operators, block, old, new, dt)

    procedure(space_operators_interface) space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(r8), intent(in) :: dt

    associate (dstate => block%dstate, dtend => block%dtend)
    call step(space_operators, block, dstate(old), dstate(old), dstate(new), dtend, dt, 1)
    end associate

  end subroutine euler

end module time_schemes_mod
