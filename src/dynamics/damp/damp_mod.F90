! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module damp_mod

  use const_mod
  use namelist_mod
  use block_mod
  use tracer_mod
  use div_damp_mod
  use vor_damp_mod
  use smag_damp_mod
  use laplace_damp_mod
  use latlon_parallel_mod

  implicit none

  private

  public damp_init
  public damp_final
  public damp_run
  public damp_update_uv
  public damp_update_w
  public damp_update_pt
  public damp_update_q

contains

  subroutine damp_init()

    call div_damp_init()
    call vor_damp_init()
    call smag_damp_init()
    call laplace_damp_init()

  end subroutine damp_init

  subroutine damp_final()

    call div_damp_final()
    call vor_damp_final()
    call smag_damp_final()
    call laplace_damp_final()

  end subroutine damp_final

  subroutine damp_run(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    integer j, m
    real(r8) c
    associate (mesh => block%mesh, u_lon => dstate%u_lon, v_lat => dstate%v_lat)
    ! This nudging of polar v helps to keep the flow neat around the poles.
    ! NOTE: DO NOT REMOVE IT!
    c = 0.8_r8
    do j = mesh%half_jms, mesh%half_jme
      if (mesh%is_south_pole(j)) then
        v_lat%d(:,j  ,:) = c * v_lat%d(:,j  ,:) + (1 - c) * v_lat%d(:,j+1,:)
      else if (mesh%is_north_pole(j+1)) then
        v_lat%d(:,j  ,:) = c * v_lat%d(:,j  ,:) + (1 - c) * v_lat%d(:,j-1,:)
      end if
    end do
    end associate

    if (use_laplace_damp) then
      call laplace_damp_run(block, dstate%u_lon, 2, 500.0_r8, block%aux%dudt_damp)
      call laplace_damp_run(block, dstate%v_lat, 2, 500.0_r8, block%aux%dvdt_damp)
      if (nonhydrostatic) call laplace_damp_run(block, dstate%w_lev, 2, 500.0_r8, block%aux%dwdt_damp)
      call laplace_damp_run(block, dstate%pt, 2, 1500.0_r8, block%aux%dptdt_damp)
      call block%aux%dptdt_damp%mul(dstate%dmg)
      do m = 1, ntracers
        call laplace_damp_run(block, tracers(block%id)%q, m, 2, 1500.0_r8, block%aux%dqdt_damp)
        call block%aux%dqdt_damp%mul(m, dstate%dmg)
      end do
    end if
    if (use_vor_damp) then
      call vor_damp_run(block, dstate, dt)
    end if
    if (use_div_damp) then
      call div_damp_run(block, dstate, dt)
    end if
    if (use_smag_damp) then
      call smag_damp_run(block, dstate, dt)
    end if

  end subroutine damp_run

  subroutine damp_update_uv(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    integer i, j, k

    associate (mesh  => block%mesh         , &
               dmg   => dstate%dmg         , & ! in
               dudt  => block%aux%dudt_damp, & ! in
               dvdt  => block%aux%dvdt_damp, & ! in
               u_lon => dstate%u_lon       , & ! out
               v_lat => dstate%v_lat       )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          u_lon%d(i,j,k) = u_lon%d(i,j,k) + dt * dudt%d(i,j,k)
        end do
      end do
    end do
    call fill_halo(u_lon)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          v_lat%d(i,j,k) = v_lat%d(i,j,k) + dt * dvdt%d(i,j,k)
        end do
      end do
    end do
    call fill_halo(v_lat)
    end associate

  end subroutine damp_update_uv

  subroutine damp_update_w(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    integer i, j, k

    associate (mesh  => block%mesh         , &
               dmg   => dstate%dmg         , & ! in
               dwdt  => block%aux%dwdt_damp, & ! in
               w_lev => dstate%w_lev       )   ! out
    do k = mesh%half_kds + 1, mesh%half_kde - 1
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          w_lev%d(i,j,k) = w_lev%d(i,j,k) + dt * dwdt%d(i,j,k)
        end do
      end do
    end do
    call fill_halo(w_lev)
    end associate

  end subroutine damp_update_w

  subroutine damp_update_pt(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    integer i, j, k

    associate (mesh  => block%mesh          , &
               dmg   => dstate%dmg          , & ! in
               dptdt => block%aux%dptdt_damp, & ! in
               pt    => dstate%pt           )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          pt%d(i,j,k) = pt%d(i,j,k) + dt * dptdt%d(i,j,k) / dmg%d(i,j,k)
        end do
      end do
    end do
    call fill_halo(pt)
    end associate

  end subroutine damp_update_pt

  subroutine damp_update_q(block, dstate, dt, m)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt
    integer, intent(in) :: m

    integer i, j, k

    associate (mesh => block%mesh         , &
               dmg  => dstate%dmg         , & ! in
               dqdt => block%aux%dqdt_damp, & ! in
               q    => tracers(block%id)%q)   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          q%d(i,j,k,m) = q%d(i,j,k,m) + dt * dqdt%d(i,j,k,m) / dmg%d(i,j,k)
        end do
      end do
    end do
    call fill_halo(q, m)
    end associate

  end subroutine damp_update_q

end module damp_mod
