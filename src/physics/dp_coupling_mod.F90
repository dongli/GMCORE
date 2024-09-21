! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module dp_coupling_mod

  use const_mod
  use namelist_mod
  use block_mod
  use physics_types_mod
  use formula_mod
  use latlon_parallel_mod
  use tracer_mod
  use simple_physics_driver_mod   , only: simple_objects => objects
#ifdef HAS_CAM
  use cam_physics_driver_mod      , only: cam_physics_d2p, cam_physics_p2d, &
                                          cam_objects => objects
#endif
  use mars_nasa_physics_driver_mod, only: mars_nasa_d2p, &
                                          mars_nasa_p2d, &
                                          mars_nasa_objects => objects
  use filter_mod
  use perf_mod

  implicit none

  private

  public dp_coupling_d2p
  public dp_coupling_p2d

contains

  subroutine dp_coupling_d2p(block, itime)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime

    if (physics_suite == 'N/A') return

    call perf_start('dp_coupling_d2p')

    select case (physics_suite)
    case ('simple_physics')
      call common_d2p(block, itime, tracers(block%id), simple_objects(block%id)%state)
#ifdef HAS_CAM
    case ('cam')
      call common_d2p(block, itime, tracers(block%id), cam_objects(block%id)%state)
      call cam_physics_d2p()
#endif
    case ('mars_nasa')
      call common_d2p(block, itime, tracers(block%id), mars_nasa_objects(block%id)%state)
      call mars_nasa_d2p()
    end select

    call perf_stop('dp_coupling_d2p')

  contains

    subroutine common_d2p(block, itime, tracers, pstate)

      type(block_type), intent(in) :: block
      integer, intent(in) :: itime
      type(tracers_type), intent(in) :: tracers
      class(physics_state_type), intent(inout) :: pstate

      integer i, j, k, m, icol

      associate (mesh   => block%mesh         , &
                 dstate => block%dstate(itime), &
                 aux    => block%aux          )
      ! Full levels
      do k = mesh%full_kds, mesh%full_kde
        icol = 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            pstate%u     (icol,k) = aux    %u     %d(i,j,k)
            pstate%v     (icol,k) = aux    %v     %d(i,j,k)
            pstate%t     (icol,k) = aux    %t     %d(i,j,k)
            pstate%pt    (icol,k) = dstate %pt    %d(i,j,k)
            pstate%p     (icol,k) = dstate %ph    %d(i,j,k)
            pstate%p_lev (icol,k) = dstate %ph_lev%d(i,j,k)
            pstate%pk    (icol,k) = dstate %ph    %d(i,j,k)**rd_o_cpd / pk0
            pstate%pk_lev(icol,k) = dstate %ph_lev%d(i,j,k)**rd_o_cpd / pk0
            pstate%dp    (icol,k) = dstate %ph_lev%d(i,j,k+1) - dstate%ph_lev%d(i,j,k)
            pstate%dp_dry(icol,k) = dstate %dmg   %d(i,j,k)
            pstate%omg   (icol,k) = aux    %omg   %d(i,j,k)
            pstate%z     (icol,k) = dstate %gz    %d(i,j,k) / g
            pstate%dz    (icol,k) = (dstate%gz_lev%d(i,j,k+1) - dstate%gz_lev%d(i,j,k)) / g
            pstate%rho   (icol,k) = dry_air_density(pstate%t(icol,k), pstate%p(icol,k))
            icol = icol + 1
          end do
        end do
      end do
      ! Tracers
      do m = 1, ntracers
        if (physics_use_wet_tracers(m)) then
          do k = mesh%full_kds, mesh%full_kde
            icol = 1
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%full_ids, mesh%full_ide
                pstate%q(icol,k,m) = wet_mixing_ratio(tracers%q%d(i,j,k,m), tracers%qm%d(i,j,k))
                icol = icol + 1
              end do
            end do
          end do
        else
          do k = mesh%full_kds, mesh%full_kde
            icol = 1
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%full_ids, mesh%full_ide
                pstate%q(icol,k,m) = tracers%q%d(i,j,k,m)
                icol = icol + 1
              end do
            end do
          end do
        end if
      end do
      ! Half levels
      do k = mesh%half_kds, mesh%half_kde
        icol = 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            pstate%p_lev(icol,k) = dstate%ph_lev%d(i,j,k)
            pstate%z_lev(icol,k) = dstate%gz_lev%d(i,j,k) / g
            icol = icol + 1
          end do
        end do
      end do
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        icol = 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            pstate%dz_lev(icol,k) = pstate%z(icol,k) - pstate%z(icol,k-1)
            icol = icol + 1
          end do
        end do
      end do
      ! Surface
      icol = 1
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          pstate%ps     (icol) = dstate%phs%d(i,j)
          pstate%wsp_bot(icol) = sqrt(aux%u%d(i,j,mesh%full_kde)**2 + aux%v%d(i,j,mesh%full_kde)**2)
          icol = icol + 1
        end do
      end do
      end associate

    end subroutine common_d2p

  end subroutine dp_coupling_d2p

  subroutine dp_coupling_p2d(block, itime)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime

    integer iblk

    call perf_start('dp_coupling_p2d')

    iblk = block%id

    select case (physics_suite)
    case ('simple_physics')
      call common_p2d(block, itime, tracers(iblk), simple_objects(iblk)%state, simple_objects(iblk)%tend)
#ifdef HAS_CAM
    case ('cam')
      call cam_physics_p2d()
      call common_p2d(block, itime, tracers(iblk), cam_objects(iblk)%state, cam_objects(iblk)%tend)
#endif
    case ('mars_nasa')
      call mars_nasa_p2d()
      call common_p2d(block, itime, tracers(iblk), mars_nasa_objects(iblk)%state, mars_nasa_objects(iblk)%tend)
    end select

    call perf_stop('dp_coupling_p2d')

  contains

    subroutine common_p2d(block, itime, tracers, pstate, ptend)

      type(block_type), intent(inout) :: block
      integer, intent(in) :: itime
      type(tracers_type), intent(in) :: tracers
      class(physics_state_type), intent(in) :: pstate
      class(physics_tend_type), intent(inout) :: ptend

      integer i, j, k, m, icol

      associate (mesh   => block%mesh          , &
                 dstate => block%dstate(itime) , &
                 dudt   => block%aux%dudt_phys , & ! out
                 dvdt   => block%aux%dvdt_phys , & ! out
                 dptdt  => block%aux%dptdt_phys, & ! out
                 dqdt   => block%aux%dqdt_phys , & ! out
                 old_q  => tracers%q           , &
                 old_qm => tracers%qm          )
      if (ptend%updated_u .and. ptend%updated_v) then
        do k = mesh%full_kds, mesh%full_kde
          icol = 1
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              dudt%d(i,j,k) = ptend%dudt(icol,k)
              dvdt%d(i,j,k) = ptend%dvdt(icol,k)
              icol = icol + 1
            end do
          end do
        end do
        if (filter_ptend) then
          call filter_run(block%big_filter, dudt)
          call filter_run(block%big_filter, dvdt)
        end if
        call fill_halo(dudt, west_halo=.false., south_halo=.false., north_halo=.false.)
        call fill_halo(dvdt, west_halo=.false., east_halo=.false., south_halo=.false.)
      end if
      if (ptend%updated_t) then
        ! Convert temperature tendency to potential temperature tendency.
        do k = mesh%full_kds, mesh%full_kde
          icol = 1
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              dptdt%d(i,j,k) = dstate%dmg%d(i,j,k) / pstate%pk(icol,k) * ( &
                (1 + rv_o_rd * old_q%d(i,j,k,idx_qv)) * ptend%dtdt(icol,k) + &
                rv_o_rd * pstate%t(icol,k) * ptend%dqdt(icol,k,idx_qv) * (1 + old_qm%d(i,j,k)))
              icol = icol + 1
            end do
          end do
        end do
        if (filter_ptend) then
          call filter_run(block%big_filter, dptdt)
        end if
      else if (ptend%updated_pt) then
        do k = mesh%full_kds, mesh%full_kde
          icol = 1
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              dptdt%d(i,j,k) = dstate%dmg%d(i,j,k) * ptend%dptdt(icol,k)
              icol = icol + 1
            end do
          end do
        end do
        if (filter_ptend) then
          call filter_run(block%big_filter, dptdt)
        end if
      end if
      do m = 1, ntracers
        ! FIXME: Actually, most physics codes directly update q, so we do not need dqdt.
        if (ptend%updated_q(m)) then
          if (physics_use_wet_tracers(m)) then
            do k = mesh%full_kds, mesh%full_kde
              icol = 1
              do j = mesh%full_jds, mesh%full_jde
                do i = mesh%full_ids, mesh%full_ide
                  dqdt%d(i,j,k,m) = dstate%dmg%d(i,j,k) * ptend%dqdt(icol,k,m) * (1 + old_qm%d(i,j,k))
                  icol = icol + 1
                end do
              end do
            end do
          else
            do k = mesh%full_kds, mesh%full_kde
              icol = 1
              do j = mesh%full_jds, mesh%full_jde
                do i = mesh%full_ids, mesh%full_ide
                  dqdt%d(i,j,k,m) = dstate%dmg%d(i,j,k) * ptend%dqdt(icol,k,m)
                  icol = icol + 1
                end do
              end do
            end do
          end if
          if (filter_ptend) then
            call filter_run(block%big_filter, dqdt, m)
          end if
        end if
      end do
      end associate

    end subroutine common_p2d

  end subroutine dp_coupling_p2d

end module dp_coupling_mod
