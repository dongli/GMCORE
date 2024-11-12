! ==============================================================================
! This file is part of GoMars since 2023.
!
! GoMars is a Martian general circulation model developed in Institute of
! Atmospheric Physics (IAP), Chinese Academy of Sciences (CAS).
!
! GMCORE is a dynamical core for atmospheric model used in GoMars.
!
! GoMars and GMCORE are distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module gomars_v1_driver_mod

  use datetime
  use const_mod
  use namelist_mod, only: restart
  use vert_coord_mod
  use gomars_v1_namelist_mod
  use gomars_v1_objects_mod
  use gomars_v1_tracers_mod
  use gomars_v1_orbit_mod
  use gomars_v1_rad_mod
  use gomars_v1_pbl_mod
  use gomars_v1_lsm_mod
  use gomars_v1_mp_mod
  use gomars_v1_damp_mod

  implicit none

  private

  public gomars_v1_init_stage2
  public gomars_v1_init_stage3
  public gomars_v1_run
  public gomars_v1_final
  public gomars_v1_d2p
  public gomars_v1_p2d
  public objects

contains

  ! ============================================================================
  ! Description:
  !
  !   This subroutine initialize vertical coordinate which may depend on the
  !   topography. The tracers are also allocated, so users should already add
  !   the necessary tracer species.
  ! ============================================================================

  subroutine gomars_v1_init_stage2(namelist_path, mesh, dt_adv, dt_phys, input_ngroup, model_root)

    character(*), intent(in) :: namelist_path
    type(physics_mesh_type), intent(in), target :: mesh(:)
    real(r8), intent(in) :: dt_adv
    real(r8), intent(in) :: dt_phys
    integer , intent(in) :: input_ngroup
    character(*), intent(in), optional :: model_root

    dt = dt_phys

    call gomars_v1_tracers_init(dt_adv)
    call gomars_v1_objects_init(mesh)
    call gomars_v1_orbit_init()
    call gomars_v1_rad_init()
    call gomars_v1_pbl_init()
    call gomars_v1_damp_init()

  end subroutine gomars_v1_init_stage2

  ! ============================================================================
  ! Description:
  !
  !   This subroutine runs some initializations that need to be after initial
  !   conditions.
  ! ============================================================================

  subroutine gomars_v1_init_stage3()

    integer iblk, icol, k, l, is, ig, n

    ! Check model top pressure. It must be lower than ptop in nasa_rad_mod.

    do iblk = 1, size(objects)
      associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state)
      n = 2 * mesh%nlev + 3
      ! Set reference pressure and interpolation coeficients (from emiss).
      state%pl(2) = ptrop * 0.5_r8
      do k = 3, n, 2
        state%pl(k) = vert_coord_calc_mg_lev((k - 1) / 2, psl)
      end do
      do k = 4, n - 1, 2
        state%pl(k) = vert_coord_calc_mg((k - 2) / 2, psl)
      end do
      do k = 3, n - 1
        state%aadj(k) = log(state%pl(k+1) / state%pl(k))
      end do
      do k = 3, n - 2, 2
        l = (k - 1) / 2
        state%badj(k) = mesh%lev(l) / state%pl(k+1) - mesh%ilev(l) / state%pl(k)
      end do
      do k = 4, n - 1, 2
        l = (k - 2) / 2
        state%badj(k) = mesh%ilev(l+1) / state%pl(k+1) - mesh%lev(l) / state%pl(k)
      end do
      ! Set some initial variables (from init1).
      if (.not. restart) then
        ! FIXME: Here is for cold run.
        do icol = 1, mesh%ncol
          state%tstrat(icol) = state%t(icol,mesh%nlev)
          state%co2ice(icol) = 0
          state%gt    (icol) = state%t(icol,mesh%nlev)
        end do
      end if
      call ini_optdst(qextv, qscatv, gv, qexti, qscati, gi, &
                      state%qxvdst, state%qsvdst, state%gvdst, &
                      state%qxidst, state%qsidst, state%gidst, &
                      state%qextrefdst)
      call ini_optcld(state%qxvcld, state%qsvcld, state%gvcld, &
                      state%qxicld, state%qsicld, state%gicld, &
                      state%qextrefcld, state%taurefcld)
      ! firstcomp3:
      if (.not. restart) then
        do icol = 1, mesh%ncol
          state%dnirflux(icol) = 1
          state%dndiffv (icol) = 0
          do is = 1, l_nspectv
            do ig = 1, l_ngauss
              state%detau(icol,is,ig) = 0.1_r8
            end do
          end do
        end do
      end if
      end associate
    end do

  end subroutine gomars_v1_init_stage3

  subroutine gomars_v1_run(time)

    type(datetime_type), intent(in) :: time

    integer iblk, icol, k, l, is, ig, m
    integer l_scavup, l_scavdn
    real(r8) ls, time_of_day, rsdist, cosz, directsol, tsat, rho

    ls = time%solar_longitude()
    time_of_day = time%time_of_day()
    call update_solar_decl_angle(ls)
    rsdist = solar_dist(ls)**2

    ! Calculate solar flux at the current Mars distance.
    do is = 1, l_nspectv
      solar(is) = solar_1au(is) * rsdist
    end do

    do iblk = 1, size(objects)
      associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state)
      do icol = 1, mesh%nlev
        cosz = solar_cos_zenith_angle(mesh%lon(icol), mesh%lat(icol), time_of_day)
        if (cosz >= 1.0e-5_r8) then
          call dsolflux(solar, cosz, gweight, fzerov, state%detau(icol,:,:), directsol)
        else
          directsol = 0
        end if
        state%dnvflux(icol) = state%dndiffv(icol) + directsol
        ! Calculate the ground temperature.
        call tempgr( &
          mesh%lat       (icol                  ), &
          state%ps       (icol                  ), &
          state%t        (icol,mesh%nlev        ), &
          state%q        (icol,mesh%nlev,iMa_vap), &
          state%qcond    (icol,          iMa_vap), &
          state%alsp     (icol                  ), &
          state%npcflag  (icol                  ), &
          state%rhouch   (icol                  ), &
          state%co2ice   (icol                  ), &
          state%fa       (icol                  ), &
          state%dnirflux (icol                  ), &
          state%dnvflux  (icol                  ), &
          state%subflux  (icol                  ), &
          state%gndice   (icol                  ), &
          state%rhosoil  (icol,:                ), &
          state%cpsoil   (icol,:                ), &
          state%scond    (icol,:                ), &
          state%stemp    (icol,:                ), &
          sthick                                 , &
          state%zin      (icol,1                ), &
          state%gt       (icol                  ), &
          state%surfalb  (icol                  )  &
        )
        ! Calculate the average cosine of the solar zenith angle.
        call solarza(                              &
          mesh%lon       (icol                  ), &
          mesh%cos_lat   (icol                  ), &
          mesh%sin_lat   (icol                  ), &
          cos_decl                               , &
          sin_decl                               , &
          time_of_day                            , &
          cosz                                     &
        )
        ! Calculate the CO2 condensation temperature at the surface.
        tsat = 3182.48_r8 / (23.3494_r8 - log(state%ps(icol) / 100.0_r8))
        if (state%gt(icol) < tsat .or. state%co2ice(icol) > 0) then
          state%gt(icol) = tsat
        end if
        ! Calculate the ice cloud optical depth, where present.
        ! FIXME: Is this finished?
        ! Call gridvel.
        call potemp1(            &
          state%ps     (icol  ), &
          state%tstrat (icol  ), &
          state%t      (icol,:), &
          state%aadj           , &
          state%badj           , &
          state%plogadj        , &
          state%pl             , &
          state%om             , &
          state%tl             , &
          state%teta             &
        )
        ! Save initial values for teta, upi, vpi and qpi.

        call coldair(            &
          state%tstrat (icol  ), &
          state%dp     (icol,:), &
          state%pl             , &
          state%tl             , &
          state%gt     (icol  ), &
          state%co2ice (icol  ), &
          state%zin    (icol,:), &
          state%dmadt  (icol  ), &
          state%atmcond(icol,:)  &
        )

        call potemp2(            &
          state%ps     (icol  ), &
          state%tstrat (icol  ), &
          state%t      (icol,:), &
          state%aadj           , &
          state%badj           , &
          state%plogadj        , &
          state%pl             , &
          state%om             , &
          state%tl             , &
          state%teta             &
        )
        ! Radiation calculation.

        ! Aerosol scavenging by CO2 snow fall.
        if (co2scav) then
          l_scavup = 0
          l_scavdn = 0
          do l = 1, mesh%nlev
            if (state%atmcond(icol,l) > 0) then
              l_scavup = l
              exit
            end if
          end do
          do l = mesh%nlev, 1, -1
            if (state%atmcond(icol,l) > 0) then
              l_scavdn = l
              exit
            end if
          end do
          if (l_scavup /= 0 .and. l_scavdn /= 0) then
            do l = l_scavup, l_scavdn
              k = 2 * l + 2
              if (l_scavdn == mesh%nlev) then
                ! If condensation occurs down to the surface, put all aerosols on the surface (in fact, only the cloud mass matters).
                rho = (state%pl(k+1) - state%pl(k-1)) / g
                state%qpig(iMa_vap) = state%qpig(iMa_vap) + scaveff * state%qpi(k,iMa_cld) * rho
                state%qpig(iMa_dt ) = state%qpig(iMa_dt ) + scaveff * state%qpi(k,iMa_dt ) * rho
                state%qpig(iMa_cor) = state%qpig(iMa_cor) + scaveff * state%qpi(k,iMa_cor) * rho
                state%srfdnflx(icol,iMa_dt ) = state%srfdnflx(icol,iMa_dt ) + scaveff * state%qpi(k,iMa_dt ) * rho / dt
                state%srfdnflx(icol,iMa_cor) = state%srfdnflx(icol,iMa_cor) + scaveff * state%qpi(k,iMa_cor) * rho / dt
                state%srfdnflx(icol,iMa_cld) = state%srfdnflx(icol,iMa_cld) + scaveff * state%qpi(k,iMa_cld) * rho / dt
              else
                ! If condensation occurs in a restricted portion, put aerosols in the highest layer unaffected by CO2 condensation.
                do m = 1, naer
                  state%qpi(2*l_scavdn+4,m) = state%qpi(2*l_scavdn+4,m) + scaveff * state%qpi(k,m) * &
                    (state%pl(k+1) - state%pl(k-1)) / (state%pl(2*l_scavdn+5) - state%pl(2*l_scavdn+3))
                end do
              end if
              do m = 1, naer
                state%qpi(k,m) = state%qpi(k,m) * (1 - scaveff)
              end do
              state%atmcond(icol,l) = 0
            end do
          end if
        end if ! co2scav
        if (active_water) then
          do l = 1, mesh%nlev
            k = 2 * l + 2
            state%qh2o(k) = mwratio * state%q(icol,l,iMa_vap)
            state%qh2o(k+1) = state%qh2o(k)
          end do
        else
          do l = 1, mesh%nlev
            k = 2 * l + 2
            state%qh2o(k) = 1.0e-7_r8
            state%qh2o(k+1) = state%qh2o(k)
          end do
        end if
        if (.not. active_dust) then

        end if
        ! Call opt_dst.
      end do
      end associate
    end do

  end subroutine gomars_v1_run

  subroutine gomars_v1_final()

    call gomars_v1_objects_final()
    call gomars_v1_rad_final()
    call gomars_v1_pbl_final()
    call gomars_v1_damp_final()

  end subroutine gomars_v1_final

  subroutine gomars_v1_d2p()

    integer iblk, icol, k

    ! Calculate ts.

    do iblk = 1, size(objects)
      associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state)
      do icol = 1, mesh%ncol

      end do
      end associate
    end do

  end subroutine gomars_v1_d2p

  subroutine gomars_v1_p2d()

  end subroutine gomars_v1_p2d

end module gomars_v1_driver_mod