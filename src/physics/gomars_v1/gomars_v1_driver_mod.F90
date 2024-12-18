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
  use gomars_v1_pbl_mod
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

    dt      = dt_phys
    dt_mp   = dt_phys / nsplit
    nlev    = mesh(1)%nlev
    nlayrad = nlev + 1
    nlevrad = nlev + 2

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
          state%tstrat    (icol) = state%t(icol,mesh%nlev)
          state%co2ice_sfc(icol) = 0
          state%tg        (icol) = state%t(icol,mesh%nlev)
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
          state%irflx_sfc_dn(icol) = 1
          state%vsdif_sfc_dn(icol) = 0
          do is = 1, nspectv
            do ig = 1, ngauss
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

    integer iblk, icol, k, l, m, substep
    real(r8) ls, time_of_day, tsat, rhodz
    real(r8) nfluxtopv, nfluxtopi, diffvt, albi, sunlte

    ls = time%solar_longitude()
    time_of_day = time%time_of_day()
    call update_solar_decl_angle(ls)
    call update_solar(ls)    

    blocks: do iblk = 1, size(objects)
      associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state, tend => objects(iblk)%tend)
      columns: do icol = 1, mesh%nlev
        ! ----------------------------------------------------------------------
        ! Calculate the direct solar flux.
        state%cosz(icol) = solar_cos_zenith_angle(mesh%lon(icol), mesh%lat(icol), time_of_day)
        if (state%cosz(icol) >= 1.0e-5_r8) then
          call dsolflux(                     &
            solar                          , & ! in
            state%cosz        (icol       ), & ! in
            gweight                        , & ! in
            fzerov                         , & ! in
            state%detau       (icol,:,:   ), & ! in
            state%solar_sfc_dn(icol       )  & ! out
          )
        else
          state%solar_sfc_dn(icol) = 0
        end if
        state%vsflx_sfc_dn(icol) = state%vsdif_sfc_dn(icol) + state%solar_sfc_dn(icol)
        ! Recalculate the average cosine of the solar zenith angle.
        call solarza(                        &
          mesh%lon            (icol       ), & ! in
          mesh%cos_lat        (icol       ), & ! in
          mesh%sin_lat        (icol       ), & ! in
          cos_decl                         , & ! in
          sin_decl                         , & ! in
          time_of_day                      , & ! in
          state%cosz          (icol       )  & ! out
        )
        ! ----------------------------------------------------------------------
        ! Calculate the ground temperature.
        call tempgr(                         &
          mesh%lat            (icol       ), & ! in
          state%ps            (icol       ), & ! in
          state%t             (icol,nlev  ), & ! in
          state%q             (icol,nlev,:), & ! in
          state%tm_sfc        (icol,     :), & ! in
          state%alsp          (icol       ), & ! in
          state%polarcap      (icol       ), & ! in
          state%rhouch        (icol       ), & ! in
          state%co2ice_sfc    (icol       ), & ! inout
          state%ht_sfc        (icol       ), & ! in
          state%irflx_sfc_dn  (icol       ), & ! in
          state%vsflx_sfc_dn  (icol       ), & ! in
          state%h2osub_sfc    (icol       ), & ! out
          state%h2oice_sfc    (icol       ), & ! out
          state%rhosoil       (icol,:     ), & ! in
          state%cpsoil        (icol,:     ), & ! in
          state%scond         (icol,:     ), & ! in
          state%stemp         (icol,:     ), & ! inout
          state%zin           (icol,1     ), & ! in
          state%tg            (icol       ), & ! out
          state%als           (icol       )  & ! out
        )
        ! ----------------------------------------------------------------------
        ! Calculate potential temperature on all levels.
        call potemp1(                        &
          state%ps            (icol       ), & ! in
          state%tstrat        (icol       ), & ! in
          state%t             (icol,:     ), & ! in
          state%aadj                       , & ! in
          state%badj                       , & ! in
          state%plogadj                    , & ! out
          state%pl                         , & ! out
          state%om                         , & ! out
          state%tl                         , & ! out
          state%teta                         & ! out
        )
        call coldair(                        &
          state%tstrat        (icol       ), & ! inout
          state%dp            (icol,:     ), & ! in
          state%pl                         , & ! in
          state%tl                         , & ! inout
          state%tg            (icol       ), & ! inout
          state%co2ice_sfc    (icol       ), & ! inout
          state%q             (icol,:,:   ), & ! inout
          state%tm_sfc        (icol,  :   ), & ! out
          state%tmflx_sfc_dn  (icol,  :   ), & ! out
          state%zin           (icol,:     ), & ! in
          state%dmadt         (icol       ), & ! out
          state%atmcond       (icol,:     )  & ! inout
        )
        call potemp2(                        &
          state%ps            (icol       ), & ! in
          state%tstrat        (icol       ), & ! in
          state%t             (icol,:     ), & ! in
          state%aadj                       , & ! in
          state%badj                       , & ! in
          state%plogadj                    , & ! out
          state%pl                         , & ! out
          state%om                         , & ! out
          state%tl                         , & ! inout
          state%teta                         & ! out
        )
        ! ----------------------------------------------------------------------
        ! Radiation calculation
        ! Transfer pressure and temperature onto radiation levels.
        call fillpt(                         &
          state%pl                         , & ! in
          state%tl                         , & ! in
          state%tg            (icol      ) , & ! in
          state%tstrat        (icol      ) , & ! in
          state%plev_rad                   , & ! out
          state%tlev_rad                   , & ! out
          state%pmid_rad                   , & ! out
          state%tmid_rad                     & ! out
        )
        if (active_water) then
          do k = 1, mesh%nlev
            l = 2 * k + 2
            state%qh2o(l  ) = mwratio * state%q(icol,k,iMa_vap)
            state%qh2o(l+1) = state%qh2o(l)
          end do
        else
          do k = 1, mesh%nlev
            l = 2 * k + 2
            state%qh2o(l  ) = 1.0e-7_r8
            state%qh2o(l+1) = state%qh2o(l)
          end do
        end if
        if (.not. active_dust) then
          ! Semi-prescribed dust.
        else
          call opt_dst(                     &
            state%q           (icol,:,:  ), & ! in
            state%pl                      , & ! in
            state%qxvdst                  , & ! out
            state%qsvdst                  , & ! out
            state%gvdst                   , & ! out
            state%qxidst                  , & ! out
            state%qsidst                  , & ! out
            state%gidst                   , & ! out
            state%qextrefdst              , & ! out
            state%taurefdst               , & ! out
            state%taudst      (icol,2    )  & ! out
          )
          do l = 1, 3
            state%taurefdst(l) = 0
            state%taucum   (l) = 0
          end do
          do l = 4, 2 * mesh%nlev + 3
            state%taucum(l) = state%taucum(l-1) + state%taurefdst(l)
          end do
        end if
        ! Fill special bottom radiation level to zero.
        state%taurefdst(2*mesh%nlev+4) = 0
        state%tausurf(icol) = state%taucum(2*mesh%nlev+3)
        if (cloudon) then
          call opt_cld(                     &
            state%q           (icol,:,:  ), & ! in
            state%pl                      , & ! in
            state%qxvcld                  , & ! out
            state%qsvcld                  , & ! out
            state%gvcld                   , & ! out
            state%qxicld                  , & ! out
            state%qsicld                  , & ! out
            state%gicld                   , & ! out
            state%qextrefcld              , & ! out
            state%taurefcld               , & ! out
            state%taucld      (icol,2    )  & ! out
          )
        else
          state%taurefcld = 0
        end if
        if (state%cosz(icol) >= 1.0e-5) then
          ! Check for ground ice. Change albedo if there is any ice.
          state%als(icol) = state%alsp(icol)
          if (state%co2ice_sfc(icol) > 0) then
            state%als(icol) = merge(alices, alicen, mesh%lat(icol) < 0)
          else if (albfeed .and. state%tm_sfc(icol,iMa_vap) > icethresh_kgm2 .and. state%polarcap(icol)) then
            state%als(icol) = icealb
          end if
          ! Calculate optical depth due to all sources in the visible bands.
          call optcv(                       &
            state%plev_rad                , & ! in
            state%pmid_rad                , & ! in
            state%tmid_rad                , & ! in
            state%qh2o                    , & ! in
            state%qxvdst                  , & ! in
            state%qsvdst                  , & ! in
            state%gvdst                   , & ! in
            state%qxvcld                  , & ! in
            state%qsvcld                  , & ! in
            state%gvcld                   , & ! in
            state%qextrefcld              , & ! in
            state%wbarv                   , & ! out
            state%cosbv                   , & ! out
            state%dtauv                   , & ! out
            state%tauv                    , & ! out
            state%taucumv                 , & ! out
            state%taugsurf                , & ! out
            state%taurefdst               , & ! inout
            state%taurefcld                 & ! inout
          )
          ! Calculate the fluxes in the visible bands.
          call sfluxv(                      &
            state%dtauv                   , & ! in
            state%tauv                    , & ! in
            state%taucumv                 , & ! in
            state%taugsurf                , & ! in
            state%cosz        (icol      ), & ! in
            state%als         (icol      ), & ! in
            state%wbarv                   , & ! in
            state%cosbv                   , & ! in
            state%fluxupv                 , & ! out
            state%fluxdnv                 , & ! out
            state%fmnetv                  , & ! out
            nfluxtopv                     , & ! out
            diffvt                        , & ! out
            state%detau       (icol,:,:  )  & ! out
          )
          state%suntot(3) = state%fmnetv(1) - nfluxtopv
          do k = 2, nlayrad
            l = 2 * k + 1
            state%suntot(l) = state%fmnetv(k) - state%fmnetv(k-1)
          end do
        else
          ! If the sun is down, no solar flux, nor downward flux.
          do k = 1, nlayrad
            l = 2 * k + 1
            state%suntot (l) = 0
            state%fluxdnv(k) = 0
          end do
          diffvt = 0
          state%fluxupv(1) = 0
          state%fluxdnv(1) = 0
          state%fluxupv(nlayrad) = 0
          state%fluxdnv(nlayrad) = 0
        end if
        state%vsflx_sfc_dn(icol) = state%fluxdnv(nlayrad)
        state%vsdif_sfc_dn(icol) = diffvt
        state%fuptopv     (icol) = state%fluxupv(1)
        state%fdntopv     (icol) = state%fluxdnv(1)
        state%fupsfcv     (icol) = state%fluxupv(nlayrad)
        state%fdnsfcv     (icol) = state%fluxdnv(nlayrad)
        ! Set up and solve for the infrared fluxes.
        ! Check for ground ice, and change albedo if there is any ice.
        albi = 1 - egognd
        if (state%co2ice_sfc(icol) > 0) then
          albi = 1 - merge(egoco2s, egoco2n, mesh%lat(icol) < 0)
        end if
        ! Calculate the optical depth due to all sources in the infrared bands.
        call optci(                         &
          state%plev_rad                  , & ! in
          state%pmid_rad                  , & ! in
          state%tmid_rad                  , & ! in
          state%qh2o                      , & ! in
          state%qxidst                    , & ! in
          state%qsidst                    , & ! in
          state%gidst                     , & ! in
          state%qextrefdst                , & ! in
          state%qxicld                    , & ! in
          state%qsicld                    , & ! in
          state%gicld                     , & ! in
          state%qextrefcld                , & ! in
          state%wbari                     , & ! out
          state%cosbi                     , & ! out
          state%dtaui                     , & ! out
          state%taui                      , & ! out
          state%taucumi                   , & ! out
          state%taugsurf                  , & ! out
          state%taurefdst                 , & ! inout
          state%taurefcld                   & ! inout
        )
        ! Calculate the fluxes in the infrared bands.
        call sfluxi(                        &
          state%plev_rad                  , & ! in
          state%tlev_rad                  , & ! in
          state%dtaui                     , & ! in
          state%taucumi                   , & ! in
          state%taugsurf                  , & ! in
          albi                            , & ! in
          state%wbari                     , & ! in
          state%cosbi                     , & ! in
          state%fluxupi                   , & ! out
          state%fluxdni                   , & ! out
          state%fmneti                    , & ! out
          nfluxtopi                         & ! out
        )
        state%irtot(3) = state%fmneti(1) - nfluxtopi
        do k = 2, nlayrad
          state%irtot(2*k+1) = state%fmneti(k) - state%fmneti(k-1)
        end do
        state%irflx_sfc_dn(icol) = state%fluxdni(nlayrad)
        state%fluxsfc     (icol) = (1 - state%als(icol)) * state%fluxdnv(nlayrad)
        state%fuptopi     (icol) = state%fluxupi(1)
        state%fupsfci     (icol) = state%fluxupi(nlayrad)
        state%fdnsfci     (icol) = state%fluxdni(nlayrad)
        ! Change atmospheric temperature for solar and infrared heating.
        ! Include heating rates in boundary layer scheme.
        ! Store total radiative heating rates in qrad and later pass them into the new boundary layer scheme.
        ! Also add in the non-LTE correction.
        do l = 2, 2 * mesh%nlev + 2, 2
          k = (l - 2) / 2
          sunlte = state%suntot(l+1) * 2.2e2_r8 * state%plev_rad(l) / (1 + 2.2e2_r8 * state%plev_rad(l))
          ! FIXME: Could we use state%rho * state%dz?
          rhodz = (state%pl(l+1) - state%pl(l-1)) / g
          state%ht_rad(icol,k) = (sunlte + state%irtot(l+1)) / (cpd * rhodz * state%om(l))
        end do
        ! Non-LTE fudge only the stratosphere.
        l = 2
        sunlte = state%suntot(l+1) * 2.2e2_r8 * state%plev_rad(l) / (1 + 2.2e2_r8 * state%plev_rad(l))
        rhodz = (state%pl(l+1) - state%pl(l-1)) / g
        state%teta(l) = state%teta(l) + dt * (sunlte + state%irtot(l+1)) / (cpd * rhodz * state%om(l))
        ! ----------------------------------------------------------------------
        ! Planetary boundary layer calculation
        state%z0(icol) = z00
        if (state%co2ice_sfc(icol) > 0) state%z0(icol) = 1.0e-4_r8
        if (state%tm_sfc(icol,iMa_vap) > 100 .or. state%polarcap(icol)) state%z0(icol) = 1.0e-4_r8
        call newpbl(                    &
          state%z0          (icol    ), &
          state%tg          (icol    ), &
          state%ht_rad      (icol,:  ), &
          state%ps          (icol    ), &
          state%ts          (icol    ), &
          state%polarcap    (icol    ), &
          state%u           (icol,:  ), &
          state%v           (icol,:  ), &
          state%pt          (icol,:  ), &
          state%pt_lev      (icol,:  ), &
          state%q           (icol,:,:), &
          state%dp_dry      (icol,:  ), &
          state%z           (icol,:  ), &
          state%dz          (icol,:  ), &
          state%z_lev       (icol,:  ), &
          state%dz_lev      (icol,:  ), &
          state%shr2        (icol,:  ), &
          state%ri          (icol,:  ), &
          state%km          (icol,:  ), &
          state%kh          (icol,:  ), &
          state%ustar       (icol    ), &
          state%tstar       (icol    ), &
          state%taux        (icol    ), &
          state%tauy        (icol    ), &
          state%ht_pbl      (icol    ), &
          state%rhouch      (icol    ), &
          state%tm_sfc      (icol,:  ), &
          state%h2osub_sfc  (icol    )  &
        )
        ! ----------------------------------------------------------------------
        ! Microphysics calculation
        if (microphysics) then
          do substep = 1, nsplit
            ! call microphys
          end do
        end if
        ! ----------------------------------------------------------------------
        ! Convection adjustment
      end do columns
      end associate
    end do blocks

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