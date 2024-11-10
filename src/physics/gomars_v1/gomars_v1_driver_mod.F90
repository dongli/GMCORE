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

  use const_mod
  use namelist_mod
  use vert_coord_mod
  use gomars_v1_namelist_mod
  use gomars_v1_objects_mod
  use gomars_v1_tracers_mod
  use gomars_v1_orbit_mod
  use gomars_v1_rad_mod
  use gomars_v1_pbl_mod
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
    call gomars_v1_damp_init(mesh(1)%nlev)

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
      call ini_optdst(l_nspectv, l_nspecti, l_nrefv, mesh%nlev, &
                      qextv, qscatv, gv, qexti, qscati, gi, &
                      state%qxvdst, state%qsvdst, state%gvdst, &
                      state%qxidst, state%qsidst, state%gidst, &
                      state%qextrefdst)
      call ini_optcld(l_nspectv, l_nspecti, mesh%nlev, &
                      state%qxvcld, state%qsvcld, state%gvcld, &
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

  subroutine gomars_v1_run()

    ! Calculate solar flux at the current Mars distance.

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