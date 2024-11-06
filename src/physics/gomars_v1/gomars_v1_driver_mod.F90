module gomars_v1_driver_mod

  use const_mod
  use gomars_v1_objects_mod
  use nasa_rad_mod

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

  subroutine gomars_v1_init_stage2(namelist_path, mesh, dt_adv, dt_phys, input_ngroup, model_root)

    character(*), intent(in) :: namelist_path
    type(physics_mesh_type), intent(in), target :: mesh(:)
    real(r8), intent(in) :: dt_adv
    real(r8), intent(in) :: dt_phys
    integer , intent(in) :: input_ngroup
    character(*), intent(in), optional :: model_root

    call gomars_v1_objects_init(mesh)
    call nasa_rad_init()

  end subroutine gomars_v1_init_stage2

  subroutine gomars_v1_init_stage3()

    integer iblk, icol, k, is, ig

    ! Check model top pressure. It must be lower than ptop in nasa_rad_mod.

    do iblk = 1, size(objects)
      associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state)
      ! ini_optdst:
      do k = 1, mesh%nlev + 1
        state%qextrefdst(k) = qextv(l_nrefv)
      end do
      do is = 1, l_nspectv
        do k = 1, mesh%nlev + 1
          state%qxvdst(k,is) = qextv (is)
          state%qsvdst(k,is) = qscatv(is)
          state%gvdst (k,is) = gv    (is)
        end do
      end do
      do is = 1, l_nspecti
        do k = 1, mesh%nlev + 1
          state%qxidst(k,is) = qexti (is)
          state%qsidst(k,is) = qscati(is)
          state%gidst (k,is) = gi    (is)
        end do
      end do
      ! ini_optcld:
      do k = 1, mesh%nlev + 1
        state%qextrefcld(k) = 1
        state%taurefcld (k) = 0
      end do
      do is = 1, l_nspectv
        do k = 1, mesh%nlev + 1
          state%qxvcld(k,is) = 0
          state%qsvcld(k,is) = 0
          state%gvcld (k,is) = 0
        end do
      end do
      do is = 1, l_nspecti
        do k = 1, mesh%nlev + 1
          state%qxicld(k,is) = 0
          state%qsicld(k,is) = 0
          state%gicld (k,is) = 0
        end do
      end do
      ! firstcomp3:
      do icol = 1, mesh%ncol
        state%dnirflux(icol) = 1
        state%dndiffv (icol) = 0
        do is = 1, l_nspectv
          do ig = 1, l_ngauss
            state%detau(icol,is,ig) = 0.1_r8
          end do
        end do
      end do
      end associate
    end do

  end subroutine gomars_v1_init_stage3

  subroutine gomars_v1_run()

    ! Calculate solar flux at the current Mars distance

  end subroutine gomars_v1_run

  subroutine gomars_v1_final()

    call gomars_v1_objects_final()
    call nasa_rad_final()

  end subroutine gomars_v1_final

  subroutine gomars_v1_d2p()

    ! Calculate ts.

  end subroutine gomars_v1_d2p

  subroutine gomars_v1_p2d()

  end subroutine gomars_v1_p2d

end module gomars_v1_driver_mod