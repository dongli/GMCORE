! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module prepare_mod

  use const_mod
  use namelist_mod
  use latlon_parallel_mod
  use process_mod
  use block_mod
  use topo_reader_mod
  use latlon_topo_mod
  use latlon_bkg_mod
  use ref_mod
  use operators_mod
  use tracer_mod

  implicit none

contains

  subroutine prepare_topo()

    integer iblk

    if (.not. use_aqua_planet) then
      call topo_reader_run(topo_file, min_lon, max_lon, min_lat, max_lat)
      if (proc%is_model()) then
        do iblk = 1, size(blocks)
          call latlon_topo_regrid(blocks(iblk))
        end do
        if (use_topo_smooth) then
          do iblk = 1, size(blocks)
            call latlon_topo_smooth(blocks(iblk))
          end do
        end if
        call ref_calc_ps()
      end if
    end if

  end subroutine prepare_topo

  subroutine prepare_tracers()

    if (ideal_dry_core) return

    select case (planet)
    case ('earth')
      select case (bkg_type)
      case ('era5')
        ! Add water vapor tracer for testing.
        if (idx_qv == 0 .and. physics_suite == 'N/A') then
          call tracer_add('moist', dt_adv, 'qv', 'Water vapor', 'kg kg-1')
          call tracer_add('moist', dt_adv, 'qc', 'Cloud water', 'kg kg-1')
          call tracer_add('moist', dt_adv, 'qi', 'Cloud ice'  , 'kg kg-1')
          call tracer_add('moist', dt_adv, 'qr', 'Rain water' , 'kg kg-1')
          call tracer_add('moist', dt_adv, 'qs', 'Snow water' , 'kg kg-1')
        end if
      end select
    end select

  end subroutine prepare_tracers

  subroutine prepare_bkg()

    integer iblk, itime, cyc

    call latlon_bkg_read(min_lon, max_lon, min_lat, max_lat)
    if (proc%is_model()) then
      ! Stage 1:
      call latlon_bkg_regrid_phs()
      call latlon_bkg_calc_ph()
      do cyc = 1, 3
        call latlon_bkg_regrid_wet_qv(mute=cyc /= 1)
        call latlon_bkg_regrid_wet_qc(mute=cyc /= 1)
        call latlon_bkg_regrid_wet_qi(mute=cyc /= 1)
        call latlon_bkg_regrid_wet_qr(mute=cyc /= 1)
        call latlon_bkg_regrid_wet_qs(mute=cyc /= 1)
        do iblk = 1, size(blocks)
          call tracer_calc_qm(blocks(iblk))
        end do
        ! Stage 2:
        call latlon_bkg_calc_mgs(mute=cyc /= 1)
        call latlon_bkg_calc_mg(mute=cyc /= 1)
        call latlon_bkg_calc_dry_qv(mute=cyc /= 1)
        call latlon_bkg_calc_dry_qc(mute=cyc /= 1)
        call latlon_bkg_calc_dry_qi(mute=cyc /= 1)
        call latlon_bkg_calc_dry_qr(mute=cyc /= 1)
        call latlon_bkg_calc_dry_qs(mute=cyc /= 1)
        do iblk = 1, size(blocks)
          call tracer_calc_qm(blocks(iblk))
          call calc_ph(blocks(iblk), blocks(iblk)%dstate(1))
        end do
      end do
      ! Stage 3:
      call latlon_bkg_regrid_u()
      call latlon_bkg_regrid_v()
      if (prepare_regrid_gz) then
        call latlon_bkg_regrid_gz()
        call latlon_bkg_calc_t()
      else
        call latlon_bkg_regrid_t()
      end if
      call latlon_bkg_calc_tv()
      call latlon_bkg_calc_pt()
      if (.not. prepare_regrid_gz) then
        do iblk = 1, size(blocks)
          do itime = lbound(blocks(iblk)%dstate, 1), ubound(blocks(iblk)%dstate, 1)
            call blocks(iblk)%dstate(itime)%gz_lev%copy(blocks(iblk)%static%gzs, k=blocks(iblk)%mesh%half_kde, with_halo=.true.)
          end do
          call calc_gz(blocks(iblk), blocks(iblk)%dstate(1))
          call fill_halo(blocks(iblk)%dstate(1)%gz_lev)
        end do
      end if
    end if

  end subroutine prepare_bkg

  subroutine prepare_run()

  end subroutine prepare_run

  subroutine prepare_final()

    call topo_reader_final()
    call latlon_bkg_final()

  end subroutine prepare_final

end module prepare_mod
