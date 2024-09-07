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
          call tracer_add('moist', dt_adv, 'one', 'Uniform one test tracer', 'kg kg-1')
        end if
      end select
    end select

  end subroutine prepare_tracers

  subroutine prepare_bkg()

    integer iblk

    call latlon_bkg_read(min_lon, max_lon, min_lat, max_lat)
    if (proc%is_model()) then
      call latlon_bkg_regrid_phs()
      call latlon_bkg_calc_ph()
      call latlon_bkg_regrid_wet_qv()
      call latlon_bkg_regrid_wet_qc()
      call latlon_bkg_regrid_wet_qi()
      call latlon_bkg_regrid_wet_qr()
      call latlon_bkg_regrid_wet_qs()
      do iblk = 1, size(blocks)
        call tracer_calc_qm(blocks(iblk))
      end do
      call latlon_bkg_calc_mgs()
      call latlon_bkg_calc_mg()
      call latlon_bkg_regrid_u()
      call latlon_bkg_regrid_v()
      call latlon_bkg_regrid_t()
      ! Change wet mixing ratio to dry mixing ratio.
      call latlon_bkg_calc_dry_qv()
      call latlon_bkg_calc_dry_qc()
      call latlon_bkg_calc_dry_qi()
      call latlon_bkg_calc_dry_qr()
      call latlon_bkg_calc_dry_qs()
      do iblk = 1, size(blocks)
        call tracer_calc_qm(blocks(iblk))
        call calc_gz_lev(blocks(iblk), blocks(iblk)%dstate(1))
        tracers(iblk)%q%d(:,:,:,ntracers) = 1.0_r8
      end do
      call latlon_bkg_calc_pt()
    end if

  end subroutine prepare_bkg

  subroutine prepare_run()

  end subroutine prepare_run

  subroutine prepare_final()

    call topo_reader_final()
    call latlon_bkg_final()

  end subroutine prepare_final

end module prepare_mod
