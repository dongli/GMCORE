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

  subroutine prepare_bkg()

    integer iblk

    call latlon_bkg_read(min_lon, max_lon, min_lat, max_lat)
    if (proc%is_model()) then
      call latlon_bkg_regrid_mgs()
      call latlon_bkg_calc_mg()
      call latlon_bkg_regrid_qv()
      call latlon_bkg_calc_ph()
      call latlon_bkg_regrid_pt()
      call latlon_bkg_regrid_u()
      call latlon_bkg_regrid_v()
      do iblk = 1, size(blocks)
        call tracer_calc_qm(blocks(iblk))
      end do
    end if

  end subroutine prepare_bkg

  subroutine prepare_run()

  end subroutine prepare_run

  subroutine prepare_final()

    call topo_reader_final()
    call latlon_bkg_final()

  end subroutine prepare_final

end module prepare_mod
