! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

program gmcore_driver

  use flogger
  use namelist_mod
  use const_mod
  use block_mod
  use latlon_parallel_mod
  use initial_mod
  use restart_mod
  use gmcore_mod
  use mountain_zonal_flow_test_mod
  use rossby_haurwitz_wave_test_mod
  use jet_zonal_flow_test_mod
  use steady_geostrophic_flow_test_mod
  use cross_pole_flow_test_mod
  use shallow_water_waves_test_mod
  use vortex_erosion_test_mod
  use splash_test_mod
  use colliding_modons_test_mod
  use steady_state_test_mod
  use rossby_haurwitz_wave_3d_test_mod
  use mountain_wave_test_mod
  use baroclinic_wave_test_mod
  use held_suarez_test_mod
  use steady_state_pgf_test_mod
  use ksp15_test_mod
  use dcmip31_test_mod
  use mars_cold_run_mod
  use tropical_cyclone_test_mod
#ifdef HAS_LAPACK
  use supercell_test_mod
#endif
  use prepare_mod

  implicit none

  interface
    subroutine set_ic_interface(block)
      import block_type
      type(block_type), intent(inout), target :: block
    end subroutine set_ic_interface
  end interface
  procedure(set_ic_interface), pointer :: set_ic

  character(256) namelist_path
  integer iblk

  call get_command_argument(1, namelist_path)
  if (namelist_path == '') then
    call log_error('You should give a namelist file path!')
  end if

  call parse_namelist(namelist_path)

  call gmcore_init_stage0()

  select case (test_case)
  case ('swm_sp')
    call splash_test_set_params()
  case ('swm_cm', 'cm')
    call colliding_modons_test_set_params()
  case ('ksp15_01', 'ksp15_02')
    call ksp15_test_set_params()
  case ('ss_pgf')
    call steady_state_pgf_test_set_params()
  case ('dcmip31')
    call dcmip31_test_set_params()
#ifdef HAS_LAPACK
  case ('sc')
    call supercell_test_set_params()
#endif
  end select

  call gmcore_init_stage1(namelist_path)

  select case (test_case)
  case ('bw')
    call baroclinic_wave_test_init(namelist_path)
  case ('mz')
    call mountain_wave_test_init()
  case ('tc')
    call tropical_cyclone_test_init()
#ifdef HAS_LAPACK
  case ('sc')
    call supercell_test_init(namelist_path)
#endif
  end select

  if (initial_file == 'N/A' .and. test_case == 'N/A' .and. .not. restart) then
    call prepare_topo()
    call prepare_tracers()
  else if (initial_file /= 'N/A') then
    call initial_read_init()
  end if

  call gmcore_init_stage2(namelist_path)

  if (restart) then
    call restart_read()
  else if (initial_file /= 'N/A') then
    call initial_read()
  else if (bkg_file /= 'N/A') then
    call prepare_bkg()
    call prepare_final() ! Release memory for preparation.
  else
    select case (test_case)
    case ('swm_mz')
      set_ic => mountain_zonal_flow_test_set_ic
    case ('swm_rh')
      set_ic => rossby_haurwitz_wave_test_set_ic
    case ('swm_jz')
      set_ic => jet_zonal_flow_test_set_ic
    case ('swm_sg')
      set_ic => steady_geostrophic_flow_test_set_ic
    case ('swm_cp')
      set_ic => cross_pole_flow_test_set_ic
    case ('swm_sw')
      set_ic => shallow_water_waves_test_set_ic
    case ('swm_vr')
      set_ic => vortex_erosion_test_set_ic
    case ('swm_sp')
      set_ic => splash_test_set_ic
    case ('swm_cm', 'cm')
      set_ic => colliding_modons_test_set_ic
    case ('ss')
      set_ic => steady_state_test_set_ic
    case ('ss_pgf')
      set_ic => steady_state_pgf_test_set_ic
    case ('rh')
      set_ic => rossby_haurwitz_wave_3d_test_set_ic
    case ('mz')
      set_ic => mountain_wave_test_set_ic
    case ('bw')
      set_ic => baroclinic_wave_test_set_ic
    case ('hs')
      set_ic => held_suarez_test_set_ic
    case ('ksp15_01')
      set_ic => ksp15_01_test_set_ic
    case ('ksp15_02')
      set_ic => ksp15_02_test_set_ic
    case ('dcmip31')
      set_ic => dcmip31_test_set_ic
    case ('mars_cold_run')
      set_ic => mars_cold_run_set_ic
    case ('tc')
      set_ic => tropical_cyclone_test_set_ic
#ifdef HAS_LAPACK
    case ('sc')
      set_ic => supercell_test_set_ic
#endif
    case default
      if (proc%is_root()) call log_error('Unknown test case ' // trim(test_case) // '!')
    end select

    if (proc%is_model()) then
      do iblk = 1, size(blocks)
        call set_ic(blocks(iblk))
      end do
    end if
  end if

  call gmcore_init_stage3()

  call gmcore_run()

  call gmcore_final()

end program gmcore_driver
