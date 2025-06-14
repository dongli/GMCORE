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

module gomars_v2_namelist_mod

  use fiona
  use flogger
  use string
  use process_mod
  use gomars_v2_const_mod, only: r8, rho_ice, ice_thresh_kgm2

  implicit none

  character(1024) :: kcoef_file            = ''
  character(1024) :: dust_optics_file      = ''
  character(1024) :: cld_optics_file       = ''
  character(1024) :: albedo_file           = ''
  character(1024) :: thermal_inertia_file  = ''
  character(1024) :: gnd_ice_file          = ''

  logical :: active_water                  = .false.
  logical :: active_dust                   = .false.
  logical :: albedo_feedback               = .false.
  logical :: water_ice_latent_heat         = .false.

  real(r8) :: ice_albedo                   = 0.4_r8
  ! Ice depth threshold required to change surface albedo (um)
  real(r8) :: ice_thresh_depth             = 5.0_r8

  ! Number of soil layers
  integer :: nlev_soil                     = 40

  namelist /gomars_v2_control/ &
    kcoef_file               , &
    dust_optics_file         , &
    cld_optics_file          , &
    albedo_file              , &
    thermal_inertia_file     , &
    gnd_ice_file             , &
    active_water             , &
    active_dust              , &
    albedo_feedback          , &
    ice_albedo               , &
    ice_thresh_depth         , &
    nlev_soil

contains

  subroutine gomars_v2_parse_namelist(file_path, model_root)

    character(*), intent(in) :: file_path
    character(*), intent(in), optional :: model_root

    integer ierr
    logical is_exist

    ice_thresh_kgm2 = ice_thresh_depth * 1.0e-6_r8 * rho_ice

    open(10, file=file_path, status='old')
    read(10, nml=gomars_v2_control, iostat=ierr)
    if (ierr /= 0 .and. ierr /= -1) then
      if (proc%is_root()) call log_error('There is error in gomars_v2_control namelist in ' // trim(file_path) // '!')
    end if
    close(10)

    if (kcoef_file == '' .and. present(model_root)) then
      kcoef_file = abspath(trim(model_root) // '/data/mars/nasa_kcoef.nc')
    end if
    inquire(file=kcoef_file, exist=is_exist)
    if (.not. is_exist) then
      if (proc%is_root()) call log_error('kcoef_file ' // trim(kcoef_file) // ' does not exist!')
    end if

    if (dust_optics_file == '' .and. present(model_root)) then
      dust_optics_file = abspath(trim(model_root) // '/data/mars/nasa_dust_optics.nc')
    end if
    inquire(file=dust_optics_file, exist=is_exist)
    if (.not. is_exist) then
      if (proc%is_root()) call log_error('dust_optics_file ' // trim(dust_optics_file) // ' does not exist!')
    end if

    if (cld_optics_file == '' .and. present(model_root)) then
      cld_optics_file = abspath(trim(model_root) // '/data/mars/nasa_cld_optics.nc')
    end if
    inquire(file=cld_optics_file, exist=is_exist)
    if (.not. is_exist) then
      if (proc%is_root()) call log_error('cld_optics_file ' // trim(cld_optics_file) // ' does not exist!')
    end if

    if (albedo_file == '' .and. present(model_root)) then
      albedo_file = abspath(trim(model_root) // '/data/mars/mgs_albedo.nc')
    end if
    inquire(file=albedo_file, exist=is_exist)
    if (.not. is_exist) then
      if (proc%is_root()) call log_error('albedo_file ' // trim(albedo_file) // ' does not exist!')
    end if

    if (thermal_inertia_file == '' .and. present(model_root)) then
      thermal_inertia_file = abspath(trim(model_root) // '/data/mars/thermal_interia_2011.nc')
    end if
    inquire(file=thermal_inertia_file, exist=is_exist)
    if (.not. is_exist) then
      if (proc%is_root()) call log_error('thermal_inertia_file ' // trim(thermal_inertia_file) // ' does not exist!')
    end if

    if (gnd_ice_file == '' .and. present(model_root)) then
      gnd_ice_file = abspath(trim(model_root) // '/data/mars/gnd_ice_map.nc')
    end if
    inquire(file=gnd_ice_file, exist=is_exist)
    if (.not. is_exist) then
      if (proc%is_root()) call log_error('gnd_ice_file ' // trim(gnd_ice_file) // ' does not exist!')
    end if

  end subroutine gomars_v2_parse_namelist

end module gomars_v2_namelist_mod
