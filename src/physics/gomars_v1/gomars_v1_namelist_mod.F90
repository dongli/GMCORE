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

module gomars_v1_namelist_mod

  use gomars_v1_const_mod

  implicit none

  real(r8) :: psf              = 701.0_r8  ! Reference surface pressure (Pa)
  real(r8) :: icealb           = 0.4_r8    ! Albedo of surface water ice when albfeed is true.
  real(r8) :: icethresh_depth  = 5.0_r8    ! Ice depth threshold required to change albedo (um)
  integer  :: nsplit           = 1         ! Split number for microphysics
  logical  :: microphysics     = .false.
  logical  :: cloudon          = .false.
  logical  :: active_dust      = .false.
  logical  :: co2scav          = .false.
  logical  :: active_water     = .false.
  logical  :: albfeed          = .false.
  logical  :: latent_heat      = .true.

  logical  :: use_ddl          = .true.
  ! Dust devil lifting efficiency
  real(r8) :: alpha_d          = 5.0e-11_r8

  character(8) :: wsl_scheme   = 'newman'  ! newman or kmh
  logical  :: use_wsl_newman   = .true.
  logical  :: use_wsl_kmh      = .false.
  real(r8) :: tau_thresh       = 0.04_r8   ! Wind stress threshold for wind stress dust lifting (Pa)
  ! Wind stress lifting efficiency
  real(r8) :: alpha_n          = 5.0e-5_r8

  namelist /gomars_v1_control/ &
    psf                      , &
    icealb                   , &
    icethresh_depth          , &
    nsplit                   , &
    cloudon                  , &
    active_dust              , &
    co2scav                  , &
    active_water             , &
    albfeed                  , &
    latent_heat              , &
    use_ddl                  , &
    alpha_d                  , &
    wsl_scheme               , &
    use_wsl_newman           , &
    use_wsl_kmh              , &
    tau_thresh               , &
    alpha_d

contains

  subroutine gomars_v1_parse_namelist(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path, status='old')
    read(10, nml=gomars_v1_control)
    close(10)

    psl = psf

    icethresh_kgm2 = icethresh_depth * rho_ice * 1.0E-6_r8

    use_wsl_newman = wsl_scheme == 'newman'
    use_wsl_kmh    = wsl_scheme == 'kmh'

  end subroutine gomars_v1_parse_namelist

end module gomars_v1_namelist_mod