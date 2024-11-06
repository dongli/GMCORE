module gomars_v1_namelist_mod

  use gomars_v1_const_mod

  implicit none

  logical :: cloudon          = .false.
  logical :: active_dust      = .false.
  logical :: co2scav          = .false.
  logical :: active_water     = .false.
  logical :: albfeed          = .false.

  namelist /gomars_v1_control/ &
    cloudon                  , &
    active_dust              , &
    co2scav                  , &
    active_water             , &
    albfeed

contains

  subroutine gomars_v1_parse_namelist(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path, status='old')
    read(10, nml=gomars_v1_control)
    close(10)

  end subroutine gomars_v1_parse_namelist

end module gomars_v1_namelist_mod