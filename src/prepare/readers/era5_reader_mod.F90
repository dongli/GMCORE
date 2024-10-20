module era5_reader_mod

  use fiona
  use flogger
  use const_mod
  use namelist_mod
  use time_mod
  use formula_mod
  use process_mod

  implicit none

  integer era5_nlon
  integer era5_nlat
  integer era5_nlev

  real(r8), allocatable, dimension(:    ) :: era5_lon
  real(r8), allocatable, dimension(:    ) :: era5_lat
  real(r8), allocatable, dimension(:    ) :: era5_lev
  real(r8), allocatable, dimension(:,:,:) :: era5_u
  real(r8), allocatable, dimension(:,:,:) :: era5_v
  real(r8), allocatable, dimension(:,:,:) :: era5_t
  real(r8), allocatable, dimension(:,:,:) :: era5_z
  real(r8), allocatable, dimension(:,:,:) :: era5_qv
  real(r8), allocatable, dimension(:,:,:) :: era5_qc
  real(r8), allocatable, dimension(:,:,:) :: era5_qi
  real(r8), allocatable, dimension(:,:,:) :: era5_qr
  real(r8), allocatable, dimension(:,:,:) :: era5_qs
  real(r8), allocatable, dimension(:,:  ) :: era5_ps
  real(r8), allocatable, dimension(:,:  ) :: era5_zs

contains

  subroutine era5_reader_run(bkg_file, min_lon, max_lon, min_lat_in, max_lat_in)

    character(*), intent(in) :: bkg_file
    real(r8), intent(in) :: min_lon
    real(r8), intent(in) :: max_lon
    real(r8), intent(in) :: min_lat_in
    real(r8), intent(in) :: max_lat_in

    integer i, j, k, k0
    character(50) time_units
    real(r8) min_lat, max_lat, time_value

    call era5_reader_final()

    if (proc%is_root()) call log_notice('Use ERA5 ' // trim(bkg_file) // ' as background.')

    min_lat = max(min_lat_in, -89.75)
    max_lat = min(max_lat_in,  89.75)

    call fiona_open_dataset('era5', file_path=bkg_file, mpi_comm=proc%comm_io, ngroups=input_ngroups)
    call fiona_set_dim('era5', 'longitude', span=[0, 360], cyclic=.true.)
    call fiona_set_dim('era5', 'latitude', span=[90, -90], flip=.true.)
    call fiona_get_dim('era5', 'level', size=era5_nlev)
    allocate(era5_lev(era5_nlev))
    call fiona_start_input('era5')
    call fiona_input('era5', 'time', time_value)
    call fiona_get_att('era5', 'time', 'units', time_units)
    call time_fast_forward(time_value, time_units, change_end_time=.true.)
    call fiona_input_range('era5', 'longitude', era5_lon, coord_range=[min_lon, max_lon]); era5_nlon = size(era5_lon)
    call fiona_input_range('era5', 'latitude' , era5_lat, coord_range=[min_lat, max_lat]); era5_nlat = size(era5_lat)
    call fiona_input      ('era5', 'level'    , era5_lev)
    call fiona_input_range('era5', 'u'        , era5_u  , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_input_range('era5', 'v'        , era5_v  , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_input_range('era5', 't'        , era5_t  , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_input_range('era5', 'z'        , era5_z  , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_input_range('era5', 'q'        , era5_qv , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_input_range('era5', 'sp'       , era5_ps , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_input_range('era5', 'zs'       , era5_zs , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
  if (fiona_has_var('era5', 'clwc')) &
    call fiona_input_range('era5', 'clwc'     , era5_qc , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
  if (fiona_has_var('era5', 'ciwc')) &
    call fiona_input_range('era5', 'ciwc'     , era5_qi , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
  if (fiona_has_var('era5', 'crwc')) &
    call fiona_input_range('era5', 'crwc'     , era5_qr , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
  if (fiona_has_var('era5', 'cswc')) &
    call fiona_input_range('era5', 'cswc'     , era5_qs , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_end_input('era5')

    do j = 1, era5_nlat
      if (abs(abs(era5_lat(j)) - 90) < 1.0e-10) then
        era5_u(:,j,:) = 0
        era5_v(:,j,:) = 0
      end if
    end do

    ! Change units.
    era5_lev = era5_lev * 100.0_r8
    era5_z   = era5_z  / g
    era5_zs  = era5_zs / g

    if (.not. allocated(era5_qc)) then
      allocate(era5_qc(era5_nlon,era5_nlat,era5_nlev))
      era5_qc = 0
    end if
    if (.not. allocated(era5_qi)) then
      allocate(era5_qi(era5_nlon,era5_nlat,era5_nlev))
      era5_qi = 0
    end if
    if (.not. allocated(era5_qr)) then
      allocate(era5_qr(era5_nlon,era5_nlat,era5_nlev))
      era5_qr = 0
    end if
    if (.not. allocated(era5_qs)) then
      allocate(era5_qs(era5_nlon,era5_nlat,era5_nlev))
      era5_qs = 0
    end if

  end subroutine era5_reader_run

  subroutine era5_reader_final()

    if (allocated(era5_lon)) deallocate(era5_lon)
    if (allocated(era5_lat)) deallocate(era5_lat)
    if (allocated(era5_lev)) deallocate(era5_lev)
    if (allocated(era5_u  )) deallocate(era5_u  )
    if (allocated(era5_v  )) deallocate(era5_v  )
    if (allocated(era5_t  )) deallocate(era5_t  )
    if (allocated(era5_z  )) deallocate(era5_z  )
    if (allocated(era5_qv )) deallocate(era5_qv )
    if (allocated(era5_qc )) deallocate(era5_qc )
    if (allocated(era5_qi )) deallocate(era5_qi )
    if (allocated(era5_qr )) deallocate(era5_qr )
    if (allocated(era5_qs )) deallocate(era5_qs )
    if (allocated(era5_ps )) deallocate(era5_ps )
    if (allocated(era5_zs )) deallocate(era5_zs )

  end subroutine era5_reader_final

end module era5_reader_mod
