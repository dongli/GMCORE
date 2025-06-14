! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module latlon_mesh_mod

  use flogger
  use string
  use const_mod, only: pi, pi2, pi05, radius, omega, inf, deg, eps
  use sphere_geometry_mod

  implicit none

  private

  public latlon_mesh_type
  public global_mesh

  type latlon_mesh_type
    ! For nesting
    integer :: id = 0
    type(latlon_mesh_type), pointer :: parent => null()
    integer lon_hw              ! Halo width along longitude
    integer lat_hw              ! Halo width along latitude
    integer lev_hw              ! Halo width along level
    integer full_nlon, half_nlon
    integer full_nlat, half_nlat
    integer full_nlev, half_nlev
    integer full_ids, full_ide
    integer full_jds, full_jde, full_jds_no_pole, full_jde_no_pole
    integer full_kds, full_kde
    integer half_ids, half_ide
    integer half_jds, half_jde
    integer half_kds, half_kde
    integer full_ims, full_ime
    integer full_jms, full_jme
    integer full_kms, full_kme
    integer half_ims, half_ime
    integer half_jms, half_jme
    integer half_kms, half_kme
    real(8) start_lon, end_lon
    real(8) start_lat, end_lat
    real(8) dlon, dlat
    real(8), allocatable, dimension(:  ) :: full_dlev
    real(8), allocatable, dimension(:  ) :: half_dlev
    real(8), allocatable, dimension(:  ) :: full_lon
    real(8), allocatable, dimension(:  ) :: half_lon
    real(8), allocatable, dimension(:  ) :: full_lat
    real(8), allocatable, dimension(:  ) :: half_lat
    real(8), allocatable, dimension(:  ) :: full_lev
    real(8), allocatable, dimension(:  ) :: half_lev
    real(8), allocatable, dimension(:  ) :: full_cos_lon
    real(8), allocatable, dimension(:  ) :: half_cos_lon
    real(8), allocatable, dimension(:  ) :: full_sin_lon
    real(8), allocatable, dimension(:  ) :: half_sin_lon
    real(8), allocatable, dimension(:  ) :: full_cos_lat
    real(8), allocatable, dimension(:  ) :: half_cos_lat
    real(8), allocatable, dimension(:  ) :: full_sin_lat
    real(8), allocatable, dimension(:  ) :: half_sin_lat
    ! For output
    real(8), allocatable, dimension(:  ) :: full_lon_deg
    real(8), allocatable, dimension(:  ) :: half_lon_deg
    real(8), allocatable, dimension(:  ) :: full_lat_deg
    real(8), allocatable, dimension(:  ) :: half_lat_deg
    ! Area for weighting
    real(8) total_area
    real(8) area_pole_cap
    real(8), allocatable, dimension(:  ) :: area_cell
    real(8), allocatable, dimension(:  ) :: area_lon
    real(8), allocatable, dimension(:  ) :: area_lon_west
    real(8), allocatable, dimension(:  ) :: area_lon_east
    real(8), allocatable, dimension(:  ) :: area_lon_north
    real(8), allocatable, dimension(:  ) :: area_lon_south
    real(8), allocatable, dimension(:  ) :: area_lat
    real(8), allocatable, dimension(:  ) :: area_lat_west
    real(8), allocatable, dimension(:  ) :: area_lat_east
    real(8), allocatable, dimension(:  ) :: area_lat_north
    real(8), allocatable, dimension(:  ) :: area_lat_south
    real(8), allocatable, dimension(:  ) :: area_vtx
    real(8), allocatable, dimension(:,:) :: area_subcell
    ! Edge length
    real(8), allocatable, dimension(:  ) :: de_lon
    real(8), allocatable, dimension(:  ) :: de_lat
    real(8), allocatable, dimension(:  ) :: le_lat
    real(8), allocatable, dimension(:  ) :: le_lon
    ! Coriolis parameters
    real(8), allocatable, dimension(:  ) :: f_lon
    real(8), allocatable, dimension(:  ) :: f_lat
    ! Weight for constructing tangential wind
    real(8), allocatable, dimension(:,:) :: tg_wgt_lon
    real(8), allocatable, dimension(:,:) :: tg_wgt_lat
  contains
    procedure :: init_global      => latlon_mesh_init_global
    procedure :: init_from_parent => latlon_mesh_init_from_parent
    procedure :: reinit           => latlon_mesh_reinit
    procedure :: common_init      => latlon_mesh_common_init
    procedure :: has_south_pole   => latlon_mesh_has_south_pole
    procedure :: has_north_pole   => latlon_mesh_has_north_pole
    procedure :: is_south_pole    => latlon_mesh_is_south_pole
    procedure :: is_north_pole    => latlon_mesh_is_north_pole
    procedure :: is_pole          => latlon_mesh_is_pole
    procedure :: clear            => latlon_mesh_clear
    final :: latlon_mesh_final
  end type latlon_mesh_type

  type(latlon_mesh_type), target :: global_mesh

contains

  subroutine latlon_mesh_init_global(this, nlon, nlat, nlev, id, lon_hw, lat_hw, lev_hw, keep_lev)

    class(latlon_mesh_type), intent(inout) :: this
    integer, intent(in)           :: nlon
    integer, intent(in)           :: nlat
    integer, intent(in), optional :: nlev
    integer, intent(in), optional :: id
    integer, intent(in), optional :: lon_hw
    integer, intent(in), optional :: lat_hw
    integer, intent(in), optional :: lev_hw
    logical, intent(in), optional :: keep_lev

    integer nlev_opt, id_opt, lon_hw_opt, lat_hw_opt, lev_hw_opt
    logical keep_lev_opt
    real(8) dlat0
    real(16) x(3), y(3), z(3)
    integer i, j, ierr

    nlev_opt   =  1; if (present(nlev  )) nlev_opt   = nlev
    id_opt     = -1; if (present(id    )) id_opt     = id
    lon_hw_opt =  3; if (present(lon_hw)) lon_hw_opt = lon_hw
    lat_hw_opt =  3; if (present(lat_hw)) lat_hw_opt = lat_hw
    lev_hw_opt =  3; if (present(lev_hw)) lev_hw_opt = lev_hw

    call this%clear(keep_lev)

    if (mod(nlon, 2) /= 0) then
      call log_error('nlon is ' // to_str(nlon) // ', but it must be an even number!', __FILE__, __LINE__)
    end if

    this%full_nlon = nlon
    this%half_nlon = nlon
    this%full_ids  = 1
    this%full_ide  = this%full_nlon
    this%half_ids  = 1
    this%half_ide  = this%half_nlon
    this%full_nlat = nlat
    this%half_nlat = nlat - 1
    this%full_jds  = 1
    this%full_jde  = this%full_nlat
    this%half_jds  = 1
    this%half_jde  = this%half_nlat
    this%full_nlev = nlev_opt
    this%half_nlev = this%full_nlev + 1
    this%full_kds  = 1
    this%full_kde  = this%full_nlev
    this%half_kds  = 1
    this%half_kde  = this%half_nlev

    this%id        = id_opt
    this%lon_hw    = lon_hw_opt
    this%lat_hw    = lat_hw_opt
    this%lev_hw    = lev_hw_opt
    this%start_lon = 0
    this%end_lon   =  pi2
    this%start_lat = -pi05
    this%end_lat   =  pi05

    call this%common_init()

    this%dlon = (this%end_lon - this%start_lon) / this%full_nlon
    this%dlat = (this%end_lat - this%start_lat) / this%half_nlat
    do i = this%full_ims, this%full_ime
      this%full_lon(i) = this%start_lon + (i - 1.0d0) * this%dlon
    end do
    do i = this%half_ims, this%half_ime
      this%half_lon(i) = this%start_lon + (i - 0.5d0) * this%dlon
    end do
    do j = this%full_jds, this%full_jde
      this%full_lat(j) = this%start_lat + (j - 1.0d0) * this%dlat
    end do
    do j = this%half_jds, this%half_jde
      this%half_lat(j) = this%start_lat + (j - 0.5d0) * this%dlat
    end do

    ! Prevent floating point error.
    if (this%has_south_pole()) this%full_lat(j) = -pi05
    if (this%has_north_pole()) this%full_lat(j) =  pi05

    do i = this%full_ims, this%full_ime
      this%full_lon_deg(i) = this%full_lon(i) * deg
      this%full_cos_lon(i) = cos(this%full_lon(i))
      this%full_sin_lon(i) = sin(this%full_lon(i))
    end do
    do i = this%half_ims, this%half_ime
      this%half_lon_deg(i) = this%half_lon(i) * deg
      this%half_cos_lon(i) = cos(this%half_lon(i))
      this%half_sin_lon(i) = sin(this%half_lon(i))
    end do
    do j = this%full_jms, this%full_jme
      this%full_lat_deg(j) = merge(this%full_lat(j) * deg, inf, this%full_lat(j) /= inf)
      if (this%full_lat(j) >= -pi05 .and. this%full_lat(j) <= pi05) then
        this%full_cos_lat(j) = cos(this%full_lat(j))
        this%full_sin_lat(j) = sin(this%full_lat(j))
      end if
    end do
    do j = this%half_jms, this%half_jme
      this%half_lat_deg(j) = merge(this%half_lat(j) * deg, inf, this%half_lat(j) /= inf)
      if (this%half_lat(j) >= -pi05 .and. this%half_lat(j) <= pi05) then
        this%half_cos_lat(j) = cos(this%half_lat(j))
        this%half_sin_lat(j) = sin(this%half_lat(j))
        if (abs(this%half_lat(j)) < eps) this%half_sin_lat(j) = eps
      end if
    end do

    do j = this%full_jds, this%full_jde
      if (this%is_south_pole(j)) then
        this%area_cell(j) = radius**2 * this%dlon * (this%half_sin_lat(j) + 1.0d0)
        this%area_subcell(2,j) = radius**2 * 0.5d0 * this%dlon * (this%half_sin_lat(j) + 1.0d0)
      else if (this%is_north_pole(j)) then
        this%area_cell(j) = radius**2 * this%dlon * (1.0 - this%half_sin_lat(j-1))
        this%area_subcell(1,j) = radius**2 * 0.5d0 * this%dlon * (1.0d0 - this%half_sin_lat(j-1))
      else
        this%area_cell(j) = radius**2 * this%dlon * (this%half_sin_lat(j) - this%half_sin_lat(j-1))
        this%area_subcell(1,j) = radius**2 * 0.5d0 * this%dlon * (this%full_sin_lat(j) - this%half_sin_lat(j-1))
        this%area_subcell(2,j) = radius**2 * 0.5d0 * this%dlon * (this%half_sin_lat(j) - this%full_sin_lat(j))
        !
        !           1,j
        !           /|
        !          / |
        !         /  |
        !        /   |
        !    1,j \   |
        !         \  |
        !          \ |
        !           \|
        !          1,j-1
        !
        call lonlat2xyz(radius, this%full_lon(1), this%full_lat(j  ), x(1), y(1), z(1))
        call lonlat2xyz(radius, this%half_lon(1), this%half_lat(j-1), x(2), y(2), z(2))
        call lonlat2xyz(radius, this%half_lon(1), this%half_lat(j  ), x(3), y(3), z(3))
        this%area_lon_west(j) = spherical_area(radius, x, y, z, ierr)
        if (ierr /= 0) then
          call log_error(sphere_geometry_error_message(ierr), __FILE__, __LINE__)
        end if
        this%area_lon_east(j) = this%area_lon_west(j)
        this%area_lon(j) = this%area_lon_west(j) + this%area_lon_east(j)
        !
        !          1,j
        !           /\
        !          /  \
        !         /    \
        !        /______\
        !    1,j          2,j
        !
        call lonlat2xyz(radius, this%half_lon(1), this%half_lat(j  ), x(1), y(1), z(1))
        call lonlat2xyz(radius, this%full_lon(1), this%full_lat(j  ), x(2), y(2), z(2))
        call lonlat2xyz(radius, this%full_lon(2), this%full_lat(j  ), x(3), y(3), z(3))
        this%area_lon_north(j) = spherical_area_with_last_small_arc(radius, x, y, z, ierr)
        if (ierr /= 0) then
          call log_error(sphere_geometry_error_message(ierr), __FILE__, __LINE__)
        end if
        !
        !    1,j          2,j
        !        --------
        !        \      /
        !         \    /
        !          \  /
        !           \/
        !         1,j-1
        !
        call lonlat2xyz(radius, this%half_lon(1), this%half_lat(j-1), x(1), y(1), z(1))
        call lonlat2xyz(radius, this%full_lon(2), this%full_lat(j  ), x(2), y(2), z(2))
        call lonlat2xyz(radius, this%full_lon(1), this%full_lat(j  ), x(3), y(3), z(3))
        this%area_lon_south(j) = spherical_area_with_last_small_arc(radius, x, y, z, ierr)
        if (ierr /= 0) then
          call log_error(sphere_geometry_error_message(ierr), __FILE__, __LINE__)
        end if
      end if
    end do

    ! Set cell areas exceed the Poles.
    if (this%has_south_pole()) then
      do j = this%full_jms, this%full_jds - 1
        this%area_cell(j) = this%area_cell(2 - j)
      end do
    end if
    if (this%has_north_pole()) then
      do j = this%full_jde + 1, this%full_jme
        this%area_cell(j) = this%area_cell(2 * this%full_jde - j)
      end do
    end if

    do j = this%half_jds, this%half_jde
      !
      !          2,j+1
      !           /|
      !          / |
      !         /  |
      !        /   |
      !    1,j \   |
      !         \  |
      !          \ |
      !           \|
      !           2,j
      !
      call lonlat2xyz(radius, this%half_lon(1), this%half_lat(j  ), x(1), y(1), z(1))
      call lonlat2xyz(radius, this%full_lon(2), this%full_lat(j  ), x(2), y(2), z(2))
      call lonlat2xyz(radius, this%full_lon(2), this%full_lat(j+1), x(3), y(3), z(3))
      this%area_lat_west(j) = spherical_area(radius, x, y, z, ierr)
      if (ierr /= 0) then
        call log_error(sphere_geometry_error_message(ierr), __FILE__, __LINE__)
      end if
      this%area_lat_east(j) = this%area_lat_west(j)
      !
      !         2,j+1
      !           /\
      !          /  \
      !         /    \
      !        /______\
      !    1,j          2,j
      !
      call lonlat2xyz(radius, this%full_lon(2), this%full_lat(j+1), x(1), y(1), z(1))
      call lonlat2xyz(radius, this%half_lon(1), this%half_lat(j  ), x(2), y(2), z(2))
      call lonlat2xyz(radius, this%half_lon(2), this%half_lat(j  ), x(3), y(3), z(3))
      this%area_lat_north(j) = spherical_area_with_last_small_arc(radius, x, y, z, ierr)
      if (ierr /= 0) then
        call log_error(sphere_geometry_error_message(ierr), __FILE__, __LINE__)
      end if
      !
      !    1,j          2,j
      !        --------
      !        \      /
      !         \    /
      !          \  /
      !           \/
      !          2,j
      !
      call lonlat2xyz(radius, this%full_lon(2), this%full_lat(j), x(1), y(1), z(1))
      call lonlat2xyz(radius, this%half_lon(2), this%half_lat(j), x(2), y(2), z(2))
      call lonlat2xyz(radius, this%half_lon(1), this%half_lat(j), x(3), y(3), z(3))
      this%area_lat_south(j) = spherical_area_with_last_small_arc(radius, x, y, z, ierr)
      if (ierr /= 0) then
        call log_error(sphere_geometry_error_message(ierr), __FILE__, __LINE__)
      end if
      this%area_lat(j) = this%area_lat_north(j) + this%area_lat_south(j)
    end do

    do j = this%half_jds, this%half_jde
      if (this%is_south_pole(j)) then
        this%area_vtx(j) = this%area_lat_west(j) + this%area_lat_east(j) + this%area_lon_south(j+1)
      else if (this%is_north_pole(j+1)) then
        this%area_vtx(j) = this%area_lat_west(j) + this%area_lat_east(j) + this%area_lon_north(j)
      else
        this%area_vtx(j) = this%area_lat_west(j) + this%area_lat_east(j) + this%area_lon_south(j+1) + this%area_lon_north(j)
      end if
    end do

    this%area_pole_cap = 2 * pi * radius**2 * (1 - cos(0.5d0 * this%dlat))

    do j = this%full_jds_no_pole, this%full_jde_no_pole
      this%de_lon(j) = radius * this%full_cos_lat(j) * this%dlon
      this%le_lon(j) = 2.0d0 * this%area_lon(j) / this%de_lon(j)
    end do

    do j = this%half_jds, this%half_jde
      this%le_lat(j) = radius * this%half_cos_lat(j) * this%dlon
      this%de_lat(j) = 2.0d0 * this%area_lat(j) / this%le_lat(j)
    end do

    do j = this%full_jds, this%full_jde
      this%f_lon(j) = 2 * omega * this%full_sin_lat(j)
    end do
    do j = this%half_jds, this%half_jde
      this%f_lat(j) = 2 * omega * this%half_sin_lat(j)
    end do

    do j = this%full_jds_no_pole, this%full_jde_no_pole
      this%tg_wgt_lon(1,j) = this%le_lat(j-1) / this%de_lon(j) * 0.25d0
      this%tg_wgt_lon(2,j) = this%le_lat(j  ) / this%de_lon(j) * 0.25d0
    end do

    do j = this%half_jds, this%half_jde
      this%tg_wgt_lat(1,j) = this%le_lon(j  ) / this%de_lat(j) * 0.25d0
      this%tg_wgt_lat(2,j) = this%le_lon(j+1) / this%de_lat(j) * 0.25d0
    end do

  end subroutine latlon_mesh_init_global

  subroutine latlon_mesh_init_from_parent(this, parent, id, ids, ide, jds, jde, keep_lev)

    class(latlon_mesh_type), intent(inout) :: this
    class(latlon_mesh_type), intent(in), target :: parent
    integer, intent(in) :: id
    integer, intent(in) :: ids
    integer, intent(in) :: ide
    integer, intent(in) :: jds
    integer, intent(in) :: jde
    logical, intent(in), optional :: keep_lev

    integer i, j

    call this%clear(keep_lev)

    this%parent => parent

    this%full_nlon = ide - ids + 1
    this%half_nlon = this%full_nlon
    this%full_ids  = ids
    this%full_ide  = ide
    this%half_ids  = ids
    this%half_ide  = ide
    this%full_nlat = jde - jds + 1
    this%full_jds  = jds
    this%full_jde  = jde
    this%half_jds  = jds
    this%half_jde  = merge(jde - 1, jde, this%has_north_pole())
    this%half_nlat = this%half_jde - this%half_jds + 1

    this%full_nlev = parent%full_nlev
    this%half_nlev = parent%half_nlev
    this%full_kds  = parent%full_kds
    this%full_kde  = parent%full_kde
    this%half_kds  = parent%half_kds
    this%half_kde  = parent%half_kde

    this%id        = id
    this%lon_hw    = parent%lon_hw
    this%lat_hw    = parent%lat_hw
    this%lev_hw    = parent%lev_hw
    this%start_lon = parent%full_lon(ids)
    this%end_lon   = parent%full_lon(ide) + parent%dlon
    this%start_lat = merge(parent%full_lat(jds), -pi05, .not. this%has_south_pole())
    this%end_lat   = merge(parent%full_lat(jde),  pi05, .not. this%has_north_pole())

    call this%common_init()

    this%full_dlev = parent%full_dlev
    this%half_dlev = parent%half_dlev

    this%dlon = parent%dlon
    do i = this%full_ims, this%full_ime
      this%full_lon(i)     = parent%full_lon(i)
      this%half_lon(i)     = parent%half_lon(i)
      this%full_lon_deg(i) = parent%full_lon_deg(i)
      this%half_lon_deg(i) = parent%half_lon_deg(i)
      this%full_sin_lon(i) = parent%full_sin_lon(i)
      this%half_sin_lon(i) = parent%half_sin_lon(i)
      this%full_cos_lon(i) = parent%full_cos_lon(i)
      this%half_cos_lon(i) = parent%half_cos_lon(i)
    end do

    this%dlat = parent%dlat
    do j = this%full_jms, this%full_jme
      this%full_lat      (j) = parent%full_lat      (j)
      this%full_lat_deg  (j) = parent%full_lat_deg  (j)
      this%full_sin_lat  (j) = parent%full_sin_lat  (j)
      this%full_cos_lat  (j) = parent%full_cos_lat  (j)
      this%area_cell     (j) = parent%area_cell     (j)
      this%area_subcell(:,j) = parent%area_subcell(:,j)
      this%area_lon_west (j) = parent%area_lon_west (j)
      this%area_lon_east (j) = parent%area_lon_east (j)
      this%area_lon_north(j) = parent%area_lon_north(j)
      this%area_lon_south(j) = parent%area_lon_south(j)
      this%area_lon      (j) = parent%area_lon      (j)
      this%le_lon        (j) = parent%le_lon        (j)
      this%de_lon        (j) = parent%de_lon        (j)
      this%f_lon         (j) = parent%f_lon         (j)
      this%tg_wgt_lon  (:,j) = parent%tg_wgt_lon  (:,j)
    end do
    do j = this%half_jms, this%half_jme
      this%half_lat      (j) = parent%half_lat      (j)
      this%half_lat_deg  (j) = parent%half_lat_deg  (j)
      this%half_sin_lat  (j) = parent%half_sin_lat  (j)
      this%half_cos_lat  (j) = parent%half_cos_lat  (j)
      this%area_vtx      (j) = parent%area_vtx      (j)
      this%area_lat_west (j) = parent%area_lat_west (j)
      this%area_lat_east (j) = parent%area_lat_east (j)
      this%area_lat_north(j) = parent%area_lat_north(j)
      this%area_lat_south(j) = parent%area_lat_south(j)
      this%area_lat      (j) = parent%area_lat      (j)
      this%le_lat        (j) = parent%le_lat        (j)
      this%de_lat        (j) = parent%de_lat        (j)
      this%f_lat         (j) = parent%f_lat         (j)
      this%tg_wgt_lat  (:,j) = parent%tg_wgt_lat  (:,j)
    end do
    this%area_pole_cap = parent%area_pole_cap

    this%full_lev = parent%full_lev
    this%half_lev = parent%half_lev

  end subroutine latlon_mesh_init_from_parent

  subroutine latlon_mesh_reinit(this, lon_hw)

    class(latlon_mesh_type), intent(inout) :: this
    integer, intent(in), optional :: lon_hw

    integer nlon, nlat, nlev, lat_hw, lev_hw

    nlon = global_mesh%full_nlon
    nlat = global_mesh%full_nlat
    nlev = global_mesh%full_nlev
    lat_hw = global_mesh%lat_hw
    lev_hw = global_mesh%lev_hw

    if (associated(this%parent)) then
      call this%init_from_parent(this%parent, this%id, this%full_ids, this%full_ide, this%full_jds, this%full_jde, keep_lev=.true.)
    else if (present(lon_hw)) then
      call this%init_global(nlon, nlat, nlev, 0, lon_hw, lat_hw, lev_hw, keep_lev=.true.)
    else
      call log_error('Logical error!', __FILE__, __LINE__)
    end if

  end subroutine latlon_mesh_reinit

  subroutine latlon_mesh_common_init(this)

    class(latlon_mesh_type), intent(inout) :: this

    this%total_area = radius**2 * (this%end_lon - this%start_lon) * (sin(this%end_lat) - sin(this%start_lat))

    this%full_jds_no_pole = merge(this%full_jds + 1, this%full_jds, this%has_south_pole())
    this%full_jde_no_pole = merge(this%full_jde - 1, this%full_jde, this%has_north_pole())

    ! Use maximum lon_hw in this process and its south and north neighbors.
    this%full_ims = this%full_ids - this%lon_hw
    this%full_ime = this%full_ide + this%lon_hw
    this%full_jms = this%full_jds - this%lat_hw
    this%full_jme = this%full_jde + this%lat_hw
    this%half_ims = this%half_ids - this%lon_hw
    this%half_ime = this%half_ide + this%lon_hw
    this%half_jms = this%half_jds - this%lat_hw
    this%half_jme = this%half_jde + this%lat_hw
    this%full_kms = this%full_kds - this%lev_hw
    this%full_kme = this%full_kde + this%lev_hw
    this%half_kms = this%half_kds - this%lev_hw
    this%half_kme = this%half_kde + this%lev_hw

    if (.not. allocated(this%full_lev)) then
      allocate(this%full_dlev        (this%full_kms:this%full_kme)); this%full_dlev           = 0
      allocate(this%half_dlev        (this%half_kms:this%half_kme)); this%half_dlev           = 0
      allocate(this%full_lev         (this%full_kms:this%full_kme)); this%full_lev            = inf
      allocate(this%half_lev         (this%half_kms:this%half_kme)); this%half_lev            = inf
    end if

    allocate(this%full_lon           (this%full_ims:this%full_ime)); this%full_lon            = inf
    allocate(this%half_lon           (this%half_ims:this%half_ime)); this%half_lon            = inf
    allocate(this%full_lat           (this%full_jms:this%full_jme)); this%full_lat            = inf
    allocate(this%half_lat           (this%half_jms:this%half_jme)); this%half_lat            = inf
    allocate(this%full_cos_lon       (this%full_ims:this%full_ime)); this%full_cos_lon        = inf
    allocate(this%half_cos_lon       (this%half_ims:this%half_ime)); this%half_cos_lon        = inf
    allocate(this%full_sin_lon       (this%full_ims:this%full_ime)); this%full_sin_lon        = inf
    allocate(this%half_sin_lon       (this%half_ims:this%half_ime)); this%half_sin_lon        = inf
    allocate(this%full_cos_lat       (this%full_jms:this%full_jme)); this%full_cos_lat        = inf
    allocate(this%half_cos_lat       (this%half_jms:this%half_jme)); this%half_cos_lat        = inf
    allocate(this%full_sin_lat       (this%full_jms:this%full_jme)); this%full_sin_lat        = inf
    allocate(this%half_sin_lat       (this%half_jms:this%half_jme)); this%half_sin_lat        = inf
    allocate(this%full_lon_deg       (this%full_ims:this%full_ime)); this%full_lon_deg        = inf
    allocate(this%half_lon_deg       (this%half_ims:this%half_ime)); this%half_lon_deg        = inf
    allocate(this%full_lat_deg       (this%full_jms:this%full_jme)); this%full_lat_deg        = inf
    allocate(this%half_lat_deg       (this%half_jms:this%half_jme)); this%half_lat_deg        = inf
    allocate(this%area_cell          (this%full_jms:this%full_jme)); this%area_cell           = 0
    allocate(this%area_lon           (this%full_jms:this%full_jme)); this%area_lon            = 0
    allocate(this%area_lon_west      (this%full_jms:this%full_jme)); this%area_lon_west       = 0
    allocate(this%area_lon_east      (this%full_jms:this%full_jme)); this%area_lon_east       = 0
    allocate(this%area_lon_north     (this%full_jms:this%full_jme)); this%area_lon_north      = 0
    allocate(this%area_lon_south     (this%full_jms:this%full_jme)); this%area_lon_south      = 0
    allocate(this%area_lat           (this%half_jms:this%half_jme)); this%area_lat            = 0
    allocate(this%area_lat_west      (this%half_jms:this%half_jme)); this%area_lat_west       = 0
    allocate(this%area_lat_east      (this%half_jms:this%half_jme)); this%area_lat_east       = 0
    allocate(this%area_lat_north     (this%half_jms:this%half_jme)); this%area_lat_north      = 0
    allocate(this%area_lat_south     (this%half_jms:this%half_jme)); this%area_lat_south      = 0
    allocate(this%area_vtx           (this%half_jms:this%half_jme)); this%area_vtx            = 0
    allocate(this%area_subcell     (2,this%full_jms:this%full_jme)); this%area_subcell        = 0
    allocate(this%de_lon             (this%full_jms:this%full_jme)); this%de_lon              = 0
    allocate(this%de_lat             (this%half_jms:this%half_jme)); this%de_lat              = 0
    allocate(this%le_lat             (this%half_jms:this%half_jme)); this%le_lat              = 0
    allocate(this%le_lon             (this%full_jms:this%full_jme)); this%le_lon              = 0
    allocate(this%f_lon              (this%full_jms:this%full_jme)); this%f_lon               = inf
    allocate(this%f_lat              (this%half_jms:this%half_jme)); this%f_lat               = inf
    allocate(this%tg_wgt_lon       (2,this%full_jms:this%full_jme)); this%tg_wgt_lon          = inf
    allocate(this%tg_wgt_lat       (2,this%half_jms:this%half_jme)); this%tg_wgt_lat          = inf

  end subroutine latlon_mesh_common_init

  pure logical function latlon_mesh_has_south_pole(this) result(res)

    class(latlon_mesh_type), intent(in) :: this

    res = this%full_jds == 1

  end function latlon_mesh_has_south_pole

  pure logical function latlon_mesh_has_north_pole(this) result(res)

    class(latlon_mesh_type), intent(in) :: this

    res = this%full_jde == global_mesh%full_nlat

  end function latlon_mesh_has_north_pole

  pure logical function latlon_mesh_is_south_pole(this, j) result(res)

    class(latlon_mesh_type), intent(in) :: this
    integer, intent(in) :: j

    ! FIXME: has_south_pole should be removed.
    res = j == 1

  end function latlon_mesh_is_south_pole

  pure logical function latlon_mesh_is_north_pole(this, j) result(res)

    class(latlon_mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j == global_mesh%full_nlat

  end function latlon_mesh_is_north_pole

  pure logical function latlon_mesh_is_pole(this, j) result(res)

    class(latlon_mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = this%is_south_pole(j) .or. this%is_north_pole(j)

  end function latlon_mesh_is_pole

  subroutine latlon_mesh_clear(this, keep_lev)

    class(latlon_mesh_type), intent(inout) :: this
    logical, intent(in), optional :: keep_lev

    logical keep_lev_opt

    if (present(keep_lev)) then
      keep_lev_opt = keep_lev
    else
      keep_lev_opt = .false.
    end if

    if (.not. keep_lev_opt) then
      if (allocated(this%full_dlev   )) deallocate(this%full_dlev     )
      if (allocated(this%half_dlev   )) deallocate(this%half_dlev     )
      if (allocated(this%full_lev    )) deallocate(this%full_lev      )
      if (allocated(this%half_lev    )) deallocate(this%half_lev      )
    end if

    if (allocated(this%full_lon      )) deallocate(this%full_lon      )
    if (allocated(this%full_lat      )) deallocate(this%full_lat      )
    if (allocated(this%half_lon      )) deallocate(this%half_lon      )
    if (allocated(this%half_lat      )) deallocate(this%half_lat      )
    if (allocated(this%full_cos_lon  )) deallocate(this%full_cos_lon  )
    if (allocated(this%half_cos_lon  )) deallocate(this%half_cos_lon  )
    if (allocated(this%full_sin_lon  )) deallocate(this%full_sin_lon  )
    if (allocated(this%half_sin_lon  )) deallocate(this%half_sin_lon  )
    if (allocated(this%full_cos_lat  )) deallocate(this%full_cos_lat  )
    if (allocated(this%half_cos_lat  )) deallocate(this%half_cos_lat  )
    if (allocated(this%full_sin_lat  )) deallocate(this%full_sin_lat  )
    if (allocated(this%half_sin_lat  )) deallocate(this%half_sin_lat  )
    if (allocated(this%full_lon_deg  )) deallocate(this%full_lon_deg  )
    if (allocated(this%half_lon_deg  )) deallocate(this%half_lon_deg  )
    if (allocated(this%full_lat_deg  )) deallocate(this%full_lat_deg  )
    if (allocated(this%half_lat_deg  )) deallocate(this%half_lat_deg  )
    if (allocated(this%area_cell     )) deallocate(this%area_cell     )
    if (allocated(this%area_lon      )) deallocate(this%area_lon      )
    if (allocated(this%area_lon_west )) deallocate(this%area_lon_west )
    if (allocated(this%area_lon_east )) deallocate(this%area_lon_east )
    if (allocated(this%area_lon_north)) deallocate(this%area_lon_north)
    if (allocated(this%area_lon_south)) deallocate(this%area_lon_south)
    if (allocated(this%area_lat      )) deallocate(this%area_lat      )
    if (allocated(this%area_lat_west )) deallocate(this%area_lat_west )
    if (allocated(this%area_lat_east )) deallocate(this%area_lat_east )
    if (allocated(this%area_lat_north)) deallocate(this%area_lat_north)
    if (allocated(this%area_lat_south)) deallocate(this%area_lat_south)
    if (allocated(this%area_vtx      )) deallocate(this%area_vtx      )
    if (allocated(this%area_subcell  )) deallocate(this%area_subcell  )
    if (allocated(this%de_lon        )) deallocate(this%de_lon        )
    if (allocated(this%de_lat        )) deallocate(this%de_lat        )
    if (allocated(this%le_lat        )) deallocate(this%le_lat        )
    if (allocated(this%le_lon        )) deallocate(this%le_lon        )
    if (allocated(this%f_lon         )) deallocate(this%f_lon         )
    if (allocated(this%f_lat         )) deallocate(this%f_lat         )
    if (allocated(this%tg_wgt_lon    )) deallocate(this%tg_wgt_lon    )
    if (allocated(this%tg_wgt_lat    )) deallocate(this%tg_wgt_lat    )

  end subroutine latlon_mesh_clear

  subroutine latlon_mesh_final(this)

    type(latlon_mesh_type), intent(inout) :: this

    call this%clear()

  end subroutine latlon_mesh_final

end module latlon_mesh_mod
