! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================
! Description:
!
!   This is module provides field types which encapsulates meta data with array.
!   By doing so, users can send the field objects into their functions without
!   the need to specify the array index ranges.
!
! Authors:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
! ==============================================================================

module latlon_field_types_mod

  use mpi
  use container
  use fiona
  use flogger
  use const_mod
  use latlon_mesh_mod
  use latlon_halo_mod
  use latlon_parallel_types_mod
  use latlon_parallel_global_mod

  implicit none

  private

  public latlon_mesh_type
  public latlon_halo_type
  public latlon_field2d_type
  public latlon_field3d_type
  public latlon_field4d_type
  public append_field

  integer, public, parameter :: field_name_len      = 32
  integer, public, parameter :: field_long_name_len = 128
  integer, public, parameter :: field_units_len     = 32
  integer, public, parameter :: field_loc_len       = 10

  integer :: field_unique_id = 0

  type latlon_field_meta_type
    character(field_name_len     ) :: name      = 'N/A'
    character(field_long_name_len) :: long_name = 'N/A'
    character(field_units_len    ) :: units     = 'N/A'
    character(field_loc_len      ) :: loc       = 'N/A'
    character(5)                   :: output    = 'N/A'
    integer :: id               = 0
    integer :: nlon             = 0
    integer :: nlat             = 0
    integer :: nlev             = 0
    logical :: full_lon         = .true.
    logical :: full_lat         = .true.
    logical :: full_lev         = .true.
    logical :: initialized      = .false.
    logical :: linked           = .false.
    logical :: restart          = .false.
    logical :: initial          = .false.
    logical :: halo_cross_pole  = .false.
    logical :: halo_diagonal    = .false.
    type(latlon_mesh_type), pointer :: mesh     => null()
    type(latlon_halo_type), pointer :: halo(:)  => null()
  contains
    procedure latlon_field_meta_init
    procedure latlon_field_meta_clear
    procedure latlon_field_meta_link
  end type

  type, extends(latlon_field_meta_type) :: latlon_field2d_type
    integer, pointer :: send_req(:) => null()
    integer, pointer :: recv_req(:) => null()
    real(r8), contiguous, pointer :: d(:,:) => null()
  contains
    procedure :: init   => latlon_field2d_init
    procedure :: clear  => latlon_field2d_clear
    procedure :: copy   => latlon_field2d_copy
    procedure :: sum    => latlon_field2d_sum
    procedure :: absmax => latlon_field2d_absmax
    procedure, private :: latlon_field2d_link_2d
    procedure, private :: latlon_field2d_link_3d
    generic :: link => latlon_field2d_link_2d, latlon_field2d_link_3d
    final latlon_field2d_final
  end type latlon_field2d_type

  type, extends(latlon_field_meta_type) :: latlon_field3d_type
    integer, pointer :: send_req(:) => null()
    integer, pointer :: recv_req(:) => null()
    real(r8), contiguous, pointer :: d(:,:,:) => null()
  contains
    procedure :: init  => latlon_field3d_init
    procedure :: clear => latlon_field3d_clear
    procedure, private :: latlon_field3d_copy_2d
    procedure, private :: latlon_field3d_copy_3d
    procedure, private :: latlon_field3d_copy_4d
    generic :: copy    => latlon_field3d_copy_2d, latlon_field3d_copy_3d, latlon_field3d_copy_4d
    procedure, private :: latlon_field3d_link_3d
    procedure, private :: latlon_field3d_link_4d
    generic :: link    => latlon_field3d_link_3d, latlon_field3d_link_4d
    procedure :: sum   => latlon_field3d_sum
    procedure :: add   => latlon_field3d_add
    procedure :: mul   => latlon_field3d_mul
    procedure :: div   => latlon_field3d_div
    procedure :: min   => latlon_field3d_min
    procedure :: max   => latlon_field3d_max
    final latlon_field3d_final
  end type latlon_field3d_type

  type, extends(latlon_field_meta_type) :: latlon_field4d_type
    character(30) :: dim4_name = ''
    character(30), allocatable :: var4_names(:)
    integer :: dim4_size = 0
    integer, pointer :: send_req(:,:) => null()
    integer, pointer :: recv_req(:,:) => null()
    real(r8), contiguous, pointer :: d(:,:,:,:) => null()
  contains
    procedure :: init  => latlon_field4d_init
    procedure :: clear => latlon_field4d_clear
    procedure, private :: latlon_field4d_copy_3d
    generic :: copy    => latlon_field4d_copy_3d
    procedure :: link  => latlon_field4d_link
    procedure :: mul   => latlon_field4d_mul
    procedure :: div   => latlon_field4d_div
    final latlon_field4d_final
  end type latlon_field4d_type

  interface append_field
    module procedure append_field2d
    module procedure append_field3d
    module procedure append_field4d
  end interface append_field

contains

  subroutine latlon_field_meta_init(this, name, long_name, units, loc, mesh, halo, &
                                    halo_cross_pole, halo_diagonal, output, restart)

    class(latlon_field_meta_type), intent(inout) :: this
    character(*), intent(in) :: name
    character(*), intent(in) :: long_name
    character(*), intent(in) :: units
    character(*), intent(in) :: loc
    type(latlon_mesh_type), intent(in), target :: mesh
    type(latlon_halo_type), intent(in), optional, target :: halo(:)
    logical, intent(in), optional :: halo_cross_pole
    logical, intent(in), optional :: halo_diagonal
    character(*), intent(in), optional :: output
    logical, intent(in), optional :: restart

    this%name      = name
    this%long_name = long_name
    this%units     = units
    this%loc       = loc
    this%mesh      => mesh
    if (present(halo           )) this%halo            => halo
    if (present(halo_cross_pole)) this%halo_cross_pole = halo_cross_pole
    if (present(halo_diagonal  )) this%halo_diagonal   = halo_diagonal
    if (present(output         )) this%output          = output
    if (present(restart        )) this%restart         = restart

    ! This id is used to avoid MPI message tag conflict.
    field_unique_id = field_unique_id + 1
    this%id         = field_unique_id

  end subroutine latlon_field_meta_init

  subroutine latlon_field_meta_clear(this)

    class(latlon_field_meta_type), intent(inout) :: this

    this%name            = 'N/A'
    this%long_name       = 'N/A'
    this%units           = 'N/A'
    this%loc             = 'N/A'
    this%output          = 'N/A'
    this%mesh            => null()
    this%halo            => null()
    this%halo_cross_pole = .false.
    this%halo_diagonal   = .false.
    this%restart         = .false.
    this%initial         = .false.

  end subroutine latlon_field_meta_clear

  subroutine latlon_field_meta_link(this, other)

    class(latlon_field_meta_type), intent(inout) :: this
    class(latlon_field_meta_type), intent(in) :: other

    if (this%loc == 'N/A') this%loc = other%loc
#ifdef NDEBUG
    if (this%loc == 'cell' .and. (other%loc /= 'cell' .and. other%loc /= 'lev')) then
      call log_error('Location ' // trim(this%loc) // ' and ' // trim(other%loc) // ' does not match!', __FILE__, __LINE__, pid=proc%id_model)
    end if
#endif
    this%mesh            => other%mesh
    this%halo            => other%halo
    this%halo_cross_pole = other%halo_cross_pole
    this%halo_diagonal   = other%halo_diagonal

  end subroutine latlon_field_meta_link

  subroutine latlon_field2d_init(this, name, long_name, units, loc, mesh, halo, &
                                 halo_cross_pole, halo_diagonal, link, output, restart)

    class(latlon_field2d_type), intent(inout) :: this
    character(*), intent(in) :: name
    character(*), intent(in) :: long_name
    character(*), intent(in) :: units
    character(*), intent(in) :: loc
    type(latlon_mesh_type), intent(in), target :: mesh
    type(latlon_halo_type), intent(in), optional, target :: halo(:)
    logical, intent(in), optional :: halo_cross_pole
    logical, intent(in), optional :: halo_diagonal
    type(latlon_field2d_type), intent(in), optional, target :: link
    character(*), intent(in), optional :: output
    logical, intent(in), optional :: restart

    call this%clear()

    call this%latlon_field_meta_init(name, long_name, units, loc, mesh, halo, &
                                     halo_cross_pole, halo_diagonal, output, restart)

    if (present(link)) then
      this%d => link%d
      this%linked = .true.
    else
      select case (loc)
      case ('cell')
        this%full_lon = .true. ; this%full_lat = .true.
        allocate(this%d(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme))
      case ('lon')
        this%full_lon = .false.; this%full_lat = .true.
        allocate(this%d(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme))
      case ('lat')
        this%full_lon = .true. ; this%full_lat = .false.
        allocate(this%d(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme))
      case ('vtx')
        this%full_lon = .false.; this%full_lat = .false.
        allocate(this%d(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme))
      end select
    end if

    this%d           = 0
    this%nlon        = merge(mesh%full_nlon, mesh%half_nlon, this%full_lon)
    this%nlat        = merge(mesh%full_nlat, mesh%half_nlat, this%full_lat)
    this%nlev        = 1
    this%initialized = .true.

    allocate(this%send_req(8)); this%send_req = MPI_REQUEST_NULL
    allocate(this%recv_req(8)); this%recv_req = MPI_REQUEST_NULL

  end subroutine latlon_field2d_init

  subroutine latlon_field2d_clear(this)

    class(latlon_field2d_type), intent(inout) :: this

    call this%latlon_field_meta_clear()

    if (this%initialized .and. .not. this%linked) then
      if (associated(this%d)) then
        deallocate(this%d)
        this%d => null()
      end if
      if (associated(this%send_req)) then
        deallocate(this%send_req)
        this%send_req => null()
      end if
      if (associated(this%recv_req)) then
        deallocate(this%recv_req)
        this%recv_req => null()
      end if
    end if

    this%full_lon     = .true.
    this%full_lat     = .true.
    this%nlon         = 0
    this%nlat         = 0
    this%nlev         = 0
    this%linked       = .false.
    this%initialized  = .false.

  end subroutine latlon_field2d_clear

  subroutine latlon_field2d_copy(this, other, with_halo)

    class(latlon_field2d_type), intent(inout) :: this
    type(latlon_field2d_type), intent(in) :: other
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    integer i, j, is, ie, js, je

    if (this%loc /= other%loc) call log_error('Location does not match!', __FILE__, __LINE__)

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    if (with_halo_opt) then
      is = merge(this%mesh%full_ims, this%mesh%half_ims, this%full_lon)
      ie = merge(this%mesh%full_ime, this%mesh%half_ime, this%full_lon)
      js = merge(this%mesh%full_jms, this%mesh%half_jms, this%full_lat)
      je = merge(this%mesh%full_jme, this%mesh%half_jme, this%full_lat)
    else
      is = merge(this%mesh%full_ids, this%mesh%half_ids, this%full_lon)
      ie = merge(this%mesh%full_ide, this%mesh%half_ide, this%full_lon)
      js = merge(this%mesh%full_jds, this%mesh%half_jds, this%full_lat)
      je = merge(this%mesh%full_jde, this%mesh%half_jde, this%full_lat)
    end if

    do j = js, je
      do i = is, ie
        this%d(i,j) = other%d(i,j)
      end do
    end do

  end subroutine latlon_field2d_copy

  real(r8) function latlon_field2d_sum(this) result(res)

    class(latlon_field2d_type), intent(in) :: this

    integer is, ie, js, je

    is = merge(this%mesh%full_ids, this%mesh%half_ids, this%full_lon)
    ie = merge(this%mesh%full_ide, this%mesh%half_ide, this%full_lon)
    js = merge(this%mesh%full_jds, this%mesh%half_jds, this%full_lat)
    je = merge(this%mesh%full_jde, this%mesh%half_jde, this%full_lat)

    res = global_sum(proc%comm_model, sum(this%d(is:ie,js:je)))

  end function latlon_field2d_sum

  real(r8) function latlon_field2d_absmax(this) result(res)

    class(latlon_field2d_type), intent(in) :: this

    integer is, ie, js, je

    is = merge(this%mesh%full_ids, this%mesh%half_ids, this%full_lon)
    ie = merge(this%mesh%full_ide, this%mesh%half_ide, this%full_lon)
    js = merge(this%mesh%full_jds, this%mesh%half_jds, this%full_lat)
    je = merge(this%mesh%full_jde, this%mesh%half_jde, this%full_lat)

    res = global_max(proc%comm_model, maxval(abs(this%d(is:ie,js:je))))

  end function latlon_field2d_absmax

  subroutine latlon_field2d_link_2d(this, other)

    class(latlon_field2d_type), intent(inout) :: this
    type(latlon_field2d_type), intent(in) :: other

    if (this%initialized .and. this%loc /= other%loc) then
      call log_error('latlon_field2d_link: cannot link fields with different loc!', __FILE__, __LINE__)
    else
      call this%latlon_field_meta_link(other)
      this%full_lon = other%full_lon
      this%full_lat = other%full_lat
      this%nlon     = other%nlon
      this%nlat     = other%nlat
      if (associated(this%send_req) .and. .not. this%linked) deallocate(this%send_req)
      if (associated(this%recv_req) .and. .not. this%linked) deallocate(this%recv_req)
      this%send_req => other%send_req
      this%recv_req => other%recv_req
    end if
    if (this%initialized .and. .not. this%linked .and. associated(this%d)) deallocate(this%d)
    this%d           => other%d
    this%linked      = .true.
    this%initialized = .true.

  end subroutine latlon_field2d_link_2d

  subroutine latlon_field2d_link_3d(this, other, i3)

    class(latlon_field2d_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: other
    integer, intent(in) :: i3

    real(r8), pointer, contiguous :: tmp(:,:)
    integer is, ie, js, je

    if (this%initialized .and. .not. (this%full_lon .eqv. other%full_lon .and. this%full_lat .eqv. other%full_lat)) then
      call log_error('latlon_field2d_link: cannot link fields with different loc!', __FILE__, __LINE__)
    else
      call this%latlon_field_meta_link(other)
      this%full_lon = other%full_lon
      this%full_lat = other%full_lat
      this%nlon     = other%nlon
      this%nlat     = other%nlat
      if (associated(this%send_req) .and. .not. this%linked) deallocate(this%send_req)
      if (associated(this%recv_req) .and. .not. this%linked) deallocate(this%recv_req)
      this%send_req => other%send_req
      this%recv_req => other%recv_req
    end if
    if (this%initialized .and. .not. this%linked .and. associated(this%d)) deallocate(this%d)
    select case (this%loc)
    case ('cell', 'lev')
      is = this%mesh%full_ims; ie = this%mesh%full_ime
      js = this%mesh%full_jms; je = this%mesh%full_jme
    case ('lon')
      is = this%mesh%half_ims; ie = this%mesh%half_ime
      js = this%mesh%full_jms; je = this%mesh%full_jme
    case ('lat')
      is = this%mesh%full_ims; ie = this%mesh%full_ime
      js = this%mesh%half_jms; je = this%mesh%half_jme
    case ('vtx')
      is = this%mesh%half_ims; ie = this%mesh%half_ime
      js = this%mesh%half_jms; je = this%mesh%half_jme
    case default
      call log_error('Unhandled branch!', __FILE__, __LINE__, pid=proc%id_model)
    end select
    ! Use a temporary array pointer to fix compile error.
    tmp                 => other%d(:,:,i3)
    this%d(is:ie,js:je) => tmp
    this%linked         = .true.
    this%initialized    = .true.

  end subroutine latlon_field2d_link_3d

  subroutine latlon_field2d_final(this)

    type(latlon_field2d_type), intent(inout) :: this

    call this%clear()

  end subroutine latlon_field2d_final

  subroutine latlon_field3d_init(this, name, long_name, units, loc, mesh, halo, &
                                 halo_cross_pole, halo_diagonal, link, output, restart)

    class(latlon_field3d_type), intent(inout) :: this
    character(*), intent(in) :: name
    character(*), intent(in) :: long_name
    character(*), intent(in) :: units
    character(*), intent(in) :: loc
    type(latlon_mesh_type), intent(in), target :: mesh
    type(latlon_halo_type), intent(in), optional, target :: halo(:)
    logical, intent(in), optional :: halo_cross_pole
    logical, intent(in), optional :: halo_diagonal
    type(latlon_field3d_type), intent(in), optional, target :: link
    character(*), intent(in), optional :: output
    logical, intent(in), optional :: restart

    call this%clear()

    call this%latlon_field_meta_init(name, long_name, units, loc, mesh, halo, &
                                     halo_cross_pole, halo_diagonal, output, restart)

    if (present(link)) then
      this%d      => link%d
      this%linked = .true.
    else
      select case (loc)
      case ('cell')
        this%full_lon = .true. ; this%full_lat = .true. ; this%full_lev = .true.
        allocate(this%d(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
      case ('lon')
        this%full_lon = .false.; this%full_lat = .true. ; this%full_lev = .true.
        allocate(this%d(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
      case ('lat')
        this%full_lon = .true. ; this%full_lat = .false.; this%full_lev = .true.
        allocate(this%d(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
      case ('vtx')
        this%full_lon = .false.; this%full_lat = .false.; this%full_lev = .true.
        allocate(this%d(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
      case ('lev')
        this%full_lon = .true. ; this%full_lat = .true. ; this%full_lev = .false.
        allocate(this%d(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
      case ('lev_lon')
        this%full_lon = .false.; this%full_lat = .true. ; this%full_lev = .false.
        allocate(this%d(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
      case ('lev_lat')
        this%full_lon = .true. ; this%full_lat = .false.; this%full_lev = .false.
        allocate(this%d(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme))
      end select
    end if

    this%d           = 0
    this%nlon        = merge(mesh%full_nlon, mesh%half_nlon, this%full_lon)
    this%nlat        = merge(mesh%full_nlat, mesh%half_nlat, this%full_lat)
    this%nlev        = merge(mesh%full_nlev, mesh%half_nlev, this%full_lev)
    this%initialized = .true.

    allocate(this%send_req(8)); this%send_req = MPI_REQUEST_NULL
    allocate(this%recv_req(8)); this%recv_req = MPI_REQUEST_NULL

  end subroutine latlon_field3d_init

  subroutine latlon_field3d_clear(this)

    class(latlon_field3d_type), intent(inout) :: this

    call this%latlon_field_meta_clear()

    if (this%initialized .and. .not. this%linked) then
      if (associated(this%d)) then
        deallocate(this%d)
        this%d => null()
      end if
      if (associated(this%send_req)) then
        deallocate(this%send_req)
        this%send_req => null()
      end if
      if (associated(this%recv_req)) then
        deallocate(this%recv_req)
        this%recv_req => null()
      end if
    end if

    this%full_lon    = .true.
    this%full_lat    = .true.
    this%full_lev    = .true.
    this%nlon        = 0
    this%nlat        = 0
    this%nlev        = 0
    this%linked      = .false.
    this%initialized = .false.

  end subroutine latlon_field3d_clear

  subroutine latlon_field3d_copy_2d(this, other, k, with_halo)

    class(latlon_field3d_type), intent(inout) :: this
    type(latlon_field2d_type), intent(in) :: other
    integer, intent(in) :: k
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    integer i, j, is, ie, js, je

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    select case (this%loc)
    case ('cell')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
      end if
    case ('lon')
      if (with_halo_opt) then
        is = this%mesh%half_ims; ie = this%mesh%half_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
      else
        is = this%mesh%half_ids; ie = this%mesh%half_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
      end if
    case ('lat')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%half_jms; je = this%mesh%half_jme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%half_jds; je = this%mesh%half_jde
      end if
    case ('lev')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
      end if
    case ('vtx')
      if (with_halo_opt) then
        is = this%mesh%half_ims; ie = this%mesh%half_ime
        js = this%mesh%half_jms; je = this%mesh%half_jme
      else
        is = this%mesh%half_ids; ie = this%mesh%half_ide
        js = this%mesh%half_jds; je = this%mesh%half_jde
      end if
    case default
      call log_error('Unhandled branch in latlon_field3d_copy!', __FILE__, __LINE__)
    end select

    do j = js, je
      do i = is, ie
        this%d(i,j,k) = other%d(i,j)
      end do
    end do

  end subroutine latlon_field3d_copy_2d

  subroutine latlon_field3d_copy_3d(this, other, with_halo)

    class(latlon_field3d_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: other
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    integer i, j, k, is, ie, js, je, ks, ke

    if (this%loc /= other%loc) call log_error('Location does not match!', __FILE__, __LINE__)

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    select case (this%loc)
    case ('cell')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case ('lon')
      if (with_halo_opt) then
        is = this%mesh%half_ims; ie = this%mesh%half_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%half_ids; ie = this%mesh%half_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case ('lat')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%half_jms; je = this%mesh%half_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%half_jds; je = this%mesh%half_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case ('lev')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
        ks = this%mesh%half_kms; ke = this%mesh%half_kme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
        ks = this%mesh%half_kds; ke = this%mesh%half_kde
      end if
    case ('vtx')
      if (with_halo_opt) then
        is = this%mesh%half_ims; ie = this%mesh%half_ime
        js = this%mesh%half_jms; je = this%mesh%half_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%half_ids; ie = this%mesh%half_ide
        js = this%mesh%half_jds; je = this%mesh%half_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case default
      call log_error('Unhandled branch in latlon_field3d_copy!', __FILE__, __LINE__)
    end select

    do k = ks, ke
      do j = js, je
        do i = is, ie
          this%d(i,j,k) = other%d(i,j,k)
        end do
      end do
    end do

  end subroutine latlon_field3d_copy_3d

  subroutine latlon_field3d_copy_4d(this, other, m, with_halo)

    class(latlon_field3d_type), intent(inout) :: this
    type(latlon_field4d_type), intent(in) :: other
    integer, intent(in) :: m
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    integer i, j, k, is, ie, js, je, ks, ke

    if (this%loc /= other%loc) call log_error('Location does not match!', __FILE__, __LINE__)

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    select case (this%loc)
    case ('cell')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case default
      call log_error('Unhandled branch in latlon_field3d_copy!', __FILE__, __LINE__)
    end select

    do k = ks, ke
      do j = js, je
        do i = is, ie
          this%d(i,j,k) = other%d(i,j,k,m)
        end do
      end do
    end do

  end subroutine latlon_field3d_copy_4d

  real(r8) function latlon_field3d_sum(this) result(res)

    class(latlon_field3d_type), intent(in) :: this

    integer is, ie, js, je, ks, ke

    is = merge(this%mesh%full_ids, this%mesh%half_ids, this%full_lon)
    ie = merge(this%mesh%full_ide, this%mesh%half_ide, this%full_lon)
    js = merge(this%mesh%full_jds, this%mesh%half_jds, this%full_lat)
    je = merge(this%mesh%full_jde, this%mesh%half_jde, this%full_lat)
    ks = merge(this%mesh%full_kds, this%mesh%half_kds, this%full_lev)
    ke = merge(this%mesh%full_kde, this%mesh%half_kde, this%full_lev)

    res = global_sum(proc%comm_model, sum(this%d(is:ie,js:je,ks:ke)))

  end function latlon_field3d_sum

  subroutine latlon_field3d_add(this, other, with_halo)

    class(latlon_field3d_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: other
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    integer i, j, k, is, ie, js, je, ks, ke

    if (this%loc /= other%loc) call log_error('Location does not match!', __FILE__, __LINE__)

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    select case (this%loc)
    case ('cell')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case ('lon')
      if (with_halo_opt) then
        is = this%mesh%half_ims; ie = this%mesh%half_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%half_ids; ie = this%mesh%half_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case ('lat')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%half_jms; je = this%mesh%half_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%half_jds; je = this%mesh%half_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case ('lev')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
        ks = this%mesh%half_kms; ke = this%mesh%half_kme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
        ks = this%mesh%half_kds; ke = this%mesh%half_kde
      end if
    case ('vtx')
      if (with_halo_opt) then
        is = this%mesh%half_ims; ie = this%mesh%half_ime
        js = this%mesh%half_jms; je = this%mesh%half_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%half_ids; ie = this%mesh%half_ide
        js = this%mesh%half_jds; je = this%mesh%half_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case default
      call log_error('Unhandled branch in latlon_field3d_add!', __FILE__, __LINE__)
    end select

    do k = ks, ke
      do j = js, je
        do i = is, ie
          this%d(i,j,k) = this%d(i,j,k) + other%d(i,j,k)
        end do
      end do
    end do

  end subroutine latlon_field3d_add

  subroutine latlon_field3d_mul(this, other, with_halo)

    class(latlon_field3d_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: other
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    integer i, j, k, is, ie, js, je, ks, ke

    if (this%loc /= other%loc) call log_error('Location does not match!', __FILE__, __LINE__)

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    select case (this%loc)
    case ('cell')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case ('lon')
      if (with_halo_opt) then
        is = this%mesh%half_ims; ie = this%mesh%half_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%half_ids; ie = this%mesh%half_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case ('lat')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%half_jms; je = this%mesh%half_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%half_jds; je = this%mesh%half_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case ('lev')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
        ks = this%mesh%half_kms; ke = this%mesh%half_kme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
        ks = this%mesh%half_kds; ke = this%mesh%half_kde
      end if
    case ('vtx')
      if (with_halo_opt) then
        is = this%mesh%half_ims; ie = this%mesh%half_ime
        js = this%mesh%half_jms; je = this%mesh%half_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%half_ids; ie = this%mesh%half_ide
        js = this%mesh%half_jds; je = this%mesh%half_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case default
      call log_error('Unhandled branch in latlon_field3d_mul!', __FILE__, __LINE__)
    end select

    do k = ks, ke
      do j = js, je
        do i = is, ie
          this%d(i,j,k) = this%d(i,j,k) * other%d(i,j,k)
        end do
      end do
    end do

  end subroutine latlon_field3d_mul

  subroutine latlon_field3d_div(this, other, with_halo)

    class(latlon_field3d_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: other
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    integer i, j, k, is, ie, js, je, ks, ke

    if (this%loc /= other%loc) call log_error('Location does not match!', __FILE__, __LINE__)

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    select case (this%loc)
    case ('cell')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case ('lon')
      if (with_halo_opt) then
        is = this%mesh%half_ims; ie = this%mesh%half_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%half_ids; ie = this%mesh%half_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case ('lat')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%half_jms; je = this%mesh%half_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%half_jds; je = this%mesh%half_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case ('lev')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
        ks = this%mesh%half_kms; ke = this%mesh%half_kme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
        ks = this%mesh%half_kds; ke = this%mesh%half_kde
      end if
    case ('vtx')
      if (with_halo_opt) then
        is = this%mesh%half_ims; ie = this%mesh%half_ime
        js = this%mesh%half_jms; je = this%mesh%half_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%half_ids; ie = this%mesh%half_ide
        js = this%mesh%half_jds; je = this%mesh%half_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case default
      call log_error('Unhandled branch in latlon_field3d_div!', __FILE__, __LINE__)
    end select

    do k = ks, ke
      do j = js, je
        do i = is, ie
          this%d(i,j,k) = this%d(i,j,k) / other%d(i,j,k)
        end do
      end do
    end do

  end subroutine latlon_field3d_div

  real(r8) function latlon_field3d_min(this) result(res)

    class(latlon_field3d_type), intent(in) :: this

    integer is, ie, js, je, ks, ke

    is = merge(this%mesh%full_ids, this%mesh%half_ids, this%full_lon)
    ie = merge(this%mesh%full_ide, this%mesh%half_ide, this%full_lon)
    js = merge(this%mesh%full_jds, this%mesh%half_jds, this%full_lat)
    je = merge(this%mesh%full_jde, this%mesh%half_jde, this%full_lat)
    ks = merge(this%mesh%full_kds, this%mesh%half_kds, this%full_lev)
    ke = merge(this%mesh%full_kde, this%mesh%half_kde, this%full_lev)

    res = global_min(proc%comm_model, minval(this%d(is:ie,js:je,ks:ke)))

  end function latlon_field3d_min

  real(r8) function latlon_field3d_max(this) result(res)

    class(latlon_field3d_type), intent(in) :: this

    integer is, ie, js, je, ks, ke

    is = merge(this%mesh%full_ids, this%mesh%half_ids, this%full_lon)
    ie = merge(this%mesh%full_ide, this%mesh%half_ide, this%full_lon)
    js = merge(this%mesh%full_jds, this%mesh%half_jds, this%full_lat)
    je = merge(this%mesh%full_jde, this%mesh%half_jde, this%full_lat)
    ks = merge(this%mesh%full_kds, this%mesh%half_kds, this%full_lev)
    ke = merge(this%mesh%full_kde, this%mesh%half_kde, this%full_lev)

    res = global_max(proc%comm_model, maxval(this%d(is:ie,js:je,ks:ke)))

  end function latlon_field3d_max

  subroutine latlon_field3d_link_3d(this, other)

    class(latlon_field3d_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: other

    if (this%initialized .and. this%loc /= other%loc) then
      call log_error('latlon_field3d_link_3d: cannot link fields with different loc!', __FILE__, __LINE__)
    else
      call this%latlon_field_meta_link(other)
      this%full_lon = other%full_lon
      this%full_lat = other%full_lat
      this%full_lev = other%full_lev
      this%nlon     = other%nlon
      this%nlat     = other%nlat
      this%nlev     = other%nlev
      if (associated(this%send_req) .and. .not. this%linked) deallocate(this%send_req)
      if (associated(this%recv_req) .and. .not. this%linked) deallocate(this%recv_req)
      this%send_req => other%send_req
      this%recv_req => other%recv_req
    end if
    if (this%initialized .and. .not. this%linked .and. associated(this%d)) deallocate(this%d)
    this%d           => other%d
    this%linked      = .true.
    this%initialized = .true.

  end subroutine latlon_field3d_link_3d

  subroutine latlon_field3d_link_4d(this, other, i4)

    class(latlon_field3d_type), intent(inout) :: this
    type(latlon_field4d_type), intent(in) :: other
    integer, intent(in) :: i4

    real(r8), pointer, contiguous :: tmp(:,:,:)
    integer is, ie, js, je, ks, ke

    if (this%initialized .and. this%loc /= other%loc) then
      call log_error('latlon_field3d_link_4d: cannot link fields with different loc!', __FILE__, __LINE__)
    else
      call this%latlon_field_meta_link(other)
      this%full_lon = other%full_lon
      this%full_lat = other%full_lat
      this%full_lev = other%full_lev
      this%nlon     = other%nlon
      this%nlat     = other%nlat
      this%nlev     = other%nlev
      if (associated(this%send_req) .and. .not. this%linked) deallocate(this%send_req)
      if (associated(this%recv_req) .and. .not. this%linked) deallocate(this%recv_req)
      this%send_req => other%send_req(:,i4)
      this%recv_req => other%recv_req(:,i4)
    end if
    if (this%initialized .and. .not. this%linked .and. associated(this%d)) deallocate(this%d)
    select case (this%loc)
    case ('cell')
      is = this%mesh%full_ims; ie = this%mesh%full_ime
      js = this%mesh%full_jms; je = this%mesh%full_jme
      ks = this%mesh%full_kms; ke = this%mesh%full_kme
    case ('lon')
      is = this%mesh%half_ims; ie = this%mesh%half_ime
      js = this%mesh%full_jms; je = this%mesh%full_jme
      ks = this%mesh%full_kms; ke = this%mesh%full_kme
    case ('lat')
      is = this%mesh%full_ims; ie = this%mesh%full_ime
      js = this%mesh%half_jms; je = this%mesh%half_jme
      ks = this%mesh%full_kms; ke = this%mesh%full_kme
    case ('lev')
      is = this%mesh%full_ims; ie = this%mesh%full_ime
      js = this%mesh%full_jms; je = this%mesh%full_jme
      ks = this%mesh%half_kms; ke = this%mesh%half_kme
    case ('vtx')
      is = this%mesh%half_ims; ie = this%mesh%half_ime
      js = this%mesh%half_jms; je = this%mesh%half_jme
      ks = this%mesh%full_kms; ke = this%mesh%full_kme
    case default
      stop 'Unhandled branch in latlon_field3d_link_4d!'
    end select
    ! Use a temporary array pointer to fix compile error.
    tmp                       => other%d(:,:,:,i4)
    this%d(is:ie,js:je,ks:ke) => tmp
    this%linked               = .true.
    this%initialized          = .true.

  end subroutine latlon_field3d_link_4d

  subroutine latlon_field3d_final(this)

    type(latlon_field3d_type), intent(inout) :: this

    call this%clear()

  end subroutine latlon_field3d_final

  subroutine latlon_field4d_init(this, name, long_name, units, loc, mesh, halo, &
                                 halo_cross_pole, halo_diagonal, dim4_name, dim4_size, var4_names, output, restart)

    class(latlon_field4d_type), intent(inout) :: this
    character(*), intent(in) :: name
    character(*), intent(in) :: long_name
    character(*), intent(in) :: units
    character(*), intent(in) :: loc
    type(latlon_mesh_type), intent(in), target :: mesh
    type(latlon_halo_type), intent(in), optional, target :: halo(:)
    logical, intent(in), optional :: halo_cross_pole
    logical, intent(in), optional :: halo_diagonal
    character(*), intent(in) :: dim4_name
    integer, intent(in) :: dim4_size
    character(*), intent(in) :: var4_names(:)
    character(*), intent(in), optional :: output
    logical, intent(in), optional :: restart

    call this%clear()

    call this%latlon_field_meta_init(name, long_name, units, loc, mesh, halo, &
                                     halo_cross_pole, halo_diagonal, output, restart)

    this%dim4_name = dim4_name
    this%dim4_size = dim4_size
    allocate(this%var4_names(dim4_size)); this%var4_names = var4_names

    select case (loc)
    case ('cell')
      this%full_lon = .true. ; this%full_lat = .true. ; this%full_lev = .true.
      allocate(this%d(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme,dim4_size))
    end select

    this%d           = 0
    this%nlon        = merge(mesh%full_nlon, mesh%half_nlon, this%full_lon)
    this%nlat        = merge(mesh%full_nlat, mesh%half_nlat, this%full_lat)
    this%nlev        = merge(mesh%full_nlev, mesh%half_nlev, this%full_lev)
    this%linked      = .false.
    this%initialized = .true.

    allocate(this%send_req(8,dim4_size)); this%send_req = MPI_REQUEST_NULL
    allocate(this%recv_req(8,dim4_size)); this%recv_req = MPI_REQUEST_NULL

  end subroutine latlon_field4d_init

  subroutine latlon_field4d_clear(this)

    class(latlon_field4d_type), intent(inout) :: this

    call this%latlon_field_meta_clear()

    if (this%initialized .and. .not. this%linked) then
      if (associated(this%d)) then
        deallocate(this%d)
        this%d => null()
      end if
      if (associated(this%send_req)) then
        deallocate(this%send_req)
        this%send_req => null()
      end if
      if (associated(this%recv_req)) then
        deallocate(this%recv_req)
        this%recv_req => null()
      end if
    end if
    if (allocated(this%var4_names)) deallocate(this%var4_names)

    this%full_lon    = .true.
    this%full_lat    = .true.
    this%full_lev    = .true.
    this%nlon        = 0
    this%nlat        = 0
    this%nlev        = 0
    this%linked      = .false.
    this%initialized = .false.

  end subroutine latlon_field4d_clear

  subroutine latlon_field4d_copy_3d(this, other, m, with_halo)

    class(latlon_field4d_type), intent(in) :: this
    type(latlon_field3d_type), intent(inout) :: other
    integer, intent(in) :: m
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    integer i, j, k, is, ie, js, je, ks, ke

    if (this%loc /= other%loc) call log_error('Location does not match!', __FILE__, __LINE__)

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    select case (this%loc)
    case ('cell')
      if (with_halo_opt) then
        is = this%mesh%full_ims; ie = this%mesh%full_ime
        js = this%mesh%full_jms; je = this%mesh%full_jme
        ks = this%mesh%full_kms; ke = this%mesh%full_kme
      else
        is = this%mesh%full_ids; ie = this%mesh%full_ide
        js = this%mesh%full_jds; je = this%mesh%full_jde
        ks = this%mesh%full_kds; ke = this%mesh%full_kde
      end if
    case default
      call log_error('Unhandled branch in latlon_field4d_copy!', __FILE__, __LINE__)
    end select

    do k = ks, ke
      do j = js, je
        do i = is, ie
          this%d(i,j,k,m) = other%d(i,j,k)
        end do
      end do
    end do

  end subroutine latlon_field4d_copy_3d

  subroutine latlon_field4d_link(this, other)

    class(latlon_field4d_type), intent(inout) :: this
    type(latlon_field4d_type), intent(in) :: other

    if (this%initialized .and. this%loc /= other%loc) then
      call log_error('latlon_field4d_link: cannot link fields with different loc!', __FILE__, __LINE__)
    else
      call this%latlon_field_meta_link(other)
      this%full_lon = other%full_lon
      this%full_lat = other%full_lat
      this%full_lev = other%full_lev
      this%nlon     = other%nlon
      this%nlat     = other%nlat
      this%nlev     = other%nlev
    end if
    if (this%initialized .and. .not. this%linked .and. associated(this%d)) deallocate(this%d)
    this%d           => other%d
    this%linked      = .true.
    this%initialized = .true.

  end subroutine latlon_field4d_link

  subroutine latlon_field4d_mul(this, idx4, other, with_halo)

    class(latlon_field4d_type), intent(inout) :: this
    integer, intent(in) :: idx4
    type(latlon_field3d_type), intent(in) :: other
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    integer i, j, k, is, ie, js, je, ks, ke

    if (this%loc /= other%loc) call log_error('Location does not match!', __FILE__, __LINE__)

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    if (with_halo_opt) then
      is = this%mesh%full_ims; ie = this%mesh%full_ime
      js = this%mesh%full_jms; je = this%mesh%full_jme
      ks = this%mesh%full_kms; ke = this%mesh%full_kme
    else
      is = this%mesh%full_ids; ie = this%mesh%full_ide
      js = this%mesh%full_jds; je = this%mesh%full_jde
      ks = this%mesh%full_kds; ke = this%mesh%full_kde
    end if

    do k = ks, ke
      do j = js, je
        do i = is, ie
          this%d(i,j,k,idx4) = this%d(i,j,k,idx4) * other%d(i,j,k)
        end do
      end do
    end do

  end subroutine latlon_field4d_mul

  subroutine latlon_field4d_div(this, idx4, other, with_halo)

    class(latlon_field4d_type), intent(inout) :: this
    integer, intent(in) :: idx4
    type(latlon_field3d_type), intent(in) :: other
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    integer i, j, k, is, ie, js, je, ks, ke

    if (this%loc /= other%loc) call log_error('Location does not match!', __FILE__, __LINE__)

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    if (with_halo_opt) then
      is = this%mesh%full_ims; ie = this%mesh%full_ime
      js = this%mesh%full_jms; je = this%mesh%full_jme
      ks = this%mesh%full_kms; ke = this%mesh%full_kme
    else
      is = this%mesh%full_ids; ie = this%mesh%full_ide
      js = this%mesh%full_jds; je = this%mesh%full_jde
      ks = this%mesh%full_kds; ke = this%mesh%full_kde
    end if

    do k = ks, ke
      do j = js, je
        do i = is, ie
          this%d(i,j,k,idx4) = this%d(i,j,k,idx4) / other%d(i,j,k)
        end do
      end do
    end do

  end subroutine latlon_field4d_div

  subroutine latlon_field4d_final(this)

    type(latlon_field4d_type), intent(inout) :: this

    call this%clear()

  end subroutine latlon_field4d_final

  subroutine append_field2d(fields, name, long_name, units, loc, mesh, halo, field, &
                            halo_cross_pole, halo_diagonal, output, restart, link)

    type(array_type), intent(inout) :: fields
    character(*), intent(in) :: name
    character(*), intent(in) :: long_name
    character(*), intent(in) :: units
    character(*), intent(in) :: loc
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)
    type(latlon_field2d_type), intent(inout), target :: field
    logical, intent(in), optional :: halo_cross_pole
    logical, intent(in), optional :: halo_diagonal
    character(*), intent(in), optional :: output
    logical, intent(in), optional :: restart
    type(latlon_field2d_type), intent(in), optional :: link

    type(latlon_field2d_type), pointer :: ptr

    call field%init(                   &
      name           =name           , &
      long_name      =long_name      , &
      units          =units          , &
      loc            =loc            , &
      mesh           =mesh           , &
      halo           =halo           , &
      link           =link           , &
      output         =output         , &
      restart        =restart        , &
      halo_cross_pole=halo_cross_pole, &
      halo_diagonal  =halo_diagonal  )
    ptr => field
    call fields%append_ptr(ptr)

  end subroutine append_field2d

  subroutine append_field3d(fields, name, long_name, units, loc, mesh, halo, field, &
                            halo_cross_pole, halo_diagonal, output, restart, link)

    type(array_type), intent(inout) :: fields
    character(*), intent(in) :: name
    character(*), intent(in) :: long_name
    character(*), intent(in) :: units
    character(*), intent(in) :: loc
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)
    type(latlon_field3d_type), intent(inout), target :: field
    logical, intent(in), optional :: halo_cross_pole
    logical, intent(in), optional :: halo_diagonal
    character(*), intent(in), optional :: output
    logical, intent(in), optional :: restart
    type(latlon_field3d_type), intent(in), optional :: link

    type(latlon_field3d_type), pointer :: ptr

    call field%init(                   &
      name           =name           , &
      long_name      =long_name      , &
      units          =units          , &
      loc            =loc            , &
      mesh           =mesh           , &
      halo           =halo           , &
      halo_cross_pole=halo_cross_pole, &
      halo_diagonal  =halo_diagonal  , &
      link           =link           , &
      output         =output         , &
      restart        =restart        )
    ptr => field
    call fields%append_ptr(ptr)

  end subroutine append_field3d

  subroutine append_field4d(fields, name, long_name, units, loc, mesh, halo, field, &
                            dim4_name, dim4_size, var4_names, &
                            halo_cross_pole, halo_diagonal, output, restart)

    type(array_type), intent(inout) :: fields
    character(*), intent(in) :: name
    character(*), intent(in) :: long_name
    character(*), intent(in) :: units
    character(*), intent(in) :: loc
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)
    type(latlon_field4d_type), intent(inout), target :: field
    character(*), intent(in) :: dim4_name
    integer, intent(in) :: dim4_size
    character(*), intent(in) :: var4_names(:)
    logical, intent(in), optional :: halo_cross_pole
    logical, intent(in), optional :: halo_diagonal
    character(*), intent(in), optional :: output
    logical, intent(in), optional :: restart

    type(latlon_field4d_type), pointer :: ptr

    call field%init(                   &
      name           =name           , &
      long_name      =long_name      , &
      units          =units          , &
      loc            =loc            , &
      mesh           =mesh           , &
      halo           =halo           , &
      halo_cross_pole=halo_cross_pole, &
      halo_diagonal  =halo_diagonal  , &
      dim4_name      =dim4_name      , &
      dim4_size      =dim4_size      , &
      var4_names     =var4_names     , &
      output         =output         , &
      restart        =restart        )
    ptr => field
    call fields%append_ptr(ptr)

  end subroutine append_field4d

end module latlon_field_types_mod
