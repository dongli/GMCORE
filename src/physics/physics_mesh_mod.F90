module physics_mesh_mod

  use const_mod

  implicit none

  private

  public physics_mesh_type

  type physics_mesh_type
    integer :: ncol = 0
    integer :: nlev = 0
    integer :: cell_start_2d(3) = 0
    integer :: cell_count_2d(3) = 0
    integer :: cell_start_3d(4) = 0
    integer :: cell_count_3d(4) = 0
    integer :: lev_start(4) = 0
    integer :: lev_count(4) = 0
    real(r8) :: ptop = 0
    real(r8) :: ztop = 0
    real(r8) :: min_lon = 0
    real(r8) :: max_lon = 0
    real(r8) :: min_lat = 0
    real(r8) :: max_lat = 0
    real(r8), allocatable, dimension(:) :: lon  ! Longitude (rad)
    real(r8), allocatable, dimension(:) :: lat  ! Latitude (rad)
    real(r8), allocatable, dimension(:) :: lev  ! Vertical coordinate at full levels (1)
    real(r8), allocatable, dimension(:) :: ilev ! Vertical coordinate at half levels (1)
    real(r8), allocatable, dimension(:) :: dlev ! Vertical coordinate interval (1)
    real(r8), allocatable, dimension(:) :: area ! Cell area (m2)
    real(r8), allocatable, dimension(:) :: sin_lat
    real(r8), allocatable, dimension(:) :: cos_lat
  contains
    procedure :: init => physics_mesh_init
    procedure :: clear => physics_mesh_clear
    final physics_mesh_final
  end type physics_mesh_type

contains

  subroutine physics_mesh_init(this, ncol, nlev, lon, lat, lev, ilev, dlev, area, ptop, ztop)

    class(physics_mesh_type), intent(inout) :: this
    integer , intent(in) :: ncol
    integer , intent(in) :: nlev
    real(r8), intent(in) :: lon (ncol)
    real(r8), intent(in) :: lat (ncol)
    real(r8), intent(in) :: lev (nlev)
    real(r8), intent(in) :: ilev(nlev+1)
    real(r8), intent(in) :: dlev(nlev)
    real(r8), intent(in) :: area(ncol)
    real(r8), intent(in), optional :: ptop
    real(r8), intent(in), optional :: ztop

    call this%clear()

    this%ncol = ncol
    this%nlev = nlev
    allocate(this%lon    (ncol  )); this%lon     = lon
    allocate(this%lat    (ncol  )); this%lat     = lat
    allocate(this%lev    (nlev  )); this%lev     = lev
    allocate(this%ilev   (nlev+1)); this%ilev    = ilev
    allocate(this%dlev   (nlev  )); this%dlev    = dlev
    allocate(this%area   (ncol  )); this%area    = area
    allocate(this%sin_lat(ncol  )); this%sin_lat = sin(lat)
    allocate(this%cos_lat(ncol  )); this%cos_lat = cos(lat)
    if (present(ptop)) this%ptop = ptop
    if (present(ztop)) this%ztop = ztop

    this%min_lon = min_lon
    this%max_lon = max_lon
    this%min_lat = min_lat
    this%max_lat = max_lat

  end subroutine physics_mesh_init

  subroutine physics_mesh_clear(this)

    class(physics_mesh_type), intent(inout) :: this

    this%ncol = 0
    this%nlev = 0
    if (allocated(this%lon    )) deallocate(this%lon    )
    if (allocated(this%lat    )) deallocate(this%lat    )
    if (allocated(this%lev    )) deallocate(this%lev    )
    if (allocated(this%ilev   )) deallocate(this%ilev   )
    if (allocated(this%dlev   )) deallocate(this%dlev   )
    if (allocated(this%area   )) deallocate(this%area   )
    if (allocated(this%sin_lat)) deallocate(this%sin_lat)
    if (allocated(this%cos_lat)) deallocate(this%cos_lat)

  end subroutine physics_mesh_clear

  subroutine physics_mesh_final(this)

    type(physics_mesh_type), intent(inout) :: this

    call this%clear()

  end subroutine physics_mesh_final

end module physics_mesh_mod
