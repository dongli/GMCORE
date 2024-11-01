module buffer

  !-----------------------------------------------------------------------
  !
  ! Purpose:
  !   LOW level handler for f90 arrays.
  !
  ! Author: J. Edwards
  !
  ! This file is used with genf90.pl to generate buffer.F90
  !
  !-----------------------------------------------------------------------

  use shr_kind_mod,     only: r8 => shr_kind_r8, r4=> shr_kind_r4, i4=> shr_kind_i4
  use cam_logfile,      only: iulog
  use cam_abortutils,   only: endrun

  implicit none
  private

  ! The maximum number of dims in a fortran array
#define MAXDIMS 7

  type buffer_field_default_type
    private
    real(r8), pointer :: data(:,:,:,:,:,:,:) => null()
  end type buffer_field_default_type

  type buffer_field_int
    private
    integer(i4), pointer :: data(:,:,:,:,:,:,:) => null()
  end type buffer_field_int

  type buffer_field_double
    private
    real(r8), pointer :: data(:,:,:,:,:,:,:) => null()
  end type buffer_field_double

  type buffer_field_real
    private
    real(r4), pointer :: data(:,:,:,:,:,:,:) => null()
  end type buffer_field_real

  integer(i4), parameter,public :: dtype_i4 = 1
  real(r8), parameter,public :: dtype_r8 = 1_r8
  real(r4), parameter,public :: dtype_r4 = 1_r4

  interface buffer_field_deallocate
    module procedure  buffer_field_deallocate_int
    module procedure  buffer_field_deallocate_double
    module procedure  buffer_field_deallocate_real
  end interface

  interface buffer_field_allocate
    module procedure  buffer_field_allocate_int
    module procedure  buffer_field_allocate_double
    module procedure  buffer_field_allocate_real
  end interface

  interface buffer_set_field
    module procedure buffer_set_field_const_int
    module procedure buffer_set_field_const_double
    module procedure buffer_set_field_const_real
    module procedure buffer_set_field_1d_int
    module procedure buffer_set_field_2d_int
    module procedure buffer_set_field_3d_int
    module procedure buffer_set_field_4d_int
    module procedure buffer_set_field_5d_int
    module procedure buffer_set_field_6d_int
    module procedure buffer_set_field_7d_int
    module procedure buffer_set_field_1d_double
    module procedure buffer_set_field_2d_double
    module procedure buffer_set_field_3d_double
    module procedure buffer_set_field_4d_double
    module procedure buffer_set_field_5d_double
    module procedure buffer_set_field_6d_double
    module procedure buffer_set_field_7d_double
    module procedure buffer_set_field_1d_real
    module procedure buffer_set_field_2d_real
    module procedure buffer_set_field_3d_real
    module procedure buffer_set_field_4d_real
    module procedure buffer_set_field_5d_real
    module procedure buffer_set_field_6d_real
    module procedure buffer_set_field_7d_real
  end interface

  interface buffer_get_field_ptr
    module procedure buffer_get_field_ptr_1d_int
    module procedure buffer_get_field_ptr_2d_int
    module procedure buffer_get_field_ptr_3d_int
    module procedure buffer_get_field_ptr_4d_int
    module procedure buffer_get_field_ptr_5d_int
    module procedure buffer_get_field_ptr_6d_int
    module procedure buffer_get_field_ptr_7d_int
    module procedure buffer_get_field_ptr_1d_double
    module procedure buffer_get_field_ptr_2d_double
    module procedure buffer_get_field_ptr_3d_double
    module procedure buffer_get_field_ptr_4d_double
    module procedure buffer_get_field_ptr_5d_double
    module procedure buffer_get_field_ptr_6d_double
    module procedure buffer_get_field_ptr_7d_double
    module procedure buffer_get_field_ptr_1d_real
    module procedure buffer_get_field_ptr_2d_real
    module procedure buffer_get_field_ptr_3d_real
    module procedure buffer_get_field_ptr_4d_real
    module procedure buffer_get_field_ptr_5d_real
    module procedure buffer_get_field_ptr_6d_real
    module procedure buffer_get_field_ptr_7d_real
  end interface

  public buffer_field_deallocate, buffer_field_allocate, buffer_set_field, buffer_get_field_ptr, buffer_field_default_type
  public buffer_field_is_alloc

contains

  subroutine buffer_field_deallocate_int(bfg, dtype)

    type(buffer_field_default_type), intent(inout) :: bfg
    integer(i4), intent(in) :: dtype

    type(buffer_field_int) b1

    b1 = transfer(bfg, b1)

    if (.not.associated(b1%data)) then
      call endrun('Attempt to deallocate unassociated array pointer!')
    end if

    deallocate(b1%data)

    nullify(bfg%data)

  end subroutine buffer_field_deallocate_int

  subroutine buffer_field_deallocate_double(bfg, dtype)

    type(buffer_field_default_type), intent(inout) :: bfg
    real(r8), intent(in) :: dtype

    type(buffer_field_double) b1

    b1 = transfer(bfg, b1)

    if (.not.associated(b1%data)) then
      call endrun('Attempt to deallocate unassociated array pointer!')
    end if

    deallocate(b1%data)

    nullify(bfg%data)

  end subroutine buffer_field_deallocate_double

  subroutine buffer_field_deallocate_real(bfg, dtype)

    type(buffer_field_default_type), intent(inout) :: bfg
    real(r4), intent(in) :: dtype

    type(buffer_field_real) b1

    b1 = transfer(bfg, b1)

    if (.not.associated(b1%data)) then
       call endrun('Attempt to deallocate unassociated array pointer!')
    end if

    deallocate(b1%data)

    nullify(bfg%data)

  end subroutine buffer_field_deallocate_real

  logical function buffer_field_is_alloc(bfg)

    type(buffer_field_default_type),intent(in) :: bfg

    buffer_field_is_alloc = associated(bfg%data)

  end function buffer_field_is_alloc

  subroutine buffer_field_allocate_int(bfg, dimsizes, dtype)

    type(buffer_field_default_type), intent(inout) :: bfg
    integer, intent(in) :: dimsizes(:)
    integer(i4), intent(in) :: dtype

    integer alldimsizes(MAXDIMS)
    integer ierr
    type(buffer_field_int) b1

    alldimsizes(:) = 1
    alldimsizes(1:size(dimsizes)) = dimsizes

    if (associated(bfg%data)) then
      call endrun('Attempt to allocate array to associated pointer!')
    end if

    allocate(b1%data(alldimsizes(1),alldimsizes(2),alldimsizes(3),alldimsizes(4), &
                     alldimsizes(5),alldimsizes(6),alldimsizes(7)), stat=ierr)

    if (ierr/=0) then
      call endrun("allocate failed")
    end if

    b1%data = 0
    bfg = transfer(b1, bfg)

  end subroutine buffer_field_allocate_int

  subroutine buffer_field_allocate_double(bfg, dimsizes, dtype)

    type(buffer_field_default_type), intent(inout) :: bfg
    integer, intent(in) :: dimsizes(:)
    real(r8), intent(in) :: dtype

    integer alldimsizes(MAXDIMS)
    integer ierr
    type(buffer_field_double) b1

    alldimsizes(:) = 1
    alldimsizes(1:size(dimsizes)) = dimsizes

    if (associated(bfg%data)) then
       call endrun('Attempt to allocate array to associated pointer!')
    end if

    allocate(b1%data(alldimsizes(1),alldimsizes(2),alldimsizes(3),alldimsizes(4), &
                     alldimsizes(5),alldimsizes(6),alldimsizes(7)), stat=ierr)

    if (ierr/=0) then
      call endrun("allocate failed")
    end if

    b1%data = 0
    bfg = transfer(b1, bfg)

  end subroutine buffer_field_allocate_double

  subroutine buffer_field_allocate_real (bfg, dimsizes, dtype)

    type(buffer_field_default_type), intent(inout) :: bfg
    integer, intent(in) :: dimsizes(:)
    real(r4), intent(in) :: dtype

    integer alldimsizes(MAXDIMS)
    integer ierr
    type(buffer_field_real) b1

    alldimsizes(:) = 1
    alldimsizes(1:size(dimsizes)) = dimsizes

    if (associated(bfg%data)) then
       call endrun('Attempt to allocate array to associated pointer!')
    end if

    allocate(b1%data(alldimsizes(1),alldimsizes(2),alldimsizes(3),alldimsizes(4), &
                     alldimsizes(5),alldimsizes(6),alldimsizes(7)), stat=ierr)

    if (ierr/=0) then
      call endrun("allocate failed")
    end if

    b1%data = 0
    bfg = transfer(b1, bfg)

  end subroutine buffer_field_allocate_real

  subroutine buffer_set_field_const_int(bfg, const, start, kount)

    type(buffer_field_default_type) :: bfg
    integer(i4), intent(in) :: const
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_int) ptr

    integer i, ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      ns = size(start)
      strt(1:ns) = start

      fin = strt + cnt - 1

      do i = 1, ns
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do
      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4), &
               strt(5):fin(5),strt(6):fin(6),strt(7):fin(7)) = const
    else
      ptr%data = const
    end if

  end subroutine buffer_set_field_const_int

  subroutine buffer_set_field_const_double(bfg, const, start,kount)

    type(buffer_field_default_type) :: bfg
    real(r8), intent(in) :: const
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_double) ptr

    integer i, ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      ns = size(start)
      strt(1:ns) = start

      fin = strt + cnt - 1

      do i = 1, ns
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4), &
               strt(5):fin(5),strt(6):fin(6),strt(7):fin(7)) = const
    else
      ptr%data = const
    end if

  end subroutine buffer_set_field_const_double

  subroutine buffer_set_field_const_real(bfg, const, start,kount)

    type(buffer_field_default_type) :: bfg
    real(r4), intent(in) :: const
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_real) ptr

    integer i, ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    if (present(start).and.present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      ns = size(start)
      strt(1:ns) = start

      fin = strt + cnt - 1

      do i = 1, ns
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4), &
               strt(5):fin(5),strt(6):fin(6),strt(7):fin(7)) = const
    else
      ptr%data = const
    end if

  end subroutine buffer_set_field_const_real

  !
  ! Given a physics_buffer chunk and an index return a pointer to a field chunk
  !

  subroutine buffer_get_field_ptr_1d_int(bfg, field, start, kount)

    type(buffer_field_default_type), intent(in) :: bfg
    integer(i4), pointer :: field(:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_int), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2),strt(3),strt(4),strt(5),strt(6),strt(7))

  end subroutine buffer_get_field_ptr_1d_int

  subroutine buffer_get_field_ptr_2d_int(bfg, field, start, kount)

    type(buffer_field_default_type), intent(in) :: bfg
    integer(i4), pointer :: field(:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_int), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3),strt(4),strt(5),strt(6),strt(7))

  end subroutine buffer_get_field_ptr_2d_int

  subroutine buffer_get_field_ptr_3d_int(bfg, field, start, kount)

    type(buffer_field_default_type), intent(in) :: bfg
    integer(i4), pointer :: field(:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_int), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4),strt(5),strt(6),strt(7))

  end subroutine buffer_get_field_ptr_3d_int

  subroutine buffer_get_field_ptr_4d_int(bfg, field, start,kount)

    type(buffer_field_default_type), intent(in) :: bfg
    integer(i4), pointer :: field(:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_int), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4),strt(5),strt(6),strt(7))

  end subroutine buffer_get_field_ptr_4d_int

  subroutine buffer_get_field_ptr_5d_int(bfg, field, start,kount)

    type(buffer_field_default_type), intent(in) :: bfg
    integer(i4), pointer :: field(:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_int), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4),strt(5):fin(5),strt(6),strt(7))

  end subroutine buffer_get_field_ptr_5d_int

  subroutine buffer_get_field_ptr_6d_int(bfg, field, start, kount)

    type(buffer_field_default_type), intent(in) :: bfg
    integer(i4), pointer :: field(:,:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_int), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4), &
                      strt(5):fin(5),strt(6):fin(6),strt(7))

  end subroutine buffer_get_field_ptr_6d_int

  subroutine buffer_get_field_ptr_7d_int(bfg, field, start, kount)

    type(buffer_field_default_type), intent(in) :: bfg
    integer(i4), pointer :: field(:,:,:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_int), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4), &
                      strt(5):fin(5),strt(6):fin(6),strt(7):fin(7))

  end subroutine buffer_get_field_ptr_7d_int

  subroutine buffer_get_field_ptr_1d_double(bfg, field, start,kount)

    type(buffer_field_default_type), intent(in) :: bfg
    real(r8), pointer :: field(:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_double), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2),strt(3),strt(4),strt(5),strt(6),strt(7))

  end subroutine buffer_get_field_ptr_1d_double

  subroutine buffer_get_field_ptr_2d_double(bfg, field, start,kount)

    type(buffer_field_default_type), intent(in) :: bfg
    real(r8), pointer :: field(:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_double), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3),strt(4),strt(5),strt(6),strt(7))

  end subroutine buffer_get_field_ptr_2d_double

  subroutine buffer_get_field_ptr_3d_double(bfg, field, start,kount)

    type(buffer_field_default_type), intent(in) :: bfg
    real(r8), pointer :: field(:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_double), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4),strt(5),strt(6),strt(7))

  end subroutine buffer_get_field_ptr_3d_double

  subroutine buffer_get_field_ptr_4d_double(bfg, field, start,kount)

    type(buffer_field_default_type), intent(in) :: bfg
    real(r8), pointer :: field(:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_double), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4),strt(5),strt(6),strt(7))

  end subroutine buffer_get_field_ptr_4d_double

  subroutine buffer_get_field_ptr_5d_double(bfg, field, start, kount)

    type(buffer_field_default_type), intent(in) :: bfg
    real(r8), pointer :: field(:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_double), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4),strt(5):fin(5),strt(6),strt(7))

  end subroutine buffer_get_field_ptr_5d_double

  subroutine buffer_get_field_ptr_6d_double(bfg, field, start, kount)

    type(buffer_field_default_type), intent(in) :: bfg
    real(r8), pointer :: field(:,:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_double), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4), &
                      strt(5):fin(5),strt(6):fin(6),strt(7))

  end subroutine buffer_get_field_ptr_6d_double

  subroutine buffer_get_field_ptr_7d_double(bfg, field, start, kount)

    type(buffer_field_default_type), intent(in) :: bfg
    real(r8), pointer :: field(:,:,:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_double), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4), &
                      strt(5):fin(5),strt(6):fin(6),strt(7):fin(7))

  end subroutine buffer_get_field_ptr_7d_double

  subroutine buffer_get_field_ptr_1d_real(bfg, field, start, kount)

    type(buffer_field_default_type), intent(in) :: bfg
    real(r4), pointer :: field(:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_real), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2),strt(3),strt(4),strt(5),strt(6),strt(7))

  end subroutine buffer_get_field_ptr_1d_real

  subroutine buffer_get_field_ptr_2d_real(bfg, field, start, kount)

    type(buffer_field_default_type), intent(in) :: bfg
    real(r4), pointer :: field(:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_real), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3),strt(4),strt(5),strt(6),strt(7))

  end subroutine buffer_get_field_ptr_2d_real

  subroutine buffer_get_field_ptr_3d_real(bfg, field, start, kount)

    type(buffer_field_default_type), intent(in) :: bfg
    real(r4), pointer :: field(:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_real), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4),strt(5),strt(6),strt(7))

  end subroutine buffer_get_field_ptr_3d_real

  subroutine buffer_get_field_ptr_4d_real(bfg, field, start, kount)

    type(buffer_field_default_type), intent(in) :: bfg
    real(r4), pointer :: field(:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_real), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4),strt(5),strt(6),strt(7))

  end subroutine buffer_get_field_ptr_4d_real

  subroutine buffer_get_field_ptr_5d_real(bfg, field, start, kount)
    type(buffer_field_default_type), intent(in) :: bfg
    real(r4), pointer :: field(:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_real), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4),strt(5):fin(5),strt(6),strt(7))

  end subroutine buffer_get_field_ptr_5d_real

  subroutine buffer_get_field_ptr_6d_real(bfg, field, start, kount)

    type(buffer_field_default_type), intent(in) :: bfg
    real(r4), pointer :: field(:,:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_real), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4), &
                      strt(5):fin(5),strt(6):fin(6),strt(7))

  end subroutine buffer_get_field_ptr_6d_real

  subroutine buffer_get_field_ptr_7d_real(bfg, field, start, kount)

    type(buffer_field_default_type), intent(in) :: bfg
    real(r4), pointer :: field(:,:,:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_real), target :: ptr
    integer ns, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)

    strt(:) = 1
    cnt = shape(ptr%data)

    if (present(start)) then
      ns = size(start)
      strt(1:ns) = start
    end if
    if (present(kount)) then
      cnt(1:ns) = kount
    end if
    fin = strt + cnt - 1

    field => ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4), &
                      strt(5):fin(5),strt(6):fin(6),strt(7):fin(7))

  end subroutine buffer_get_field_ptr_7d_real

  subroutine buffer_set_field_1d_int(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    integer(i4), intent(in) :: field(:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_int) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2),strt(3),strt(4),strt(5),strt(6),strt(7)) = field
    else
      ptr%data(:,1,1,1,1,1,1) = field
    end if

  end subroutine buffer_set_field_1d_int

  subroutine buffer_set_field_2d_int(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    integer(i4), intent(in) :: field(:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_int) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3),strt(4),strt(5),strt(6),strt(7)) = field
    else
      ptr%data(:,:,1,1,1,1,1) = field
    end if

  end subroutine buffer_set_field_2d_int

  subroutine buffer_set_field_3d_int(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    integer(i4), intent(in) :: field(:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_int) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4),strt(5),strt(6),strt(7)) = field
    else
      ptr%data(:,:,:,1,1,1,1) = field
    end if

  end subroutine buffer_set_field_3d_int

  subroutine buffer_set_field_4d_int(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    integer(i4), intent(in) :: field(:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_int) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg, ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4),strt(5),strt(6),strt(7)) = field
    else
      ptr%data(:,:,:,:,1,1,1) = field
    end if

  end subroutine buffer_set_field_4d_int

  subroutine buffer_set_field_5d_int(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    integer(i4), intent(in) :: field(:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_int) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4),strt(5):fin(5),strt(6),strt(7)) = field
    else
      ptr%data(:,:,:,:,:,1,1) = field
    end if

  end subroutine buffer_set_field_5d_int

  subroutine buffer_set_field_6d_int(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    integer(i4), intent(in) :: field(:,:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_int) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4), &
               strt(5):fin(5),strt(6):fin(6),strt(7)) = field
    else
      ptr%data(:,:,:,:,:,:,1) = field
    end if

  end subroutine buffer_set_field_6d_int

  subroutine buffer_set_field_7d_int(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    integer(i4), intent(in) :: field(:,:,:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_int) ptr

    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

       ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4), &
                strt(5):fin(5),strt(6):fin(6),strt(7):fin(7)) = field
    else
      ptr%data = field
    end if

  end subroutine buffer_set_field_7d_int

  subroutine buffer_set_field_1d_double(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    real(r8), intent(in) :: field(:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_double) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2),strt(3),strt(4),strt(5),strt(6),strt(7)) = field
    else
      ptr%data(:,1,1,1,1,1,1) = field
    end if

  end subroutine buffer_set_field_1d_double

  subroutine buffer_set_field_2d_double(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    real(r8), intent(in) :: field(:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_double) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3),strt(4),strt(5),strt(6),strt(7)) = field
    else
      ptr%data(:,:,1,1,1,1,1) = field
    end if

  end subroutine buffer_set_field_2d_double

  subroutine buffer_set_field_3d_double(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    real(r8), intent(in) :: field(:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_double) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4),strt(5),strt(6),strt(7)) = field
    else
      ptr%data(:,:,:,1,1,1,1) = field
    end if

  end subroutine buffer_set_field_3d_double

  subroutine buffer_set_field_4d_double(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    real(r8), intent(in) :: field(:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_double) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4),strt(5),strt(6),strt(7)) = field
    else
      ptr%data(:,:,:,:,1,1,1) = field
    end if

  end subroutine buffer_set_field_4d_double

  subroutine buffer_set_field_5d_double(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    real(r8), intent(in) :: field(:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_double) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4),strt(5):fin(5),strt(6),strt(7)) = field
    else
      ptr%data(:,:,:,:,:,1,1) = field
    end if

  end subroutine buffer_set_field_5d_double

  subroutine buffer_set_field_6d_double(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    real(r8), intent(in) :: field(:,:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_double) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4), &
               strt(5):fin(5),strt(6):fin(6),strt(7)) = field
    else
      ptr%data(:,:,:,:,:,:,1) = field
    end if

  end subroutine buffer_set_field_6d_double

  subroutine buffer_set_field_7d_double(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    real(r8), intent(in) :: field(:,:,:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_double) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4), &
               strt(5):fin(5),strt(6):fin(6),strt(7):fin(7)) = field
    else
      ptr%data = field
    end if

  end subroutine buffer_set_field_7d_double

  subroutine buffer_set_field_1d_real(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    real(r4), intent(in) :: field(:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_real) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2),strt(3),strt(4),strt(5),strt(6),strt(7)) = field
    else
      ptr%data(:,1,1,1,1,1,1) = field
    end if

  end subroutine buffer_set_field_1d_real

  subroutine buffer_set_field_2d_real(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    real(r4), intent(in) :: field(:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_real) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3),strt(4),strt(5),strt(6),strt(7)) = field
    else
      ptr%data(:,:,1,1,1,1,1) = field
    end if

  end subroutine buffer_set_field_2d_real

  subroutine buffer_set_field_3d_real(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    real(r4), intent(in) :: field(:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_real) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4),strt(5),strt(6),strt(7)) = field
    else
      ptr%data(:,:,:,1,1,1,1) = field
    end if

  end subroutine buffer_set_field_3d_real

  subroutine buffer_set_field_4d_real(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    real(r4), intent(in) :: field(:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_real) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4),strt(5),strt(6),strt(7)) = field
    else
      ptr%data(:,:,:,:,1,1,1) = field
    end if

  end subroutine buffer_set_field_4d_real

  subroutine buffer_set_field_5d_real(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    real(r4), intent(in) :: field(:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_real) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4),strt(5):fin(5),strt(6),strt(7)) = field
    else
      ptr%data(:,:,:,:,:,1,1) = field
    end if

  end subroutine buffer_set_field_5d_real

  subroutine buffer_set_field_6d_real(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    real(r4), intent(in) :: field(:,:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)
    type(buffer_field_real) ptr

    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4), &
               strt(5):fin(5),strt(6):fin(6),strt(7)) = field
    else
      ptr%data(:,:,:,:,:,:,1) = field
    end if

  end subroutine buffer_set_field_6d_real

  subroutine buffer_set_field_7d_real(bfg, field, start, kount)

    type(buffer_field_default_type) :: bfg
    real(r4), intent(in) :: field(:,:,:,:,:,:,:)
    integer, intent(in), optional :: start(:), kount(:)

    type(buffer_field_real) ptr
    integer i, nc, strt(7), fin(7), cnt(7)

    ptr = transfer(bfg,ptr)
    if (present(start) .and. present(kount)) then
      strt(:) = 1
      cnt = shape(ptr%data)

      nc = size(start)
      strt(1:nc) = start
      fin = strt + cnt - 1

      do i = 1, nc
        fin(i) = strt(i) + kount(i) - 1
        if (strt(i) < 1 .or. fin(i) > cnt(i)) then
          call endrun('Start plus kount exceeds dimension bounds!')
        end if
      end do

      ptr%data(strt(1):fin(1),strt(2):fin(2),strt(3):fin(3),strt(4):fin(4), &
               strt(5):fin(5),strt(6):fin(6),strt(7):fin(7)) = field
    else
      ptr%data = field
    end if

  end subroutine buffer_set_field_7d_real

end module buffer
