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

module gomars_v2_spectra_mod

  use gomars_v2_const_mod

  implicit none

  private

  public gomars_v2_spectra_init
  public gomars_v2_spectra_final
  public spec_vs
  public spec_ir

  type spectra_type
    integer :: n = 0                           ! Number of spectral intervals
    real(r8), allocatable, dimension(:) :: wn  ! Wavenumbers at spectral interval centers
    real(r8), allocatable, dimension(:) :: wl  ! Wavelengths at spectral interval centers
    real(r8), allocatable, dimension(:) :: dwn ! Delta wavenumber of each spectral interval (cm-1)
    real(r8), allocatable, dimension(:) :: bwn ! Wavenumbers at spectral interval edges
  contains
    procedure :: init  => spectra_init
    procedure :: clear => spectra_clear
    final spectra_final
  end type spectra_type

  integer, parameter :: nspec_vs = 7  ! Number of visible spectral intervals
  integer, parameter :: nspec_ir = 5  ! Number of IR spectral intervals
  type(spectra_type) spec_vs
  type(spectra_type) spec_ir

contains

  subroutine gomars_v2_spectra_init()

    call spec_vs%init(nspec_vs, [ &
       2222.22_r8, & ! ->   4.50 microns
       3087.37_r8, & ! ->   3.24 microns
       4030.63_r8, & ! ->   2.48 microns
       5370.57_r8, & ! ->   1.86 microns
       7651.11_r8, & ! ->   1.31 microns
      12500.00_r8, & ! ->   0.80 microns
      25000.00_r8, & ! ->   0.40 microns
      41666.67_r8  & ! ->   0.24 microns
    ])

    call spec_ir %init(nspec_ir, [ &
        10.000_r8, & ! -> 1000.00 microns
       166.667_r8, & ! ->   60.00 microns
       416.667_r8, & ! ->   24.00 microns
       833.333_r8, & ! ->   12.00 microns
      1250.000_r8, & ! ->    8.00 microns
      2222.222_r8  & ! ->    4.50 microns
    ])

  end subroutine gomars_v2_spectra_init

  subroutine gomars_v2_spectra_final()

    call spec_vs%clear()
    call spec_ir %clear()

  end subroutine gomars_v2_spectra_final

  subroutine spectra_init(this, n, bwn)

    class(spectra_type), intent(inout) :: this
    integer , intent(in) :: n        ! Bin number
    real(r8), intent(in) :: bwn(n+1) ! Wavenumber at bin edges (cm-1)

    integer i

    call this%clear()

    allocate(this%wn (n  ))
    allocate(this%wl (n  ))
    allocate(this%dwn(n  ))
    allocate(this%bwn(n+1))

    this%n = n
    this%bwn = bwn
    do i = 1, n
      this%wn(i) = 0.5_r8 * (bwn(i) + bwn(i+1))
      this%wl(i) = 1.0e4_r8 / this%wn(i)
      this%dwn(i) = bwn(i+1) - bwn(i)
    end do

  end subroutine spectra_init

  subroutine spectra_clear(this)

    class(spectra_type), intent(inout) :: this

    this%n = 0
    if (allocated(this%wn )) deallocate(this%wn )
    if (allocated(this%wl )) deallocate(this%wl )
    if (allocated(this%dwn)) deallocate(this%dwn)
    if (allocated(this%bwn)) deallocate(this%bwn)

  end subroutine spectra_clear

  subroutine spectra_final(this)

    type(spectra_type), intent(inout) :: this

    call this%clear()

  end subroutine spectra_final

end module gomars_v2_spectra_mod
