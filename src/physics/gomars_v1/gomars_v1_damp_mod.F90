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

module gomars_v1_damp_mod

  use gomars_v1_const_mod
  use vert_coord_mod

  implicit none

  private

  public gomars_v1_damp_init
  public gomars_v1_damp_run
  public gomars_v1_damp_final

  ! Number of top levels for Rayleigh friction
  integer , parameter :: lray       = 3
  integer , parameter :: k1         = 4
  real(r8), parameter :: alfray     = 2
  ! Relaxation time in days
  real(r8), parameter :: trefr      = 0.5_r8

  real(r8), allocatable :: rayk(:)

contains

  subroutine gomars_v1_damp_init(nlev)

    integer, intent(in) :: nlev

    integer k
    real(r8) pl(nlev)

    call gomars_v1_damp_final()

    ! Set a reference pressure profile.
    do k = 1, nlev
      pl(k) = vert_coord_calc_mg(k, psl)
    end do

    allocate(rayk(nlev))

    rayk = 0

    do k = 1, lray
      rayk(k) = (pl(k1) / pl(k))**alfray / (trefr * mars_sol_seconds)
    end do

  end subroutine gomars_v1_damp_init

  subroutine gomars_v1_damp_run(nlev, upi, vpi, om, teta)

    integer , intent(in   ) :: nlev
    real(r8), intent(inout), dimension(2*nlev+3) :: upi
    real(r8), intent(inout), dimension(2*nlev+3) :: vpi
    real(r8), intent(in   ), dimension(2*nlev+3) :: om
    real(r8), intent(inout), dimension(2*nlev+3) :: teta

    integer k, l
    real(r8) rkei, rkef

    do k = 1, lray
      l = 2 * k + 2
      rkei = 0.5_r8 * (upi(l)**2 + vpi(l)**2)
      upi(l) = upi(l) - rayk(k) * upi(l) * dt
      vpi(l) = vpi(l) - rayk(k) * vpi(l) * dt
      rkef = 0.5_r8 * (upi(l)**2 + vpi(l)**2)
      teta(l) = teta(l) + (rkei - rkef) / (om(l) * cpd)
    end do

  end subroutine gomars_v1_damp_run

  subroutine gomars_v1_damp_final()

    if (allocated(rayk)) deallocate(rayk)

  end subroutine gomars_v1_damp_final

end module gomars_v1_damp_mod