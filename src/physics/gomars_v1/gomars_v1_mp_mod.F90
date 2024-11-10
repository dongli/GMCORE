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

module gomars_v1_mp_mod

  use gomars_v1_const_mod
  use gomars_v1_tracers_mod

  implicit none

  private

  public mp_init
  public mp_final

  real(r8), allocatable :: aerdens(:), stdv(:)

contains

  subroutine mp_init()

    call mp_final()

    allocate(aerdens(ntracers))
    allocate(stdv   (ntracers))

    aerdens(iMa_dt ) = dpden_dt
    aerdens(iNb_dt ) = dpden_dt
    aerdens(iMa_cld) = dpden_ice
    aerdens(iNb_cld) = dpden_ice
    aerdens(iMa_cor) = dpden_ice

    stdv(iMa_dt ) = dev_dt
    stdv(iNb_dt ) = dev_dt
    stdv(iMa_cld) = dev_ice
    stdv(iNb_cld) = dev_ice
    stdv(iMa_cor) = dev_ice

  end subroutine mp_init

  subroutine mp_final()

    if (allocated(aerdens)) deallocate(aerdens)
    if (allocated(stdv   )) deallocate(stdv   )

  end subroutine mp_final

end module gomars_v1_mp_mod