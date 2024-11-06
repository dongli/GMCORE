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
! Description:
!
!   This module calculates ustar and tstar from bulk Richardson number.
! ==============================================================================

module gomars_v2_sfc_mod

  use gomars_v2_const_mod
  use gomars_v2_namelist_mod
  use gomars_v2_types_mod

  implicit none

  private

  public gomars_v2_sfc_run

contains

  subroutine gomars_v2_sfc_run(state)

    type(gomars_v2_state_type), intent(inout) :: state

    real(r8) dpt, rib, lnz
    integer icol, nlev

    associate (mesh    => state%mesh   , &
               pt      => state%pt     , & ! in
               z       => state%z      , & ! in
               tg      => state%tg     , & ! in
               wsp_bot => state%wsp_bot, & ! in
               z0      => state%z0     , & ! in
               fm      => state%fm     , & ! out
               fh      => state%fh     , & ! out
               cdm     => state%cdm    , & ! out
               cdh     => state%cdh    , & ! out
               ustar   => state%ustar  , & ! out
               tstar   => state%tstar  )   ! out
    nlev = mesh%nlev
    do icol = 1, mesh%ncol
      dpt = pt(icol,nlev) - tg(icol)
      ! FIXME: Check if z(icol,nlev) is about 5 m.
      ! Calculate bulk Richardson number.
      rib = g * z(icol,nlev) * dpt / (pt(icol,nlev) * wsp_bot(icol)**2 + 1.0e-9_r8)
      ! Calculate stability functions (m for momentum, h for heat).
      if (rib >= 0) then
        fm(icol) = 1.0_r8 / (1 + (10 * rib / sqrt(1 + 5 * rib)))
        fh(icol) = 1.0_r8 / (1 + (15 * rib / sqrt(1 + 5 * rib)))
      else
        fm(icol) = sqrt(1 - 16 * rib)
        fh(icol) = sqrt(1 - 64 * rib)
      end if
      lnz = log(z(icol,nlev) / z0(icol))
      ! Calculate drag coefficients.
      cdm(icol) = fm(icol) * (ka / lnz)**2
      cdh(icol) = sqrt(fh(icol)) * ka / lnz
      ustar = sqrt(cdm(icol)) * wsp_bot(icol)
      tstar = cdh(icol) * dpt
    end do
    end associate

  end subroutine gomars_v2_sfc_run

end module gomars_v2_sfc_mod
