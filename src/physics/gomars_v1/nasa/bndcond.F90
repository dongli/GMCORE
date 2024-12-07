subroutine bndcond(tg, z, lnzz, z0, psi, cdm, cdh, ustar, tstar)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_pbl_mod

  implicit none

  real(r8), intent(in ) :: tg
  real(r8), intent(in ) :: z  (nlev)
  real(r8), intent(in ) :: lnzz
  real(r8), intent(in ) :: z0
  real(r8), intent(in ) :: psi(2*nlev+1,nvar)
  real(r8), intent(out) :: cdm
  real(r8), intent(out) :: cdh
  real(r8), intent(out) :: ustar
  real(r8), intent(out) :: tstar

  real(r8) dpt
  real(r8) wsp
  real(r8) rib
  real(r8) fm
  real(r8) fh

  ! For now, use psi at the bottom full level. Later interpolate to the surface.
  dpt = psi(2*nlev,3) - tg
  ! Calculate the bulk Richardson number.
  wsp = sqrt(psi(2*nlev,1)**2 + psi(2*nlev,2)**2)
  rib = g * dpt * z(nlev) / (psi(2*nlev,3) * wsp**2 + 1.0e-9_r8)
  ! Calculate the stability functions.
  if (rib >= 0) then
    fm = 1.0_r8 / (1 + 10 * rib / sqrt(1 + 5 * rib))
    fh = 1.0_r8 / (1 + 15 * rib / sqrt(1 + 5 * rib))
  else
    fm = sqrt(1 - 16 * rib)
    fh = sqrt(1 - 64 * rib)
  end if
  ! Calculate the drag coefficients.
  cdm = fm * (ka / lnzz)**2
  cdh = sqrt(fh) * ka / lnzz
  ! Calculate the friction velocity.
  ustar = sqrt(cdm) * wsp
  ! Calculate the temperature scale.
  tstar = cdh * dpt

end subroutine bndcond