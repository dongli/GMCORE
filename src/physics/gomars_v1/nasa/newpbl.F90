subroutine newpbl(z0, ps, pl, teta, om, u, v, q, z, z_lev)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_tracers_mod
  use gomars_v1_pbl_mod

  implicit none

  real(r8), intent(in) :: z0
  real(r8), intent(in) :: ps
  real(r8), intent(in) :: pl   (2*nlev+3)
  real(r8), intent(in) :: teta (2*nlev+3)
  real(r8), intent(in) :: om   (2*nlev+3)
  real(r8), intent(in) :: u    (  nlev  )
  real(r8), intent(in) :: v    (  nlev  )
  real(r8), intent(in) :: q    (  nlev  ,ntracers)
  real(r8), intent(in) :: z    (  nlev  )
  real(r8), intent(in) :: z_lev(  nlev+1)

  integer k, l
  real(r8) lnzz
  real(r8) psi  (2*nlev+1,nvar)
  real(r8) rho  (2*nlev+1)

  do k = 1, nlev
    l = 2 * k
    psi(l,1) = u   (k)
    psi(l,2) = v   (k)
    psi(l,3) = teta(l+2)
    psi(l,4) = q   (k,iMa_vap)
  end do
  do k = 1, nlev + 1
    l = 2 * k - 1
    psi(l,3) = teta(l+2)
  end do

  lnzz = log(z(nlev) / z0)

  do k = 1, 2 * nlev + 1
    rho(k) = pl(k+2) / rd / (om(k+2) * teta(k+2))
  end do

end subroutine newpbl