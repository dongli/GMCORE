subroutine sfluxv(dtauv, tauv, taucumv, taugsurf, sol, cosz, als, wbarv, cosbv, &
                  fluxupv, fluxdnv, fmnetv, nfluxtopv, diffvt, detau)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod, only: r8, nlev, nspectv, ngauss, nlayrad, nlevrad
  use gomars_v1_rad_mod

  implicit none

  real(r8), intent(in ) :: dtauv   (nlayrad ,nspectv,ngauss)
  real(r8), intent(in ) :: tauv    (nlevrad ,nspectv,ngauss)
  real(r8), intent(in ) :: taucumv (2*nlev+3,nspectv,ngauss)
  real(r8), intent(in ) :: taugsurf(         nspectv,ngauss-1)
  real(r8), intent(in ) :: sol     (         nspectv       )
  real(r8), intent(in ) :: cosz
  real(r8), intent(in ) :: als
  real(r8), intent(in ) :: wbarv   (nlayrad ,nspectv,ngauss)
  real(r8), intent(in ) :: cosbv   (nlayrad ,nspectv,ngauss)
  real(r8), intent(out) :: fluxupv (nlayrad)
  real(r8), intent(out) :: fluxdnv (nlayrad)
  real(r8), intent(out) :: fmnetv  (nlayrad)
  real(r8), intent(out) :: nfluxtopv
  real(r8), intent(out) :: diffvt
  real(r8), intent(out) :: detau   (         nspectv,ngauss)

  integer n, k, l, s, g
  real(r8) fzero
  real(r8) btop
  real(r8) bsfc

  n = 2 * nlev + 3

  fluxupv   = 0
  fluxdnv   = 0
  fmnetv    = 0
  nfluxtopv = 0
  diffvt    = 0

  ! Calculate the net flux in each spectral interval.
  do s = 1, nspectv
    fzero = fzerov(s)
    if (fzerov(s) < 0.99_r8) then
      do g = 1, ngauss - 1
        call getdetau(    &
          dtauv  (1,s,g), &
          tauv   (1,s,g), &
          taucumv(1,s,g), &
          wbarv  (1,s,g), &
          cosbv  (1,s,g), &
          detau  (  s,g)  &
        )
        if (taugsurf(s,g) < tlimits) then
          fzero = fzero + (1 - fzerov(s)) * gweight(g)
        else
          ! Set up the top and bottom boundary conditions.
          btop = 0
          bsfc = als * cosz * sol(s) * exp(-min(detau(s,g) / cosz, maxexp))
          ! We can now solve for the coefficients of the two-stream problem.
          ! call gfluxv
        end if
      end do
    else

    end if
  end do

end subroutine sfluxv