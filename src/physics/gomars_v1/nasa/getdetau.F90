subroutine getdetau(dtdel, tdel, taucumin, wdel, cdel, detau)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod

  implicit none

  real(r8), intent(in ) :: dtdel    (nlayrad)
  real(r8), intent(in ) :: tdel     (nlevrad)
  real(r8), intent(in ) :: taucumin (2*nlev+3)
  real(r8), intent(in ) :: wdel     (nlayrad)
  real(r8), intent(in ) :: cdel     (nlayrad)
  real(r8), intent(out) :: detau

  integer k, l
  real(r8) factor
  real(r8) tau   (nlevrad)
  real(r8) taucum(2*nlev+3)

  ! Delta-Eddington scaling
  factor    = 1 - wdel(1) * cdel(1)**2
  tau   (1) = tdel(1) * factor
  taucum(1) = 0
  taucum(2) = taucumin(2) * factor
  taucum(3) = taucum(2) + (taucumin(3) - taucumin(2)) * factor

  do l = 1, nlayrad - 1
    factor      = 1 - wdel(l) * cdel(l)**2
    tau   (l+1) = tau(l) + dtdel(l) * factor
    k           = 2 * l + 2
    taucum(k  ) = tau(l+1)
    taucum(k+1) = taucum(k) + (taucumin(k+1) - taucumin(k)) * factor
  end do

  l         = nlayrad
  factor    = 1 - wdel(l) * cdel(l)**2
  tau(l+1)  = tau(l) + dtdel(l) * factor
  k         = 2 * l + 1
  taucum(k) = tau(l+1)
  detau     = taucum(k)

end subroutine getdetau