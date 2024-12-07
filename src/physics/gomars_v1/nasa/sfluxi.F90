subroutine sfluxi( &
  plev, tlev, &
  dtaui, taucumi, taugsurf, &
  albi, wbari, cosbi, &
  fluxupi, fluxdni, fmneti, nfluxtopi)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_rad_mod

  implicit none

  real(r8), intent(in ) :: plev     (2*nlev+3)
  real(r8), intent(in ) :: tlev     (2*nlev+3)
  real(r8), intent(in ) :: dtaui    (nlayrad ,nspecti,ngauss  )
  real(r8), intent(in ) :: taucumi  (2*nlev+3,nspecti,ngauss  )
  real(r8), intent(in ) :: taugsurf (         nspecti,ngauss-1)
  real(r8), intent(in ) :: albi
  real(r8), intent(in ) :: wbari    (nlayrad ,nspecti,ngauss  )
  real(r8), intent(in ) :: cosbi    (nlayrad ,nspecti,ngauss  )
  real(r8), intent(out) :: fluxupi  (nlayrad)
  real(r8), intent(out) :: fluxdni  (nlayrad)
  real(r8), intent(out) :: fmneti   (nlayrad)
  real(r8), intent(out) :: nfluxtopi

  integer k, l, is, ig
  real(r8) ttop
  real(r8) tsfc
  real(r8) nts
  real(r8) ntt
  real(r8) btop
  real(r8) bsfc
  real(r8) fzero
  real(r8) ftopup
  real(r8) fmupi(nlevrad)
  real(r8) fmdni(nlevrad)

  nfluxtopi = 0
  fluxupi   = 0
  fluxdni   = 0
  fmneti    = 0

  ttop = tlev(2)
  tsfc = tlev(2*nlev+3)

  ntt = ttop * 10 - 499
  nts = tsfc * 10 - 499

  do is = 1, nspecti
    bsfc = (1 - albi) * planckir(is,nts)
    fzero = fzeroi(is)
    if (fzero < 0.99) then
      do ig = 1, ngauss - 1
        if (taugsurf(is,ig) < tlimiti) then
          fzero = fzero + (1 - fzeroi(is)) * gweight(ig)
        else
          btop = (1 - exp(-dtaui(1,is,ig) * plev(2) / (plev(4) - plev(2)) / ubari)) * planckir(is,ntt)
          ! call gfluxi
          nfluxtopi = nfluxtopi + ftopup * dwni(is) * gweight(ig) * (1 - fzeroi(is))
          do l = 1, nlevrad - 1
            fluxupi(l) = fluxupi(l) + fmupi(l) * dwni(is) * gweight(ig) * (1 - fzeroi(is))
            fluxdni(l) = fluxdni(l) + fmdni(l) * dwni(is) * gweight(ig) * (1 - fzeroi(is))
            fmneti (l) = fmneti (l) + (fmupi(l) - fmdni(l)) * dwni(is) * gweight(ig) * (1 - fzeroi(is))
          end do
        end if
      end do
    else
      ! Special 17th Gauss-point
      ig = ngauss
      btop = (1 - exp(-dtaui(1,is,ig) * plev(2) / (plev(4) - plev(2)) / ubari)) * planckir(is,ntt)
      ! call gfluxi
      nfluxtopi = nfluxtopi + ftopup * dwni(is) * fzero
      do l = 1, nlevrad - 1
        fluxupi(l) = fluxupi(l) + fmupi(l) * dwni(is) * fzero
        fluxdni(l) = fluxdni(l) + fmdni(l) * dwni(is) * fzero
        fmneti (l) = fmneti (l) + (fmupi(l) - fmdni(l)) * dwni(is) * fzero
      end do
    end if
  end do

end subroutine sfluxi