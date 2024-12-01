subroutine sfluxv(dtauv, tauv, taucumv, taugsurf, sol, cosz, als, wbarv, cosbv, &
                  fluxupv, fluxdnv, fmnetv, nfluxtopv, diffvt, detau)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_rad_mod

  implicit none

  real(r8), intent(in ) :: dtauv   (nlayrad ,nspectv,ngauss  )
  real(r8), intent(in ) :: tauv    (nlevrad ,nspectv,ngauss  )
  real(r8), intent(in ) :: taucumv (2*nlev+3,nspectv,ngauss  )
  real(r8), intent(in ) :: taugsurf(         nspectv,ngauss-1)
  real(r8), intent(in ) :: sol     (         nspectv         )
  real(r8), intent(in ) :: cosz
  real(r8), intent(in ) :: als
  real(r8), intent(in ) :: wbarv   (nlayrad ,nspectv,ngauss  )
  real(r8), intent(in ) :: cosbv   (nlayrad ,nspectv,ngauss  )
  real(r8), intent(out) :: fluxupv (nlayrad)
  real(r8), intent(out) :: fluxdnv (nlayrad)
  real(r8), intent(out) :: fmnetv  (nlayrad)
  real(r8), intent(out) :: nfluxtopv
  real(r8), intent(out) :: diffvt
  real(r8), intent(out) :: detau   (         nspectv,ngauss  )

  integer n, k, l, is, ig
  real(r8) fzero
  real(r8) btop
  real(r8) bsfc
  real(r8) fluxup
  real(r8) fluxdn
  real(r8) fmupv(nlayrad)
  real(r8) fmdnv(nlayrad)
  real(r8) diffv

  n = 2 * nlev + 3

  fluxupv   = 0
  fluxdnv   = 0
  fmnetv    = 0
  nfluxtopv = 0
  diffvt    = 0

  ! Calculate the net flux in each spectral interval.
  do is = 1, nspectv
    fzero = fzerov(is)
    if (fzerov(is) < 0.99_r8) then
      do ig = 1, ngauss - 1
        call getdetau(    &
          dtauv  (1,is,ig), &
          tauv   (1,is,ig), &
          taucumv(1,is,ig), &
          wbarv  (1,is,ig), &
          cosbv  (1,is,ig), &
          detau  (  is,ig)  &
        )
        if (taugsurf(is,ig) < tlimits) then
          fzero = fzero + (1 - fzerov(is)) * gweight(ig)
        else
          ! Set up the top and bottom boundary conditions.
          btop = 0
          bsfc = als * cosz * sol(is) * exp(-min(detau(is,ig) / cosz, maxexp))
          ! We can now solve for the coefficients of the two-stream problem.
          call gfluxv( &
            dtauv  (1,is,ig), &
            tauv   (1,is,ig), &
            taucumv(1,is,ig), &
            wbarv  (1,is,ig), &
            cosbv  (1,is,ig), &
            cosz          , &
            sol           , &
            als           , &
            btop          , &
            bsfc          , &
            fmupv         , &
            fmdnv         , &
            diffv         , &
            fluxup        , &
            fluxdn        , &
            detau    (is,ig)  &
          )
          ! Calculate the cumulative visible net flux.
          nfluxtopv = nfluxtopv + (fluxup - fluxdn) * gweight(ig) * (1 - fzerov(is))
          do l = 1, nlayrad
            fmnetv (l) = fmnetv (l) + (fmupv(l) - fmdnv(l)) * gweight(ig) * (1 - fzerov(is))
            fluxupv(l) = fluxupv(l) + fmupv(l) * gweight(ig) * (1 - fzerov(is))
            fluxdnv(l) = fluxdnv(l) + fmdnv(l) * gweight(ig) * (1 - fzerov(is))
          end do
          ! The diffuse component of the downward solar flux.
          diffvt = diffvt + diffv * gweight(ig) * (1 - fzerov(is))
        end if
      end do
    else
      ! Special 17th Gauss point
      ig = ngauss
      call getdetau(    &
        dtauv  (1,is,ig), &
        tauv   (1,is,ig), &
        taucumv(1,is,ig), &
        wbarv  (1,is,ig), &
        cosbv  (1,is,ig), &
        detau  (  is,ig)  &
      )
      btop = 0
      bsfc = als * cosz * sol(is) * exp(-min(detau(is,ig) / cosz, maxexp))
      call gfluxv(      &
        dtauv  (1,is,ig), &
        tauv   (1,is,ig), &
        taucumv(1,is,ig), &
        wbarv  (1,is,ig), &
        cosbv  (1,is,ig), &
        cosz          , &
        sol           , &
        als           , &
        btop          , &
        bsfc          , &
        fmupv         , &
        fmdnv         , &
        diffv         , &
        fluxup        , &
        fluxdn        , &
        detau    (is,ig)  &
      )
      nfluxtopv = nfluxtopv + (fluxup - fluxdn) * gweight(ig) * fzero
      do l = 1, nlayrad
        fmnetv (l) = fmnetv (l) + (fmupv(l) - fmdnv(l)) * gweight(ig) * fzero
        fluxupv(l) = fluxupv(l) + fmupv(l) * gweight(ig) * fzero
        fluxdnv(l) = fluxdnv(l) + fmdnv(l) * gweight(ig) * fzero
      end do
      diffvt = diffvt + diffv * gweight(ig) * fzero
    end if
  end do

end subroutine sfluxv