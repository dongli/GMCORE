subroutine optci( &
  plev, pmid, tmid, qh2o, &
  qxidst, qsidst, gidst, qextrefdst, &
  qxicld, qsicld, gicld, qextrefcld, &
  wbari, cosbi, dtaui, taui, taucumi, &
  taugsurf, taurefdst, taurefcld)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_rad_mod

  implicit none

  real(r8), intent(in   ) :: plev      (2*nlev+3)
  real(r8), intent(in   ) :: pmid      (2*nlev+3)
  real(r8), intent(in   ) :: tmid      (2*nlev+3)
  real(r8), intent(in   ) :: qh2o      (2*nlev+3)
  real(r8), intent(in   ) :: qxidst    (2*nlev+4,nspecti)
  real(r8), intent(in   ) :: qsidst    (2*nlev+4,nspecti)
  real(r8), intent(in   ) :: gidst     (2*nlev+4,nspecti)
  real(r8), intent(in   ) :: qextrefdst(2*nlev+4)
  real(r8), intent(in   ) :: qxicld    (2*nlev+4,nspecti)
  real(r8), intent(in   ) :: qsicld    (2*nlev+4,nspecti)
  real(r8), intent(in   ) :: gicld     (2*nlev+4,nspecti)
  real(r8), intent(in   ) :: qextrefcld(2*nlev+4)
  real(r8), intent(  out) :: wbari     (nlayrad ,nspecti,ngauss)
  real(r8), intent(  out) :: cosbi     (nlayrad ,nspecti,ngauss)
  real(r8), intent(  out) :: dtaui     (nlayrad ,nspecti,ngauss)
  real(r8), intent(  out) :: taui      (nlevrad ,nspecti,ngauss)
  real(r8), intent(  out) :: taucumi   (2*nlev+3,nspecti,ngauss)
  real(r8), intent(  out) :: taugsurf  (         nspecti,ngauss-1)
  real(r8), intent(inout) :: taurefdst (2*nlev+4)
  real(r8), intent(inout) :: taurefcld (2*nlev+4)


  integer n, k, l, is, ig
  integer idx_t          (2*nlev+3)
  integer idx_p          (2*nlev+3)
  integer idx_h2o        (2*nlev+3)
  real(r8) dtauki        (2*nlev+4,nspecti,ngauss)
  real(r8) taureflk      (2*nlev+4,nspecti)
  real(r8) taucldk       (2*nlev+4,nspecti)
  real(r8) taurefdst_save(2*nlev+4)
  real(r8) taurefcld_save(2*nlev+4)
  real(r8) wratio        (2*nlev+3)
  real(r8) lcoef       (4,2*nlev+3)
  real(r8) tdst          (2*nlev+3,nspecti)
  real(r8) tcld          (2*nlev+3,nspecti)
  real(r8) ans
  real(r8) kcoef(4)
  real(r8) tauac

  n = 2 * nlev + 3

  do is = 1, nspecti
    do ig = 1, ngauss
      dtauki(n+1,is,ig) = 0
    end do
    taureflk(n+1,is) = 0
    taucldk (n+1,is) = 0
  end do

  taurefdst_save = taurefdst
  taurefcld_save = taurefcld

  do ig = 1, ngauss - 1
    do is = 1, nspecti
      taugsurf(is,ig) = 0
    end do
  end do

  do k = 2, n
    call tpindex(pmid(k), tmid(k), qh2o(k), lcoef(:,k), idx_t(k), idx_p(k), idx_h2o(k), wratio(k))
    taurefdst(k) = taurefdst(k) / qextrefdst(k)
    taurefcld(k) = taurefcld(k) / qextrefcld(k)
    do is = 1, nspecti
      tdst(k,is) = taurefdst(k) * qxidst(k,is)
      tcld(k,is) = taurefcld(k) * qxicld(k,is)
    end do
  end do

  do k = 2, n
    do is = 1, nspecti
      do ig = 1, ngauss - 1
        kcoef(1) = co2i(idx_t(k)  ,idx_p(k)  ,idx_h2o(k)  ,is,ig) + wratio(k) * &
                  (co2i(idx_t(k)  ,idx_p(k)  ,idx_h2o(k)+1,is,ig) - &
                   co2i(idx_t(k)  ,idx_p(k)  ,idx_h2o(k)  ,is,ig))
        kcoef(2) = co2i(idx_t(k)  ,idx_p(k)+1,idx_h2o(k)  ,is,ig) + wratio(k) * &
                  (co2i(idx_t(k)  ,idx_p(k)+1,idx_h2o(k)+1,is,ig) - &
                   co2i(idx_t(k)  ,idx_p(k)+1,idx_h2o(k)  ,is,ig))
        kcoef(3) = co2i(idx_t(k)+1,idx_p(k)+1,idx_h2o(k)  ,is,ig) + wratio(k) * &
                  (co2i(idx_t(k)+1,idx_p(k)+1,idx_h2o(k)+1,is,ig) - &
                   co2i(idx_t(k)+1,idx_p(k)+1,idx_h2o(k)  ,is,ig))
        kcoef(4) = co2i(idx_t(k)+1,idx_p(k)  ,idx_h2o(k)  ,is,ig) + wratio(k) * &
                  (co2i(idx_t(k)+1,idx_p(k)  ,idx_h2o(k)+1,is,ig) - &
                   co2i(idx_t(k)+1,idx_p(k)  ,idx_h2o(k)  ,is,ig))
        ! Interpolate the CO2 k-coefficients to the requested temperature and pressure.
        ans = (lcoef(1,k) * kcoef(1) + lcoef(2,k) * kcoef(2) + &
               lcoef(3,k) * kcoef(3) + lcoef(4,k) * kcoef(4)) * cmk * (plev(k) - plev(k-1))
        taugsurf(is,ig) = taugsurf(is,ig) + ans
        dtauki(k,is,ig) = tdst(k,is) + tcld(k,is) + ans
      end do
      ! Now fill in the "clear" part of the spectrum (ig = ngauss), which holds
      ! continuum opacity only.
      dtauki(k,is,ngauss) = tdst(k,is) + tcld(k,is)
    end do
  end do

  do is = 1, nspecti
    do k = 2, n - 1
      taureflk(k,is) = taurefdst(k) * qsidst(k,is)
      taucldk (k,is) = taurefcld(k) * qsicld(k,is)
    end do
  end do

  do is = 1, nspecti
    ig = ngauss
    do l = 1, nlayrad
      k = 2 * l + 1
      dtaui(l,is,ig) = dtauki(k,is,ig) + dtauki(k+1,is,ig) + 1.0d-50
      if (dtaui(l,is,ig) > 1.0e-9_r8) then
        wbari(l,is,ig) = (taureflk(k,is) + taureflk(k+1,is) + taucldk(k,is) + taucldk(k+1,is)) / dtaui(l,is,ig)
      else
        wbari(l,is,ig) = 0
        dtaui(l,is,ig) = 1.0e-9_r8
      end if
      tauac = taureflk(k,is) + taureflk(k+1,is) + taucldk(k,is) + taucldk(k+1,is)
      if (tauac > 0) then
        cosbi(l,is,ig) = (gidst(k,is) * taureflk(k,is) + gidst(k+1,is) * taureflk(k+1,is)  + &
                          gicld(k,is) * taucldk (k,is) + gicld(k+1,is) * taucldk (k+1,is)) / &
                         (taureflk(k,is) + taureflk(k+1,is) + taucldk(k,is) + taucldk(k+1,is))
      else
        cosbi(l,is,ig) = 0
      end if
    end do
    do ig = 1, ngauss - 1
      do l = 1, nlayrad
        k = 2 * l + 1
        dtaui(l,is,ig) = dtauki(k,is,ig) + dtauki(k+1,is,ig) + 1.0d-50
        if (dtaui(l,is,ig) > 1.0e-9_r8) then
          wbari(l,is,ig) = (taureflk(k,is) + taureflk(k+1,is) + taucldk(k,is) + taucldk(k+1,is)) / dtaui(l,is,ig)
        else
          wbari(l,is,ig) = 0
          dtaui(l,is,ig) = 1.0e-9_r8
        end if
        cosbi(l,is,ig) = cosbi(l,is,ngauss)
      end do
    end do
  end do

  do is = 1, nspecti
    ig = ngauss
    taui(1,is,ig) = 0
    do l = 1, nlayrad
      taui(l+1,is,ig) = taui(l,is,ig) + dtaui(l,is,ig)
    end do
    taucumi(1,is,ig) = 0
    do k = 2, n
      taucumi(k,is,ig) = taucumi(k-1,is,ig) + dtaui(k,is,ig)
    end do
    do ig = 1, ngauss - 1
      taui(1,is,ig) = 0
      do l = 1, nlayrad
        taui(l+1,is,ig) = taui(l,is,ig) + dtaui(l,is,ig)
      end do
      taucumi(1,is,ig) = 0
      do k = 2, n
        taucumi(k,is,ig) = taucumi(k-1,is,ig) + dtaui(k,is,ig)
      end do
    end do
  end do

  taurefdst = taurefdst_save
  taurefcld = taurefcld_save

end subroutine optci