subroutine optcv( &
  plev          , &
  pmid          , &
  tmid          , &
  qh2o          , &
  qxvdst        , &
  qsvdst        , &
  gvdst         , &
  qxvcld        , &
  qsvcld        , &
  gvcld         , &
  qextrefcld    , &
  wbarv         , &
  cosbv         , &
  dtauv         , &
  tauv          , &
  taucumv       , &
  taugsurf      , &
  taurefdst     , &
  taurefcld     )

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
  real(r8), intent(in   ) :: qxvdst    (2*nlev+4,nspectv)
  real(r8), intent(in   ) :: qsvdst    (2*nlev+4,nspectv)
  real(r8), intent(in   ) :: gvdst     (2*nlev+4,nspectv)
  real(r8), intent(in   ) :: qxvcld    (2*nlev+4,nspectv)
  real(r8), intent(in   ) :: qsvcld    (2*nlev+4,nspectv)
  real(r8), intent(in   ) :: gvcld     (2*nlev+4,nspectv)
  real(r8), intent(in   ) :: qextrefcld(2*nlev+4)
  real(r8), intent(  out) :: wbarv     (nlayrad ,nspectv,ngauss  )
  real(r8), intent(  out) :: cosbv     (nlayrad ,nspectv,ngauss  )
  real(r8), intent(  out) :: dtauv     (nlayrad ,nspectv,ngauss  )
  real(r8), intent(  out) :: tauv      (nlevrad ,nspectv,ngauss  )
  real(r8), intent(  out) :: taucumv   (2*nlev+3,nspectv,ngauss  )
  real(r8), intent(  out) :: taugsurf  (         nspectv,ngauss-1)
  real(r8), intent(inout) :: taurefdst (2*nlev+4)
  real(r8), intent(inout) :: taurefcld (2*nlev+4)

  integer n, k, l, is, ig
  integer idx_t  (  2*nlev+3)
  integer idx_p  (  2*nlev+3)
  integer idx_h2o(  2*nlev+3)
  real(r8) wratio(  2*nlev+3)
  real(r8) lcoef (4,2*nlev+3)
  real(r8) taurefdst_save(2*nlev+4)
  real(r8) taurefcld_save(2*nlev+4)
  real(r8) tray  (2*nlev+3,nspectv)
  real(r8) tdst  (2*nlev+3,nspectv)
  real(r8) tcld  (2*nlev+3,nspectv)
  real(r8) dtaukv(2*nlev+4,nspectv,ngauss)
  real(r8) ans
  real(r8) trayaer ! Tau Rayleigh scattering plus aerosol opacity
  real(r8) kcoef(4)

  n = 2 * nlev + 3

  taurefdst_save = taurefdst
  taurefcld_save = taurefcld

  ! Determine the total gas opacity throughout the column.
  do ig = 1, ngauss - 1
    do is = 1, nspectv
      taugsurf(is,ig) = 0
    end do
  end do

  do k = 2, n
    call tpindex(pmid(k), tmid(k), qh2o(k), lcoef(:,k), idx_t(k), idx_p(k), idx_h2o(k), wratio(k))
    taurefdst(k) = taurefdst(k) / qxvdst(k,nrefv)
    taurefcld(k) = taurefcld(k) / qextrefcld(k)
    do is = 1, nspectv
      tray(k,is) = tauray(is) * (plev(k) - plev(k-1))
      tdst(k,is) = taurefdst(k) * qxvdst(k,is)
      tcld(k,is) = taurefcld(k) * qxvcld(k,is)
    end do
  end do

  do k = 2, n
    do is = 1, nspectv
      trayaer = tray(k,is) + tdst(k,is) + tcld(k,is)
      do ig = 1, ngauss - 1
        ! Interpolate between water mixing ratios.
        ! wratio is zero if the requested water amount is equal to, or outside
        ! the range of the reference values.
        kcoef(1) = co2v(idx_t(k)  ,idx_p(k)  ,idx_h2o(k)  ,is,ig) + wratio(k) * &
                  (co2v(idx_t(k)  ,idx_p(k)  ,idx_h2o(k)+1,is,ig) - &
                   co2v(idx_t(k)  ,idx_p(k)  ,idx_h2o(k)  ,is,ig))
        kcoef(2) = co2v(idx_t(k)  ,idx_p(k)+1,idx_h2o(k)  ,is,ig) + wratio(k) * &
                  (co2v(idx_t(k)  ,idx_p(k)+1,idx_h2o(k)+1,is,ig) - &
                   co2v(idx_t(k)  ,idx_p(k)+1,idx_h2o(k)  ,is,ig))
        kcoef(3) = co2v(idx_t(k)+1,idx_p(k)+1,idx_h2o(k)  ,is,ig) + wratio(k) * &
                  (co2v(idx_t(k)+1,idx_p(k)+1,idx_h2o(k)+1,is,ig) - &
                   co2v(idx_t(k)+1,idx_p(k)+1,idx_h2o(k)  ,is,ig))
        kcoef(4) = co2v(idx_t(k)+1,idx_p(k)  ,idx_h2o(k)  ,is,ig) + wratio(k) * &
                  (co2v(idx_t(k)+1,idx_p(k)  ,idx_h2o(k)+1,is,ig) - &
                   co2v(idx_t(k)+1,idx_p(k)  ,idx_h2o(k)  ,is,ig))
        ! Interpolate the CO2 k-coefficients to the requested temperature and pressure.
        ans = (lcoef(1,k) * kcoef(1) + lcoef(2,k) * kcoef(2) + &
               lcoef(3,k) * kcoef(3) + lcoef(4,k) * kcoef(4)) * cmk * (plev(k) - plev(k-1))
        taugsurf(is,ig) = taugsurf(is,ig) + ans
        dtaukv(k,is,ig) = trayaer + ans
      end do
      ! Now fill in the "clear" part of the spectrum (ig = ngauss), which holds
      ! continuum opacity only.
      dtaukv(k,is,ngauss) = trayaer
    end do
  end do

  do is = 1, nspectv
    ! First, the special "clear" channel
    ig = ngauss
    do l = 1, nlayrad - 1
      k = 2 * l + 1
      dtauv(l,is,ig) = dtaukv(k,is,ig) + dtaukv(k+1,is,ig)
      cosbv(l,is,ig) = (gvdst(k  ,is) * taurefdst(k  ) * qsvdst(k  ,is)  + &
                        gvdst(k+1,is) * taurefdst(k+1) * qsvdst(k+1,is)  + &
                        gvcld(k  ,is) * taurefcld(k  ) * qsvcld(k  ,is)  + &
                        gvcld(k+1,is) * taurefcld(k+1) * qsvcld(k+1,is)) / &
                       (tray(k,is) + tray(k+1,is) + &
                        taurefdst(k  ) * qsvdst(k  ,is) + &
                        taurefdst(k+1) * qsvdst(k+1,is) + &
                        taurefcld(k  ) * qsvcld(k  ,is) + &
                        taurefcld(k+1) * qsvcld(k+1,is))
      wbarv(l,is,ig) = (taurefdst(k) * qsvdst(k,is) + taurefdst(k+1) * qsvdst(k+1,is) + &
                        taurefcld(k) * qsvcld(k,is) + taurefcld(k+1) * qsvcld(k+1,is) + &
                        (tray(k,is) + tray(k+1,is)) * 0.9999_r8) / dtauv(l,is,ig)
    end do
    ! Special bottom layer
    l = nlayrad
    k = 2 * l + 1
    dtauv(l,is,ig) = dtaukv(k,is,ig)
    cosbv(l,is,ig) = (gvdst(k,is) * taurefdst(k) * qsvdst(k,is)  + &
                      gvcld(k,is) * taurefcld(k) * qsvcld(k,is)) / &
                     (tray(k,is) + taurefdst(k) * qsvdst(k,is) + taurefcld(k) * qsvcld(k,is))
    wbarv(l,is,ig) = (taurefdst(k) * qsvdst(k,is) + taurefcld(k) * qsvcld(k,is) + &
                      tray(k,is) * 0.9999_r8) / dtauv(l,is,ig)
    ! Now the other Gauss points, if needed
    do ig = 1, ngauss - 1
      do l = 1, nlayrad - 1
        k = 2 * l + 1
        dtauv(l,is,ig) = dtaukv(k,is,ig) + dtaukv(k+1,is,ig)
        cosbv(l,is,ig) = cosbv(l,is,ngauss)
        wbarv(l,is,ig) = (taurefdst(k) * qsvdst(k,is) + taurefdst(k+1) * qsvdst(k+1,is) + &
                          taurefcld(k) * qsvcld(k,is) + taurefcld(k+1) * qsvcld(k+1,is) + &
                          (tray(k,is) + tray(k+1,is)) * 0.9999_r8) / dtauv(l,is,ig)
      end do
      ! Special bottom layer
      l = nlayrad
      k = 2 * l + 1
      dtauv(l,is,ig) = dtaukv(k,is,ig)
      cosbv(l,is,ig) = cosbv(l,is,ngauss)
      wbarv(l,is,ig) = (taurefdst(k) * qsvdst(k,is) + taurefcld(k) * qsvcld(k,is) + &
                        tray(k,is) * 0.9999_r8) / dtauv(l,is,ig)
    end do
  end do

  ! Calculate the total extinction optical depths.
  do is = 1, nspectv
    ig = ngauss
    tauv(1,is,ig) = 0
    do l = 1, nlayrad
      tauv(l+1,is,ig) = tauv(l,is,ig) + dtauv(l,is,ig)
    end do
    taucumv(1,is,ig) = 0
    do k = 2, n
      taucumv(k,is,ig) = taucumv(k-1,is,ig) + dtaukv(k,is,ig)
    end do
    do ig = 1, ngauss - 1
      tauv(1,is,ig) = 0
      do l = 1, nlayrad
        tauv(l+1,is,ig) = tauv(l,is,ig) + dtauv(l,is,ig)
      end do
      taucumv(1,is,ig) = 0
      do k = 2, n
        taucumv(k,is,ig) = taucumv(k-1,is,ig) + dtaukv(k,is,ig)
      end do
    end do
  end do

  ! Restore the original values of tauref and taurefcld.
  taurefdst = taurefdst_save
  taurefcld = taurefcld_save

end subroutine optcv