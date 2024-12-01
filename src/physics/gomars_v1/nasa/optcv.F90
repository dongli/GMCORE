subroutine optcv( &
  pl, tl, qh2o, &
  qxvdst, qsvdst, gvdst, &
  qxvcld, qsvcld, gvcld, qextrefcld, &
  wbarv, cosbv, dtauv, tauv, taucumv, &
  taugsurf, taurefdst, taurefcld)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod, only: r8, nlev, nspectv, ngauss, nlayrad, nlevrad, cmk, nrefv
  use gomars_v1_rad_mod, only: co2v, tgasref, pfgasref, tauray

  implicit none

  real(r8), intent(in   ) :: pl        (2*nlev+3)
  real(r8), intent(in   ) :: tl        (2*nlev+3)
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

  integer n, k, l, s, g
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
  real(r8) dp, ans
  real(r8) trayaer ! Tau Rayleigh scattering plus aerosol opacity
  real(r8) kcoef(4)

  n = 2 * nlev + 3

  taurefdst_save = taurefdst
  taurefcld_save = taurefcld

  ! Determine the total gas opacity throughout the column.
  do g = 1, ngauss - 1
    do s = 1, nspectv
      taugsurf(s,g) = 0
    end do
  end do

  do k = 2, n
    call tpindex(tl(k), pl(k), qh2o(k), lcoef(:,k), idx_t(k), idx_p(k), idx_h2o(k), wratio(k))
    taurefdst(k) = taurefdst(k) / qxvdst(k,nrefv)
    taurefcld(k) = taurefcld(k) / qextrefcld(k)
    do s = 1, nspectv
      tray(k,s) = tauray(s) * (pl(k) - pl(k-1))
      tdst(k,s) = taurefdst(k) * qxvdst(k,s)
      tcld(k,s) = taurefcld(k) * qxvcld(k,s)
    end do
  end do

  do k = 2, n
    do s = 1, nspectv
      trayaer = tray(k,s) + tdst(k,s) + tcld(k,s)
      do g = 1, ngauss - 1
        ! Interpolate between water mixing ratios.
        ! wratio is zero if the requested water amount is equal to, or outside
        ! the range of the reference values.
        kcoef(1) = co2v(idx_t(k)  ,idx_p(k)  ,idx_h2o(k)  ,s,g) + wratio(k) * &
                  (co2v(idx_t(k)  ,idx_p(k)  ,idx_h2o(k)+1,s,g) - &
                   co2v(idx_t(k)  ,idx_p(k)  ,idx_h2o(k)  ,s,g))
        kcoef(2) = co2v(idx_t(k)  ,idx_p(k)+1,idx_h2o(k)  ,s,g) + wratio(k) * &
                  (co2v(idx_t(k)  ,idx_p(k)+1,idx_h2o(k)+1,s,g) - &
                   co2v(idx_t(k)  ,idx_p(k)+1,idx_h2o(k)  ,s,g))
        kcoef(3) = co2v(idx_t(k)+1,idx_p(k)+1,idx_h2o(k)  ,s,g) + wratio(k) * &
                  (co2v(idx_t(k)+1,idx_p(k)+1,idx_h2o(k)+1,s,g) - &
                   co2v(idx_t(k)+1,idx_p(k)+1,idx_h2o(k)  ,s,g))
        kcoef(4) = co2v(idx_t(k)+1,idx_p(k)  ,idx_h2o(k)  ,s,g) + wratio(k) * &
                  (co2v(idx_t(k)+1,idx_p(k)  ,idx_h2o(k)+1,s,g) - &
                   co2v(idx_t(k)+1,idx_p(k)  ,idx_h2o(k)  ,s,g))
        ! Interpolate the CO2 k-coefficients to the requested temperature and pressure.
        ans = (lcoef(1,k) * kcoef(1) + lcoef(2,k) * kcoef(2) + &
               lcoef(3,k) * kcoef(3) + lcoef(4,k) * kcoef(4)) * cmk * (pl(k) - pl(k-1))
        taugsurf(s,g) = taugsurf(s,g) + ans
        dtaukv(k,s,g) = trayaer + ans
      end do
      ! Now fill in the "clear" part of the spectrum (g = ngauss), which holds
      ! continuum opacity only.
      dtaukv(k,s,ngauss) = trayaer
    end do
  end do

  do s = 1, nspectv
    ! First, the special "clear" channel
    g = ngauss
    do l = 1, nlayrad - 1
      k = 2 * l + 1
      dtauv(l,s,g) = dtaukv(k,s,g) + dtaukv(k+1,s,g)
      cosbv(l,s,g) = (gvdst(k  ,s) * taurefdst(k  ) * qsvdst(k  ,s)  + &
                      gvdst(k+1,s) * taurefdst(k+1) * qsvdst(k+1,s)  + &
                      gvcld(k  ,s) * taurefcld(k  ) * qsvcld(k  ,s)  + &
                      gvcld(k+1,s) * taurefcld(k+1) * qsvcld(k+1,s)) / &
                     (tray(k,s) + tray(k+1,s) + &
                      taurefdst(k  ) * qsvdst(k  ,s) + &
                      taurefdst(k+1) * qsvdst(k+1,s) + &
                      taurefcld(k  ) * qsvcld(k  ,s) + &
                      taurefcld(k+1) * qsvcld(k+1,s))
      wbarv(l,s,g) = (taurefdst(k) * qsvdst(k,s) + taurefdst(k+1) * qsvdst(k+1,s) + &
                      taurefcld(k) * qsvcld(k,s) + taurefcld(k+1) * qsvcld(k+1,s) + &
                      (tray(k,s) + tray(k+1,s)) * 0.9999_r8) / dtauv(l,s,g)
    end do
    ! Special bottom layer
    l = nlayrad
    k = 2 * l + 1
    dtauv(l,s,g) = dtaukv(k,s,g)
    cosbv(l,s,g) = (gvdst(k,s) * taurefdst(k) * qsvdst(k,s)  + &
                    gvcld(k,s) * taurefcld(k) * qsvcld(k,s)) / &
                   (tray(k,s) + taurefdst(k) * qsvdst(k,s) + taurefcld(k) * qsvcld(k,s))
    wbarv(l,s,g) = (taurefdst(k) * qsvdst(k,s) + taurefcld(k) * qsvcld(k,s) + &
                    tray(k,s) * 0.9999_r8) / dtauv(l,s,g)
    ! Now the other Gauss points, if needed
    do g = 1, ngauss - 1
      do l = 1, nlayrad - 1
        k = 2 * l + 1
        dtauv(l,s,g) = dtaukv(k,s,g) + dtaukv(k+1,s,g)
        cosbv(l,s,g) = cosbv(l,s,ngauss)
        wbarv(l,s,g) = (taurefdst(k) * qsvdst(k,s) + taurefdst(k+1) * qsvdst(k+1,s) + &
                        taurefcld(k) * qsvcld(k,s) + taurefcld(k+1) * qsvcld(k+1,s) + &
                        (tray(k,s) + tray(k+1,s)) * 0.9999_r8) / dtauv(l,s,g)
      end do
      ! Special bottom layer
      l = nlayrad
      k = 2 * l + 1
      dtauv(l,s,g) = dtaukv(k,s,g)
      cosbv(l,s,g) = cosbv(l,s,ngauss)
      wbarv(l,s,g) = (taurefdst(k) * qsvdst(k,s) + taurefcld(k) * qsvcld(k,s) + &
                      tray(k,s) * 0.9999_r8) / dtauv(l,s,g)
    end do
  end do

  ! Calculate the total extinction optical depths.
  do s = 1, nspectv
    g = ngauss
    tauv(1,s,g) = 0
    do l = 1, nlayrad
      tauv(l+1,s,g) = tauv(l,s,g) + dtauv(l,s,g)
    end do
    taucumv(1,s,g) = 0
    do k = 2, n
      taucumv(k,s,g) = taucumv(k-1,s,g) + dtaukv(k,s,g)
    end do
    do g = 1, ngauss - 1
      tauv(1,s,g) = 0
      do l = 1, nlayrad
        tauv(l+1,s,g) = tauv(l,s,g) + dtauv(l,s,g)
      end do
      taucumv(1,s,g) = 0
      do k = 2, n
        taucumv(k,s,g) = taucumv(k-1,s,g) + dtaukv(k,s,g)
      end do
    end do
  end do

  ! Restore the original values of tauref and taurefcld.
  taurefdst = taurefdst_save
  taurefcld = taurefcld_save

end subroutine optcv