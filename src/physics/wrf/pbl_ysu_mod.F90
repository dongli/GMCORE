module pbl_ysu_mod

  !---------------------------------------------------------------------------
  !
  ! This code is a revised vertical diffusion package ("ysupbl") with a
  ! nonlocal turbulent mixing in the PBL after "mrfpbl". The ysupbl (Hong et
  ! al. 2006) is based on the study of Noh et al. (2003) and accumulated
  ! realism of the behavior of Troen and Mahrt (1986) concept implemented by
  ! Hong and Pan (1996). The major ingredient of the ysupbl is the inclusion
  ! of an explicit treatment of the entrainment processes at the entrainment
  ! layer. This routine uses an implicit approach for vertical flux divergence
  ! and does not require "miter" timesteps. It includes vertical diffusion in
  ! the stable atmosphere and moist vertical diffusion in clouds.
  !
  ! mrfpbl:
  !
  !   Coded by Song-You Hong (NCEP), implemented by Jimy Dudhia (NCAR)
  !   fall 1996
  !
  ! ysupbl:
  !
  !   Coded by Song-You Hong (Yonsei University) and implemented by
  !   Song-You Hong (Yonsei University) and Jimy Dudhia (NCAR)
  !   summer 2002
  !
  ! Further modifications :
  !
  !   An enhanced stable layer mixing, april 2008:
  !     ==> Increase PBL height when surface is stable (Hong 2010).
  !
  !   Pressure-level diffusion, april 2009:
  !     ==> Negligible differences.
  !
  !   Implicit forcing for momentum with clean up, july 2009:
  !     ==> Prevents model blowup when surface layer is too low.
  !
  !   Increase of lamda, maximum (30, 0.1 x del z) feb 2010:
  !     ==> Prevents model blowup when delz is extremely large.
  !
  !   Revised Prandtl number at surface, Peggy Lemone, feb 2010:
  !     ==> Increase kh, decrease mixing due to counter-gradient term.
  !
  !   Revised thermal, Shin et al. mon. wea. rev., Song-You Hong, aug 2011:
  !     ==> Reduce the thermal strength when z1 < 0.1 * h.
  !
  !   Revised Prandtl number for free convection, Dudhia, mar 2012:
  !     ==> pr0 = 1 + bke (=0.272) when neutral, kh is reduced.
  !
  !   Minimum kzo = 0.01, lo = min (30m,delz), Hong, mar 2012:
  !     ==> Weaker mixing when stable, and LES resolution in vertical
  !
  !   gz1oz0 is removed, and psim psih are ln(z1/z0)-psim,h, Hong, mar 2012:
  !     ==> Consider thermal z0 when differs from mechanical z0.
  !
  !   A bug fix in wscale computation in stable PBL, Sukanta Basu, jun 2012:
  !     ==> wscale becomes small with height, and less mixing in stable PBL.
  !
  !   Revision in background diffusion (kzo), jan 2016:
  !     ==> kzo = 0.1 for momentum and = 0.01 for mass to account for.
  !
  !    Internal wave mixing of Large et al. (1994), Song-You Hong, feb 2016:
  !     ==> Alleviate superious excessive mixing when delz is large
  !
  ! References:
  !
  !    Hong (2010) quart. j. roy. met. soc
  !    Hong, Noh, and Dudhia (2006), mon. wea. rev.
  !    Hong and Pan (1996), mon. wea. rev.
  !    Noh, Chun, Hong, and Raasch (2003), boundary layer met.
  !    Troen and Mahrt (1986), boundary layer met.
  !
  !---------------------------------------------------------------------------

  use const_mod, only: r8, cp => cpd, g, rd, rv, rd_o_rv, rv_o_rd, rovcp => rd_o_cpd, xlv => lv, karman => ka

  implicit none

  logical :: pbl_vdiff_qc = .true.
  logical :: pbl_vdiff_qi = .true.
  real(r8), parameter :: rcl      = 1.0_r8
  real(r8), parameter :: kzm_min  = 0.1_r8
  real(r8), parameter :: kzh_min  = 0.01_r8
  real(r8), parameter :: kz_min   = 0.01_r8
  real(r8), parameter :: kz_max   = 1000
  real(r8), parameter :: ri_min   = -100
  real(r8), parameter :: rlam     = 30
  real(r8), parameter :: prmin    = 0.25_r8
  real(r8), parameter :: prmax    = 4
  real(r8), parameter :: brcr_ub  = 0                ! Critical bulk Richardson number
  real(r8), parameter :: brcr_sb  = 0.25_r8
  real(r8), parameter :: cori     = 1.0e-4_r8
  real(r8), parameter :: afac     = 6.8_r8
  real(r8), parameter :: bfac     = 6.8_r8
  real(r8), parameter :: pfac     = 2
  real(r8), parameter :: pfac_q   = 2
  real(r8), parameter :: phifac   = 8
  real(r8), parameter :: sfcfrac  = 0.1_r8
  real(r8), parameter :: d1       = 0.02_r8
  real(r8), parameter :: d2       = 0.05_r8
  real(r8), parameter :: d3       = 0.001_r8
  real(r8), parameter :: h1       = 1.0_r8 / 3.0_r8
  real(r8), parameter :: h2       = 0.6666667_r8
  real(r8), parameter :: zf_min   = 1.0e-8_r8
  real(r8), parameter :: aphi5    = 5
  real(r8), parameter :: aphi16   = 16
  real(r8), parameter :: tmin     = 1.0e-2_r8
  real(r8), parameter :: gamcrt   = 3
  real(r8), parameter :: gamcrq   = 2.0e-3_r8
  real(r8), parameter :: xka      = 2.4e-5_r8
  integer , parameter :: imvdif   = 1

contains

  subroutine pbl_ysu_run( &
    ncol                , &
    nlev                , &
    dt                  , &
    u                   , &
    v                   , &
    t                   , &
    qv                  , &
    qc                  , &
    qi                  , &
    p                   , &
    p_lev               , &
    pk                  , &
    dz                  , &
    ps                  , &
    znt                 , &
    ustar               , &
    hpbl                , &
    psim                , &
    psih                , &
    xland               , &
    hfx                 , &
    qfx                 , &
    wspd                , &
    br                  , &
    kzm                 , &
    kzh                 , &
    dudt_pbl            , &
    dvdt_pbl            , &
    dtdt_pbl            , &
    dqvdt_pbl           , &
    dqcdt_pbl           , &
    dqidt_pbl           , &
    kpbl                , &
    wstar               , &
    delta               , &
    u10                 , &
    v10                 , &
    uos                 , &
    vos                 , &
    dptdt_rad           , &
    ysu_topdown_pblmix  , &
    ctopo               , &
    ctopo2              )

    integer , intent(in   ) :: ncol
    integer , intent(in   ) :: nlev
    real(r8), intent(in   ) :: dt
    real(r8), intent(in   ), dimension(ncol,nlev  ) :: u
    real(r8), intent(in   ), dimension(ncol,nlev  ) :: v
    real(r8), intent(in   ), dimension(ncol,nlev  ) :: t
    real(r8), intent(in   ), dimension(ncol,nlev  ) :: qv
    real(r8), intent(in   ), dimension(ncol,nlev  ) :: qc
    real(r8), intent(in   ), dimension(ncol,nlev  ) :: qi
    real(r8), intent(in   ), dimension(ncol,nlev  ) :: p
    real(r8), intent(in   ), dimension(ncol,nlev+1) :: p_lev
    real(r8), intent(in   ), dimension(ncol,nlev  ) :: pk
    real(r8), intent(in   ), dimension(ncol,nlev  ) :: dz
    real(r8), intent(in   ), dimension(ncol       ) :: ps
    real(r8), intent(inout), dimension(ncol       ) :: znt          ! Roughness length (m)
    real(r8), intent(inout), dimension(ncol       ) :: ustar        ! Friction velocity scale (m s-1)
    real(r8), intent(inout), dimension(ncol       ) :: hpbl         ! 这个和pblh是什么区别？
    real(r8), intent(in   ), dimension(ncol       ) :: psim         ! Stability function for momentum
    real(r8), intent(in   ), dimension(ncol       ) :: psih         ! Stability function for heat
    real(r8), intent(in   ), dimension(ncol       ) :: xland        ! Land mask (1 for land, 2 for water)
    real(r8), intent(in   ), dimension(ncol       ) :: hfx          ! Upward heat flux at surface (W m-1)
    real(r8), intent(in   ), dimension(ncol       ) :: qfx          ! upward moisture flux at surface (kg m-2 s-1)
    real(r8), intent(in   ), dimension(ncol       ) :: wspd         ! Wind speed at lowest model level (m s-1)
    real(r8), intent(in   ), dimension(ncol       ) :: br           ! Bulk Richardson number in surface layer
    real(r8), intent(inout), dimension(ncol,nlev  ) :: kzm          ! Exchange coefficients of momentum (m2 s-1)
    real(r8), intent(inout), dimension(ncol,nlev  ) :: kzh          ! Exchange coefficients of heat and scalars (K m s-1)
    real(r8), intent(inout), dimension(ncol,nlev  ) :: dudt_pbl
    real(r8), intent(inout), dimension(ncol,nlev  ) :: dvdt_pbl
    real(r8), intent(inout), dimension(ncol,nlev  ) :: dtdt_pbl
    real(r8), intent(inout), dimension(ncol,nlev  ) :: dqvdt_pbl
    real(r8), intent(inout), dimension(ncol,nlev  ) :: dqcdt_pbl
    real(r8), intent(inout), dimension(ncol,nlev  ) :: dqidt_pbl
    integer , intent(  out), dimension(ncol       ) :: kpbl
    real(r8), intent(  out), dimension(ncol       ) :: wstar
    real(r8), intent(  out), dimension(ncol       ) :: delta
    real(r8), intent(inout), dimension(ncol       ) :: u10
    real(r8), intent(inout), dimension(ncol       ) :: v10
    real(r8), intent(in   ), dimension(ncol       ) :: uos          ! Sea surface zonal current (m s-1)
    real(r8), intent(in   ), dimension(ncol       ) :: vos          ! Sea surface meridional current (m s-1)
    real(r8), intent(in   ), dimension(ncol,ncol  ) :: dptdt_rad    ! Tendency from radiation
    integer , intent(in   ) :: ysu_topdown_pblmix
    real(r8), intent(in   ), dimension(ncol       ), optional :: ctopo
    real(r8), intent(in   ), dimension(ncol       ), optional :: ctopo2

    real(r8), dimension(ncol,nlev+1) :: z_lev
    real(r8), dimension(ncol,nlev  ) :: z
    real(r8), dimension(ncol,nlev  ) :: th
    real(r8), dimension(ncol,nlev  ) :: thv
    real(r8), dimension(ncol,nlev  ) :: thli
    real(r8), dimension(ncol,nlev  ) :: dp
    real(r8), dimension(ncol,nlev  ) :: dz_lev
    real(r8), dimension(ncol,nlev  ) :: kzml
    real(r8), dimension(ncol,nlev  ) :: kzhl
    real(r8), dimension(ncol,nlev  ) :: kzq
    real(r8), dimension(ncol       ) :: rhos
    real(r8), dimension(ncol       ) :: govrth
    real(r8), dimension(ncol       ) :: thermal
    real(r8), dimension(ncol       ) :: ws0          ! Mixed-layer velocity scale (m s-1)
    real(r8), dimension(ncol       ) :: hgamt
    real(r8), dimension(ncol       ) :: hgamq
    real(r8), dimension(ncol       ) :: brdn
    real(r8), dimension(ncol       ) :: brup
    real(r8), dimension(ncol       ) :: phim
    real(r8), dimension(ncol       ) :: phih
    real(r8), dimension(ncol       ) :: prpbl
    real(r8), dimension(ncol       ) :: wspd1
    real(r8), dimension(ncol       ) :: thermalli
    real(r8), dimension(ncol,nlev  ) :: f1
    real(r8), dimension(ncol,nlev  ) :: f2
    real(r8), dimension(ncol,nlev  ) :: rhox2
    real(r8), dimension(ncol,nlev  ) :: r1
    real(r8), dimension(ncol,nlev  ) :: r2
    real(r8), dimension(ncol,nlev  ) :: ad
    real(r8), dimension(ncol,nlev  ) :: au
    real(r8), dimension(ncol,nlev  ) :: cu
    real(r8), dimension(ncol,nlev  ) :: al
    real(r8), dimension(ncol       ) :: brcr
    real(r8), dimension(ncol       ) :: sflux
    real(r8), dimension(ncol       ) :: zol1
    real(r8), dimension(ncol       ) :: brcr_sbro
    real(r8), dimension(ncol,nlev  ) :: zfac3
    real(r8), dimension(ncol,nlev  ) :: entfac
    real(r8), dimension(ncol       ) :: ustar3
    real(r8), dimension(ncol       ) :: wstar3
    real(r8), dimension(ncol       ) :: wstar3_2
    real(r8), dimension(ncol       ) :: hgamu
    real(r8), dimension(ncol       ) :: hgamv
    real(r8), dimension(ncol       ) :: we
    real(r8), dimension(ncol       ) :: hfxpbl
    real(r8), dimension(ncol       ) :: qfxpbl
    real(r8), dimension(ncol       ) :: ufxpbl
    real(r8), dimension(ncol       ) :: vfxpbl
    real(r8), dimension(ncol,nlev  ) :: tke
    real(r8), dimension(ncol,nlev  ) :: fric
    real(r8), dimension(ncol       ) :: pblh
    real(r8), dimension(ncol,nlev,3) :: r3, f3

    logical , dimension(ncol       ) :: pblflg
    logical , dimension(ncol       ) :: sfcflg
    logical , dimension(ncol       ) :: stable
    logical , dimension(ncol       ) :: cloudflg

    real(r8) dt2, rdt
    real(r8) ws1, ws2, zfac
    real(r8) hol1, gamfac, vpert, ss, ri, tmp, dp_lev, rvls, ent_eff, rcldb
    real(r8) radsum, radflux
    real(r8) prnum, prnum0, prnumfac, prfac, prfac2
    real(r8) vconvlim, vconv
    real(r8) qmean, tmean
    real(r8) alph, chi, zk, dk, sri
    real(r8) brint, dtodsd, dtodsu
    real(r8) dsdzt, dsdzq, dsdz2, rlamdz
    real(r8) utend, vtend, ttend, qtend
    real(r8) cont, conq, conw, conwrc, conpr
    real(r8) wm2, wm3, bfxpbl, dthvx, bfx0
    real(r8) hfx0, qfx0 ! 没啥用
    real(r8) dux, dvx, dsdzu, dsdzv
    real(r8) shear, buoy
    real(r8) ep1, ep2
    real(r8) templ, temps
    logical definebrup
    integer i, k, kk

    cont   = cp / g
    conq   = xlv / g
    conw   = 1.0_r8 / g
    conwrc = conw * sqrt(rcl)
    conpr  = bfac * karman * sfcfrac
    ep1    = rv_o_rd - 1
    ep2    = rd_o_rv

    do k = 1, nlev
      do i = 1, ncol
        th  (i,k) = t(i,k) / pk(i,k)
        thli(i,k) = (t(i,k) - xlv * qc(i,k) / cp - 2.834e6_r8 * qi(i,k) / cp) / pk(i,k)
      end do
    end do

    do k = 1, nlev
      do i = 1, ncol
        thv(i,k) = th(i,k) * (1 + ep1 * qv(i,k)) ! 虚位温？
      end do
    end do

    do i = 1, ncol
      rhos(i) = ps(i) / (rd * t(i,1) * (1 + ep1 * qv(i,1))) ! 地表空气密度
      govrth(i) = g / th(i,1)
    end do

    ! Compute the height of full- and half-sigma levels above ground
    ! level, and the layer thicknesses.
    do i = 1, ncol
      z_lev(i,1) = 0
    end do

    do k = 1, nlev
      do i = 1, ncol
        z_lev(i,k+1) = dz(i,k) + z_lev(i,k) ! 半层距地面高度
        rhox2(i,k) = p(i,k) / (rd * t(i,k) * (1 + ep1 * qv(i,k))) ! 整层空气密度
      end do
    end do

    do k = 1, nlev
      do i = 1, ncol
        z (i,k) = 0.5_r8 * (z_lev(i,k) + z_lev(i,k+1)) ! 整层距地面高度
        dp(i,k) = p_lev(i,k) - p_lev(i,k+1) ! 气压厚度
      end do
    end do

    do i = 1, ncol
      dz_lev(i,1) = z(i,1) ! dzi是整层厚度，dza是半层厚度
    end do

    do k = 2, nlev
      do i = 1, ncol
        dz_lev(i,k) = z(i,k) - z(i,k-1)
      end do
    end do

    ! Initialize vertical tendencies and
    dudt_pbl  = 0
    dvdt_pbl  = 0
    dtdt_pbl  = 0
    dqvdt_pbl = 0
    dqcdt_pbl = 0
    dqidt_pbl = 0

    do i = 1, ncol
      wspd1(i) = sqrt((u(i,1) - uos(i))**2 + (v(i,1) - vos(i))**2) + 1.0e-9_r8
    end do

    ! Compute vertical diffusion

    dt2 = 2 * dt
    rdt = 1.0_r8 / dt2

    do i = 1, ncol
      hfxpbl  (i) = 0
      qfxpbl  (i) = 0
      ufxpbl  (i) = 0
      vfxpbl  (i) = 0
      hgamu   (i) = 0
      hgamv   (i) = 0
      delta   (i) = 0
      wstar3_2(i) = 0
    end do

    do i = 1, ncol
      hgamt    (i) = 0
      hgamq    (i) = 0
      kpbl     (i) = 1
      hpbl     (i) = z_lev(i,1)
      thermal  (i) = thv  (i,1)
      thermalli(i) = thli (i,1)
      pblflg   (i) = .true.
      sfcflg   (i) = .true.
      sflux    (i) = hfx(i) / rhos(i) / cp + qfx(i) / rhos(i) * ep1 * th(i,1)
      if (br(i) > 0) sfcflg(i) = .false.
    end do

    ! Compute the first guess of pbl height
    do i = 1, ncol
      stable(i) = .false.
      brup  (i) = br(i)
      brcr  (i) = brcr_ub
    end do

    do k = 2, nlev
      do i = 1, ncol
        if (.not. stable(i)) then
          brdn  (i) = brup(i)
          brup  (i) = (thv(i,k) - thermal(i)) * (g * z(i,k) / thv(i,1)) / max(u(i,k)**2 + v(i,k)**2, 1.0_r8)
          kpbl  (i) = k
          stable(i) = brup(i) > brcr(i)
        end if
      end do
    end do

    do i = 1, ncol
      k = kpbl(i)
      if (brdn(i) >= brcr(i)) then
        brint = 0
      else if (brup(i) <= brcr(i)) then
        brint = 1
      else
        brint = (brcr(i) - brdn(i)) / (brup(i) - brdn(i))
      end if
      hpbl(i) = z(i,k-1) + brint * (z(i,k) - z(i,k-1))
      if (hpbl(i) < z_lev(i,2)) kpbl(i) = 1
      if (kpbl(i) <= 1) pblflg(i) = .false.
    end do

    do i = 1, ncol
      zol1(i) = max(br(i) * psim(i)**2 / psih(i), ri_min)
      if (sfcflg(i)) then
        zol1(i) = min(zol1(i), -zf_min)
      else
        zol1(i) = max(zol1(i),  zf_min)
      end if
      hol1 = zol1(i) * hpbl(i) / z(i,1) * sfcfrac
      if (sfcflg(i)) then
        phim  (i) = (1 - aphi16 * hol1)**(-0.25_r8)
        phih  (i) = (1 - aphi16 * hol1)**(-0.5_r8)
        bfx0      = max(sflux(i), 0.0_r8)
        hfx0      = max(hfx(i) / rhos(i) / cp, 0.0_r8)
        qfx0      = max(ep1 * th(i,1) * qfx(i) / rhos(i), 0.0_r8)
        wstar3(i) = (govrth(i) * bfx0 * hpbl(i))
        wstar (i) = wstar3(i)**h1
      else
        phim  (i) = 1 + aphi5 * hol1
        phih  (i) = phim(i)
        wstar (i) = 0
        wstar3(i) = 0
      end if
      ustar3(i) = ustar(i)**3
      ws0(i) = (ustar3(i) + phifac * karman * wstar3(i) * 0.5_r8)**h1
      ws0(i) = min(ws0(i), ustar(i) * aphi16)
      ws0(i) = max(ws0(i), ustar(i) / aphi5)
    end do

    ! Compute the surface variables for pbl height estimation under unstable conditions
    do i = 1, ncol
      if (sfcflg(i) .and. sflux(i) > 0) then
        gamfac       = bfac / rhos(i) / ws0(i)
        hgamt    (i) = min(gamfac * hfx(i) / cp, gamcrt)
        hgamq    (i) = min(gamfac * qfx(i), gamcrq)
        vpert        = (hgamt(i) + ep1 * th(i,1) * hgamq(i)) / bfac * afac
        thermal  (i) = thermal  (i) + max(vpert, 0.0_r8) * min(z(i,1) / (sfcfrac * hpbl(i)), 1.0_r8)
        thermalli(i) = thermalli(i) + max(vpert, 0.0_r8) * min(z(i,1) / (sfcfrac * hpbl(i)), 1.0_r8)
        hgamt    (i) = max(hgamt(i), 0.0_r8)
        hgamq    (i) = max(hgamq(i), 0.0_r8)
        brint        = -15.9_r8 * ustar(i) * ustar(i) / wspd(i) * wstar3(i) / (ws0(i)**4)
        hgamu    (i) = brint * u(i,1)
        hgamv    (i) = brint * v(i,1)
      else
        pblflg   (i) = .false.
      end if
    end do

    ! Enhance the PBL height by considering the thermal
    do i = 1, ncol
      if (pblflg(i)) then
        kpbl(i) = 1
        hpbl(i) = z_lev(i,1)
      end if
    end do

    do i = 1, ncol
      if (pblflg(i)) then
        stable(i) = .false.
        brup  (i) = br(i)
        brcr  (i) = brcr_ub
      end if
    end do

    do k = 2, nlev
      do i = 1, ncol
        if (.not.stable(i) .and. pblflg(i)) then
          brdn  (i) = brup(i)
          brup  (i) = (thv(i,k) - thermal(i)) * (g * z(i,k) / thv(i,1)) / max(u(i,k)**2 + v(i,k)**2, 1.0_r8)
          kpbl  (i) = k
          stable(i) = brup(i) > brcr(i)
        end if
      end do
    end do

    ! Enhance PBL by theta-li
    if (ysu_topdown_pblmix == 1)then
      do i = 1, ncol
        definebrup = .false.
        do k = kpbl(i), nlev - 1
          tmp = (thli(i,k) - thermalli(i)) * (g * z(i,k) / thli(i,1)) / max(u(i,k)**2 + v(i,k)**2, 1.0_r8)
          stable(i) = tmp >= brcr(i)
          if (definebrup) then
            kpbl(i)    = k
            brup(i)    = tmp
            definebrup = .false.
          end if
          if (.not.stable(i)) then ! Overwrite brup brdn values
            brdn(i)    = tmp
            definebrup = .true.
            pblflg(i)  = .true.
          end if
        end do
      end do
    end if

    do i = 1, ncol
      if (pblflg(i)) then
        k = kpbl(i)
        if (brdn(i) >= brcr(i)) then
          brint = 0.
        else if (brup(i) <= brcr(i)) then
          brint = 1
        else
          brint = (brcr(i) - brdn(i)) / (brup(i) - brdn(i))
        end if
        hpbl(i) = z(i,k-1) + brint * (z(i,k) - z(i,k-1))
        if (hpbl(i) < z_lev(i,2)) kpbl(i) = 1
        if (kpbl(i) <= 1) pblflg(i) = .false.
      end if
    end do

    ! Stable boundary layer
    do i = 1, ncol
      if ((.not.sfcflg(i)) .and. hpbl(i) < z_lev(i,2)) then
        brup  (i) = br(i)
        stable(i) = .false.
      else
        stable(i) = .true.
      end if
    end do

    do i = 1, ncol
      if ((.not.stable(i)) .and. ((xland(i)-1.5) >= 0)) then
        brcr_sbro(i) = min(0.16_r8 * (1.0e-7_r8 * sqrt(u10(i)**2 + v10(i)**2) / (cori * znt(i)))**(-0.18_r8), 0.3_r8)
      end if
    end do

    do i = 1, ncol
      if (.not.stable(i)) then
        if (xland(i) - 1.5 >= 0) then
          brcr(i) = brcr_sbro(i)
        else
          brcr(i) = brcr_sb
        end if
      end if
    end do

    do k = 2, nlev
      do i = 1, ncol
        if (.not.stable(i)) then
          brdn(i) = brup(i)
          brup(i) = (thv(i,k) - thermal(i)) * (g * z(i,k) / thv(i,1)) / max(u(i,k)**2 + v(i,k)**2, 1.0_r8)
          kpbl(i) = k
          stable(i) = brup(i).gt.brcr(i)
        end if
      end do
    end do

    do i = 1, ncol
      if (.not. sfcflg(i) .and. hpbl(i) < z_lev(i,2)) then
        k = kpbl(i)
        if (brdn(i) >= brcr(i)) then
          brint = 0
        else if (brup(i) <= brcr(i)) then
          brint = 1
        else
          brint = (brcr(i) - brdn(i)) / (brup(i) - brdn(i))
        end if
        hpbl(i) = z(i,k-1) + brint * (z(i,k) - z(i,k-1))
        if (hpbl(i) < z_lev(i,2)) kpbl(i) = 1
        if (kpbl(i) <= 1) pblflg(i) = .false.
      end if
    end do

    ! Estimate the entrainment parameters.
    do i = 1, ncol
      cloudflg(i) = .false. 
      if (pblflg(i)) then
        k = kpbl(i) - 1
        wm3       = wstar3(i) + 5 * ustar3(i)
        wm2       = wm3**h2
        bfxpbl    = -0.15_r8 * thv(i,1) / g * wm3 / hpbl(i)
        dthvx     = max(thv(i,k+1) - thv(i,k), tmin)
        we    (i) = max(bfxpbl / dthvx, -sqrt(wm2))
        if (qc(i,k) + qi(i,k) > 0.01e-3_r8 .and. ysu_topdown_pblmix == 1) then
          if (kpbl(i) >= 2) then
            cloudflg(i) = .true. 
            templ = thli(i,k) * (p_lev(i,k+1) / 100000)**rovcp
            ! rvls is ws at full level
            rvls  = 100 * 6.112_r8 * exp(17.67_r8 * (templ - 273.16_r8) / (templ - 29.65_r8)) * (ep2 / p_lev(i,k+1))
            temps = templ + ((qv(i,k) + qc(i,k)) - rvls) / (cp / xlv + ep2 * xlv * rvls / (rd * templ**2))
            rvls  = 100 * 6.112_r8 * exp(17.67_r8 * (temps - 273.15_r8) / (temps - 29.65_r8)) * (ep2 / p_lev(i,k+1))
            rcldb = max(qv(i,k) + qc(i,k) - rvls, 0.0_r8)
            ! entrainment efficiency
            dthvx = (thli(i,k+2) + th(i,k+2) * ep1 * (qv(i,k+2) + qc(i,k+2))) &
                  - (thli(i,k  ) + th(i,k  ) * ep1 * (qv(i,k  ) + qc(i,k  )))
            dthvx = max(dthvx, 0.1_r8)
            ent_eff  = 0.2_r8 * 8 * xlv / cp * rcldb / (pk(i,k) * dthvx) + 0.2_r8
            radsum   = 0
            do kk = 1, kpbl(i) - 1
              ! Converts theta/s to temp/s
              radflux = dptdt_rad(i,kk) * pk(i,kk)
              ! Converts temp/s to W/m^2
              radflux = radflux * cp / g * (p_lev(i,kk) - p_lev(i,kk+1))
              if (radflux < 0.0) radsum = abs(radflux) + radsum
            end do
            radsum = max(radsum, 0.0_r8)

            ! Recompute entrainment from sfc thermals
            bfx0        = max(max(sflux(i), 0.0_r8) - radsum / rhox2(i,k) / cp, 0.0_r8)
            bfx0        = max(sflux(i), 0.0_r8) ! 啥意思？上一行不作数？
            wm3         = (govrth(i) * bfx0 * hpbl(i)) + 5 * ustar3(i)
            wm2         = wm3**h2
            bfxpbl      = -0.15_r8 * thv(i,1) / g * wm3 / hpbl(i)
            dthvx       = max(thv(i,k+1) - thv(i,k), tmin)
            we    (i)   = max(bfxpbl / dthvx, -sqrt(wm2))

            ! Entrainment from PBL top thermals
            bfx0        = max(radsum / rhox2(i,k) / cp - max(sflux(i), 0.0_r8), 0.0_r8)
            bfx0        = max(radsum / rhox2(i,k) / cp, 0.0_r8) ! 啥意思？上一行不作数？
            wm3         = (g / thv(i,k) * bfx0 * hpbl(i)) ! this is wstar3(i)
            wm2         = wm2 + wm3**h2
            bfxpbl      = -ent_eff * bfx0
            dthvx       = max(thv(i,k+1)-thv(i,k),0.1)
            we    (i)   = we(i) + max(bfxpbl / dthvx, -sqrt(wm3**h2))

            ! wstar3_2
            bfx0        = max(radsum / rhox2(i,k) / cp, 0.0_r8)
            wstar3_2(i) = (g / thv(i,k) * bfx0 * hpbl(i))
            ! Recompute hgamt 
            ws0   (i)   = (ustar3(i) + phifac * karman * (wstar3(i) + wstar3_2(i)) * 0.5_r8)**h1
            ws0   (i)   = min(ws0(i), ustar(i) * aphi16)
            ws0   (i)   = max(ws0(i), ustar(i) / aphi5)
            gamfac      = bfac / rhos(i) / ws0(i)
            hgamt (i)   = min(gamfac * hfx(i) / cp, gamcrt)
            hgamq (i)   = min(gamfac * qfx(i), gamcrq)
            gamfac      = bfac / rhox2(i,k) / ws0(i)
            hgamt (i)   = max(hgamt(i), 0.0_r8) + max(min(gamfac * radsum / cp, gamcrt), 0.0_r8)
            brint       = -15.9_r8 * ustar(i) * ustar(i) / wspd(i) * (wstar3(i) + wstar3_2(i)) / (ws0(i)**4)
            hgamu (i)   = brint * u(i,1)
            hgamv (i)   = brint * v(i,1)
          end if
        end if
        prpbl (i) = 1
        hfxpbl(i) = we(i) * max(th(i,k+1) - th(i,k), tmin)
        qfxpbl(i) = we(i) * min(qv(i,k+1) - qv(i,k), 0.0_r8)

        dux = u(i,k+1) - u(i,k)
        dvx = v(i,k+1) - v(i,k)
        if (dux > tmin) then
          ufxpbl(i) = max(prpbl(i) * we(i) * dux, -ustar(i) * ustar(i))
        else if (dux < -tmin) then
          ufxpbl(i) = min(prpbl(i) * we(i) * dux,  ustar(i) * ustar(i))
        else
          ufxpbl(i) = 0
        end if
        if (dvx > tmin) then
          vfxpbl(i) = max(prpbl(i) * we(i) * dvx, -ustar(i) * ustar(i))
        else if (dvx < -tmin) then
          vfxpbl(i) = min(prpbl(i) * we(i) * dvx,  ustar(i) * ustar(i))
        else
          vfxpbl(i) = 0
        end if
        delta(i) = min(d1 * hpbl(i) + d2 * wm2 / govrth(i) * d3 * hpbl(i), 100.0_r8)
      end if
    end do

    do k = 1, nlev
      do i = 1, ncol
        if (pblflg(i) .and. k >= kpbl(i))then
          entfac(i,k) = ((z_lev(i,k+1) - hpbl(i)) / delta(i))**2
        else
          entfac(i,k) = 1.0e30
        end if
      end do
    end do

    ! Compute diffusion coefficients below PBL.
    do k = 1, nlev
      do i = 1, ncol
        if (k < kpbl(i)) then
          zfac = min(max((1 - (z_lev(i,k+1) - z(i,1)) / (hpbl(i) - z(i,1))), zf_min), 1.0_r8)
          zfac3(i,k) = (1 - zfac)**3
          ws1 = (ustar3(i) + phifac * karman * wstar3(i) * (1 - zfac))**h1
          ws2 = (phifac * karman * wstar3_2(i) * zfac)**h1
          if (sfcflg(i)) then
            prfac    = conpr
            prfac2   = 15.9_r8 * (wstar3(i) + wstar3_2(i)) / ustar3(i) / (1 + 4 * karman * (wstar3(i) + wstar3_2(i)) / ustar3(i))
            prnumfac = -3 * (max(z_lev(i,k+1) - sfcfrac * hpbl(i), 0.0_r8))**2 / hpbl(i)**2
          else
            prfac    = 0
            prfac2   = 0
            prnumfac = 0
            ws1 = max(ustar(i) / (1 + aphi5 * zol1(i) * z_lev(i,k+1) / z(i,1)), 0.001_r8)
          end if
          prnum0 = max(min((phih(i) / phim(i) + prfac), prmax), prmin)
          kzm(i,k) = ws1 * karman *           z_lev(i,k+1)  *      zfac **pfac + &
                     ws2 * karman *(hpbl(i) - z_lev(i,k+1)) * (1 - zfac)**pfac
          ! Do not include kzm at kpbl-1 since it changes entrainment
          if (k == kpbl(i) - 1 .and. cloudflg(i) .and. we(i) < 0) then
            kzm(i,k) = 0
          end if
          prnum    = 1 + (prnum0 - 1) * exp(prnumfac)
          kzq(i,k) = kzm(i,k) / prnum * zfac**(pfac_q - pfac)
          prnum0   = prnum0 / (1 + prfac2 * karman * sfcfrac)
          prnum    = 1 + (prnum0 - 1) * exp(prnumfac)
          kzh(i,k) = kzm(i,k) / prnum
          kzm(i,k) = min(kzm(i,k) + kzm_min, kz_max)
          kzh(i,k) = min(kzh(i,k) + kzh_min, kz_max)
          kzq(i,k) = min(kzq(i,k) + kzh_min, kz_max)
        end if
      end do
    end do

    ! Compute diffusion coefficients over PBL (free atmosphere)
    do k = 1, nlev - 1
      do i = 1, ncol
        if (k >= kpbl(i)) then
          ss = ((u(i,k+1) - u(i,k))**2 + (v(i,k+1) - v(i,k))**2) / dz_lev(i,k+1)**2 + 1.0e-9_r8
          ri = g / (0.5_r8 * (thv(i,k+1) + thv(i,k))) * (thv(i,k+1) - thv(i,k)) / (ss * dz_lev(i,k+1))
          if (imvdif == 1 .and. pbl_vdiff_qc .and. pbl_vdiff_qi)then
            if (qc(i,k) + qi(i,k) > 0.01e-3_r8 .and. (qc(i,k+1) + qi(i,k+1)) > 0.01e-3) then
              ! In cloud
              qmean = 0.5_r8 * (qv(i,k) + qv(i,k+1))
              tmean = 0.5_r8 * (t (i,k) + t (i,k+1))
              alph  = xlv * qmean / rd / tmean
              chi   = xlv * xlv * qmean / cp / rv / tmean / tmean
              ri    = (1 + alph) * (ri - g * g / ss / tmean / cp * ((chi - alph) / (1 + chi)))
            end if
          end if
          zk = karman * z_lev(i,k+1)
          rlamdz = min(max(0.1_r8 * dz_lev(i,k+1), rlam), 300.0_r8)
          rlamdz = min(dz_lev(i,k+1), rlamdz)
          dk     = (zk * rlamdz / (rlamdz + zk))**2 * sqrt(ss)
          if (ri < 0) then
            ! Unstable regime
            ri        = max(ri, ri_min)
            sri       = sqrt(-ri)
            kzm (i,k) = dk * (1 + 8 * (-ri) / (1 + 1.746_r8 * sri))
            kzh(i,k) = dk * (1 + 8 * (-ri) / (1 + 1.286_r8 * sri))
          else
            ! Stable regime
            kzh(i,k) = dk / (1 + 5 * ri)**2
            prnum     = min(1.0 + 2.1 * ri, prmax)
            kzm(i,k) = kzh(i,k) * prnum
          end if
          kzm (i,k) = min(kzm(i,k) + kzm_min, kz_max)
          kzh (i,k) = min(kzh(i,k) + kzh_min, kz_max)
          kzml(i,k) = kzm(i,k)
          kzhl(i,k) = kzh(i,k)
        end if
      end do
    end do

    ! Compute tridiagonal matrix elements for heat
    do k = 1, nlev
      do i = 1, ncol
        au(i,k) = 0
        al(i,k) = 0
        ad(i,k) = 0
        f1(i,k) = 0
      end do
    end do

    do i = 1, ncol
      ad(i,1) = 1
      f1(i,1) = th(i,1) - 300 + hfx(i) / cont / dp(i,1) * dt2
    end do

    do k = 1, nlev - 1
      do i = 1, ncol
        dtodsd = dt2 / dp(i,k  )
        dtodsu = dt2 / dp(i,k+1)
        dp_lev = p(i,k) - p(i,k+1)
        tmp    = dp_lev * kzh(i,k) / dz_lev(i,k+1)
        if (pblflg(i) .and. k < kpbl(i)) then
          dsdzt     = tmp * (-hgamt(i) / hpbl(i) - hfxpbl(i) * zfac3(i,k) / kzh(i,k))
          f1(i,k)   = f1(i,k) + dtodsd * dsdzt
          f1(i,k+1) = th(i,k+1) - 300 - dtodsu * dsdzt
        else if (pblflg(i) .and. k >= kpbl(i) .and. entfac(i,k) < 4.6_r8) then
          kzh(i,k)  = -we(i) * dz_lev(i,kpbl(i)) * exp(-entfac(i,k))
          kzh(i,k)  = sqrt(kzh(i,k) * kzhl(i,k))
          kzh(i,k)  = min(max(kzh(i,k), kzh_min), kz_max)
          f1(i,k+1) = th(i,k+1) - 300
        else
          f1(i,k+1) = th(i,k+1) - 300
        end if
        tmp       = dp_lev * kzh(i,k) / dz_lev(i,k+1)
        dsdz2     = tmp / dz_lev(i,k+1)
        au(i,k  ) = -dtodsd * dsdz2
        al(i,k  ) = -dtodsu * dsdz2
        ad(i,k  ) = ad(i,k) - au(i,k)
        ad(i,k+1) = 1 - al(i,k)
      end do
    end do

    ! Copies here to avoid duplicate input args for tridin
    do k = 1, nlev
      do i = 1, ncol
       cu(i,k) = au(i,k)
       r1(i,k) = f1(i,k)
      end do
    end do

    call tridin_ysu(ncol, nlev, 1, al, ad, cu, r1, au, f1)

    ! Recover tendencies of heat
    do k = nlev, 1, -1
      do i = 1, ncol
        ttend = (f1(i,k) - th(i,k) + 300) * rdt * pk(i,k)
        dtdt_pbl(i,k) = dtdt_pbl(i,k) + ttend
      end do
    end do

    ! Compute tridiagonal matrix elements for moisture, clouds, and gases
    do k = 1, nlev
      do i = 1, ncol
        au(i,k) = 0
        al(i,k) = 0
        ad(i,k) = 0
        f3(i,k,:) = 0
      end do
    end do

    do i = 1, ncol
      ad(i,1) = 1
      f3(i,1,1) = qv(i,1) + qfx(i) * g / dp(i,1) * dt2
    end do

    if (pbl_vdiff_qc) then
      do i = 1, ncol
        f3(i,1,2) = qc(i,1)
      end do
    end if
    if (pbl_vdiff_qi) then
      do i = 1, ncol
        f3(i,1,3) = qi(i,1)
      end do
    end if

    do k = 1, nlev - 1
      do i = 1, ncol
        if (k >= kpbl(i)) then
          kzq(i,k) = kzh(i,k)
        end if
      end do
    end do

    do k = 1, nlev - 1
      do i = 1, ncol
        dtodsd = dt2 / dp(i,k  )
        dtodsu = dt2 / dp(i,k+1)
        dp_lev = p(i,k) - p(i,k+1)
        if (pblflg(i) .and. k < kpbl(i)) then
          dsdzq = dp_lev * kzq(i,k) / dz_lev(i,k+1) * (-qfxpbl(i) * zfac3(i,k) / kzq(i,k))
          f3(i,k,1) = f3(i,k,1) + dtodsd * dsdzq
          f3(i,k+1,1) = qv(i,k+1) - dtodsu * dsdzq
        else if (pblflg(i) .and. k >= kpbl(i) .and. entfac(i,k) < 4.6) then
          kzq(i,k) = -we(i) * dz_lev(i,kpbl(i)) * exp(-entfac(i,k))
          kzq(i,k) = sqrt(kzq(i,k) * kzhl(i,k))
          kzq(i,k) = min(max(kzq(i,k), kzh_min), kz_max)
          f3(i,k+1,1) = qv(i,k+1)
        else
          f3(i,k+1,1) = qv(i,k+1)
        end if
        dsdz2     = dp_lev * kzq(i,k) / dz_lev(i,k+1)**2
        au(i,k)   = -dtodsd * dsdz2
        al(i,k)   = -dtodsu * dsdz2
        ad(i,k)   = ad(i,k) - au(i,k)
        ad(i,k+1) = 1 - al(i,k)
      end do
    end do

    if (pbl_vdiff_qc) then
      do k = 1, nlev - 1
        do i = 1, ncol
          f3(i,k+1,2) = qc(i,k+1)
        end do
      end do
    end if
    if (pbl_vdiff_qi) then
      do k = 1, nlev - 1
        do i = 1, ncol
          f3(i,k+1,3) = qi(i,k+1)
        end do
      end do
    end if

    ! Copies here to avoid duplicate input args for tridin.
    do k = 1, nlev
      do i = 1, ncol
        cu(i,k) = au(i,k)
        r3(i,k,:) = f3(i,k,:)
      end do
    end do

    ! Solve tridiagonal problem for moisture, clouds, and gases.
    call tridin_ysu(ncol, nlev, 3, al, ad, cu, r3, au, f3)

    ! Recover tendencies of heat and moisture
    do k = nlev, 1, -1
      do i = 1, ncol
       qtend = (f3(i,k,1) - qv(i,k)) * rdt
       dqvdt_pbl(i,k) = dqvdt_pbl(i,k) + qtend
      end do
    end do

    if (pbl_vdiff_qc) then
      do k = nlev, 1, -1
        do i = 1, ncol
          qtend = (f3(i,k,2) - qc(i,k)) * rdt
          dqcdt_pbl(i,k) = dqcdt_pbl(i,k) + qtend
        end do
      end do
    end if
    if (pbl_vdiff_qi) then
      do k = nlev, 1, -1
        do i = 1, ncol
          qtend = (f3(i,k,3) - qi(i,k)) * rdt
          dqidt_pbl(i,k) = dqidt_pbl(i,k) + qtend
        end do
      end do
    end if

    ! Compute tridiagonal matrix elements for momentum.
    do k = 1, nlev
      do i = 1, ncol
        au(i,k) = 0
        al(i,k) = 0
        ad(i,k) = 0
        f1(i,k) = 0
        f2(i,k) = 0
      end do
    end do

    ! paj: ctopo=1 if topo_wind=0 (default)
    do i = 1, ncol
      do k = 1, nlev - 1
        shear = kzm(i,k) * ((-hgamu(i) / hpbl(i) + (u(i,k+1) - u(i,k)) / dz_lev(i,k+1)) * (u(i,k+1) - u(i,k)) / dz_lev(i,k+1) &
                         +  (-hgamv(i) / hpbl(i) + (v(i,k+1) - v(i,k)) / dz_lev(i,k+1)) * (v(i,k+1) - v(i,k)) / dz_lev(i,k+1))
        buoy  = kzh(i,k) * g / th(i,k) * (-hgamt(i) / hpbl(i) + (th(i,k+1) - th(i,k)) / dz_lev(i,k+1))

        zk = karman * z_lev(i,k+1)
        if (k >= kpbl(i)) then
          ! Over PBL
          rlamdz = min(dz_lev(i,k+1), min(max(0.1_r8 * dz_lev(i,k+1), rlam), 300.0_r8))
        else
          ! In PBL
          rlamdz = 150
        end if
        tke(i,k) = 16.6_r8 * zk * rlamdz / (rlamdz + zk) * (shear - buoy)
        if (tke(i,k) <= 0) then
          tke(i,k) = 0
        else
          tke(i,k) = tke(i,k)**0.66_r8
        end if
      end do

      ! Hybrid PBLH of MYNN
      call get_pblh(nlev, thv(i,:), tke(i,:), z_lev(i,:), dz(i,:), xland(i), pblh(i))

      ! End of paj TKE
      
      ! Compute vconv
      ! Use Beljaars over land
      if (xland(i) < 1.5_r8) then
        vconv = (g / thv(i,1) * pblh(i) * max(sflux(i), 0.0_r8))**0.33_r8
      else
        ! For water there is no topo effect so vconv not needed
        vconv = 0
      end if
      ! raquel
      ! ctopo stability correction
      fric(i,1) = ustar(i)**2 / wspd1(i) * rhos(i) * g / dp(i,1) * dt2 * (wspd1(i) / wspd(i))**2
      if (present(ctopo)) then
        vconvlim = min(0.9_r8 * vconv + 1.5_r8 * (max((pblh(i) - 500) / 1000.0_r8, 0.0_r8)), 1.0_r8)
        ad(i,1) = 1 + fric(i,1) * vconvlim + ctopo(i) * fric(i,1) * (1 - vconvlim)
      else
       ad(i,1) = 1 + fric(i,1)
      end if
      f1(i,1) = u(i,1) + uos(i) * ustar(i)**2 * rhos(i) * g / dp(i,1) * dt2 / wspd1(i) * (wspd1(i) / wspd(i))**2
      f2(i,1) = v(i,1) + vos(i) * ustar(i)**2 * rhos(i) * g / dp(i,1) * dt2 / wspd1(i) * (wspd1(i) / wspd(i))**2
    end do

    do k = 1, nlev - 1
      do i = 1, ncol
        dtodsd = dt2 / dp(i,k  )
        dtodsu = dt2 / dp(i,k+1)
        dp_lev = p(i,k) - p(i,k+1)
        tmp    = dp_lev * kzm(i,k) / dz_lev(i,k+1)
        if (pblflg(i) .and. k < kpbl(i)) then
          dsdzu     = tmp * (-hgamu(i) / hpbl(i) - ufxpbl(i) * zfac3(i,k) / kzm(i,k))
          dsdzv     = tmp * (-hgamv(i) / hpbl(i) - vfxpbl(i) * zfac3(i,k) / kzm(i,k))
          f1(i,k  ) = f1(i,k  ) + dtodsd * dsdzu
          f1(i,k+1) = u(i,k+1) - dtodsu * dsdzu
          f2(i,k  ) = f2(i,k  ) + dtodsd * dsdzv
          f2(i,k+1) = v(i,k+1) - dtodsu * dsdzv
        else if (pblflg(i) .and. k >= kpbl(i) .and. entfac(i,k) < 4.6_r8) then
          kzm(i,k) = prpbl(i) * kzh(i,k)
          kzm(i,k) = sqrt(kzm(i,k) * kzml(i,k))
          kzm(i,k) = min(max(kzm(i,k), kzm_min), kz_max)
          f1(i,k+1) = u(i,k+1)
          f2(i,k+1) = v(i,k+1)
        else
          f1(i,k+1) = u(i,k+1)
          f2(i,k+1) = v(i,k+1)
        end if
        tmp       = dp_lev * kzm(i,k) / dz_lev(i,k+1)
        dsdz2     = tmp / dz_lev(i,k+1)
        au(i,k  ) = -dtodsd * dsdz2
        al(i,k  ) = -dtodsu * dsdz2
        ad(i,k  ) = ad(i,k) - au(i,k)
        ad(i,k+1) = 1 - al(i,k)
      end do
    end do

    ! Copies here to avoid duplicate input args for tridin
    do k = 1, nlev
      do i = 1, ncol
       cu(i,k) = au(i,k)
       r1(i,k) = f1(i,k)
       r2(i,k) = f2(i,k)
      end do
    end do

    ! Solve tridiagonal problem for momentum
    call tridi1n(ncol, nlev, al, ad, cu, r1, r2, au, f1, f2)

    ! Recover tendencies of momentum
    do k = nlev, 1, -1
      do i = 1, ncol
        utend = (f1(i,k) - u(i,k)) * rdt
        vtend = (f2(i,k) - v(i,k)) * rdt
        dudt_pbl(i,k) = dudt_pbl(i,k) + utend
        dvdt_pbl(i,k) = dvdt_pbl(i,k) + vtend
      end do
    end do

    ! paj: ctopo2=1 if topo_wind=0 (default)

    do i = 1, ncol
      if (present(ctopo) .and. present(ctopo2)) then ! mchen for NMM
        u10(i) = ctopo2(i) * u10(i) + (1 - ctopo2(i)) * u(i,1)
        v10(i) = ctopo2(i) * v10(i) + (1 - ctopo2(i)) * v(i,1)
      end if
    end do

    ! End of vertical diffusion

  end subroutine pbl_ysu_run

  subroutine tridi1n(ncol, nlev, cl, cm, cu, r1, r2, au, f1, f2)

    integer , intent(in   ) :: ncol
    integer , intent(in   ) :: nlev
    real(r8), intent(in   ), dimension(ncol,2:nlev+1) :: cl
    real(r8), intent(in   ), dimension(ncol,1:nlev  ) :: cm
    real(r8), intent(in   ), dimension(ncol,1:nlev  ) :: cu
    real(r8), intent(in   ), dimension(ncol,1:nlev  ) :: r1
    real(r8), intent(in   ), dimension(ncol,1:nlev  ) :: r2
    real(r8), intent(inout), dimension(ncol,1:nlev  ) :: au
    real(r8), intent(inout), dimension(ncol,1:nlev  ) :: f1
    real(r8), intent(inout), dimension(ncol,1:nlev  ) :: f2

    real(r8) fk
    integer i, k

    do i = 1, ncol
      fk = 1.0_r8 / cm(i,1)
      au(i,1) = fk * cu(i,1)
      f1(i,1) = fk * r1(i,1)
    end do

    do i = 1, ncol
      fk = 1.0_r8 / cm(i,1)
      f2(i,1) = fk * r2(i,1)
    end do

    do k = 2, nlev - 1
      do i = 1, ncol
        fk = 1.0_r8 / (cm(i,k) - cl(i,k) * au(i,k-1))
        au(i,k) = fk * cu(i,k)
        f1(i,k) = fk * (r1(i,k) - cl(i,k) * f1(i,k-1))
        f2(i,k) = fk * (r2(i,k) - cl(i,k) * f2(i,k-1))
      end do
    end do

    do i = 1, ncol
      fk = 1.0_r8 / (cm(i,nlev) - cl(i,nlev) * au(i,nlev-1))
      f1(i,nlev) = fk * (r1(i,nlev) - cl(i,nlev) * f1(i,nlev-1))
      f2(i,nlev) = fk * (r2(i,nlev) - cl(i,nlev) * f2(i,nlev-1))
    end do

    do k = nlev - 1, 1, -1
      do i = 1, ncol
        f1(i,k) = f1(i,k) - au(i,k) * f1(i,k+1)
        f2(i,k) = f2(i,k) - au(i,k) * f2(i,k+1)
      end do
    end do

  end subroutine tridi1n
   
  subroutine tridin_ysu(ncol, nlev, nt, cl, cm, cu, r2, au, f2)

    integer , intent(in   ) :: ncol
    integer , intent(in   ) :: nlev
    integer , intent(in   ) :: nt
    real(r8), intent(in   ), dimension(ncol,2:nlev+1 ) :: cl
    real(r8), intent(in   ), dimension(ncol,1:nlev   ) :: cm
    real(r8), intent(in   ), dimension(ncol,1:nlev   ) :: cu
    real(r8), intent(in   ), dimension(ncol,1:nlev,nt) :: r2
    real(r8), intent(inout), dimension(ncol,1:nlev   ) :: au
    real(r8), intent(inout), dimension(ncol,1:nlev,nt) :: f2

    real(r8) fk
    integer i, k, it

    do it = 1, nt
      do i = 1, ncol
        fk = 1.0_r8 / cm(i,1)
        au(i,1   ) = fk * cu(i,1   )
        f2(i,1,it) = fk * r2(i,1,it)
      end do
    end do

    do it = 1, nt
      do k = 2, nlev - 1
        do i = 1, ncol
          fk = 1.0_r8 / (cm(i,k) - cl(i,k) * au(i,k-1))
          au(i,k   ) = fk * cu(i,k)
          f2(i,k,it) = fk * (r2(i,k,it) - cl(i,k) * f2(i,k-1,it))
        end do
      end do
    end do

    do it = 1, nt
      do i = 1, ncol
        fk = 1.0_r8 / (cm(i,nlev) - cl(i,nlev) * au(i,nlev-1))
        f2(i,nlev,it) = fk * (r2(i,nlev,it) - cl(i,nlev) * f2(i,nlev-1,it))
      end do
    end do

    do it = 1, nt
      do k = nlev - 1, 1, -1
        do i = 1, ncol
          f2(i,k,it) = f2(i,k,it) - au(i,k) * f2(i,k+1,it)
        end do
      end do
    end do

  end subroutine tridin_ysu

  subroutine ysuinit(rublten,rvblten,rthblten,rqvblten,                       &
                      rqcblten,rqiblten,p_qi,p_first_scalar,                   &
                      restart, allowed_to_read,                                &
                      ids, ide, jds, jde, kds, kde,                            &
                      ims, ime, jms, jme, kms, kme,                            &
                      its, ite, jts, jte, kts, kte                 )

   logical , intent(in)          :: restart, allowed_to_read
   integer , intent(in)          ::  ids, ide, jds, jde, kds, kde,             &
                                     ims, ime, jms, jme, kms, kme,             &
                                     its, ite, jts, jte, kts, kte
   integer , intent(in)          ::  p_qi,p_first_scalar
   real(r8) , dimension( ims:ime , kms:kme , jms:jme ), intent(out) ::             &
                                                                      rublten, &
                                                                      rvblten, &
                                                                     rthblten, &
                                                                     rqvblten, &
                                                                     rqcblten, &
                                                                     rqiblten
   integer :: i, j, k, itf, jtf, ktf

   jtf = min0(jte,jde-1)
   ktf = min0(kte,kde-1)
   itf = min0(ite,ide-1)

   if(.not.restart)then
     do j = jts,jtf
       do k = kts,ktf
         do i = its,itf
            rublten(i,k,j) = 0.
            rvblten(i,k,j) = 0.
            rthblten(i,k,j) = 0.
            rqvblten(i,k,j) = 0.
            rqcblten(i,k,j) = 0.
         enddo
       enddo
     enddo
   endif

   if (p_qi .ge. p_first_scalar .and. .not.restart) then
     do j = jts,jtf
       do k = kts,ktf
         do i = its,itf
           rqiblten(i,k,j) = 0.
         enddo
       enddo
     enddo
   endif

  end subroutine ysuinit

  subroutine get_pblh(nlev, thetav, tke, zi, dz, landsea, pblh)

    !---------------------------------------------------------------
    !             NOTES ON THE PBLH FORMULATION
    !
    ! The 1.5-theta-increase method defines PBL heights as the level at
    ! which the potential temperature first exceeds the minimum potential
    ! temperature within the boundary layer by 1.5 K. When applied to
    ! observed temperatures, this method has been shown to produce PBL-
    ! height estimates that are unbiased relative to profiler-based
    ! estimates (Nielsen-Gammon et al. 2008). However, their study did not
    ! include LLJs. Banta and Pichugina (2008) show that a TKE-based
    ! threshold is a good estimate of the PBL height in LLJs. Therefore,
    ! a hybrid definition is implemented that uses both methods, weighting
    ! the TKE-method more during stable conditions (PBLH < 400 m).
    ! A variable TKE threshold (TKEeps) is used since no hard-wired
    ! value could be found to work best in all conditions.
    !---------------------------------------------------------------

      integer , intent(in ) :: nlev
      real(r8), intent(in ), dimension(nlev  ) :: thetav
      real(r8), intent(in ), dimension(nlev  ) :: tke
      real(r8), intent(in ), dimension(nlev+1) :: zi
      real(r8), intent(in ), dimension(nlev  ) :: dz
      real(r8), intent(in ) :: landsea
      real(r8), intent(out) :: pblh

      real(r8) pblh_tke, tke_lim, tke_km1, wgt, tke_max, tke_eps, thv_min
      real(r8) delt_thv   ! delta theta-v; dependent on land/sea point
      real(r8), parameter :: sbl_lim  = 200 !theta-v pbl lower limit of trust (m).
      real(r8), parameter :: sbl_damp = 400 !damping range for averaging with tke-based pblh (m).
      integer i, j, k, kthv, ktke

      ! Find max TKE and min thetav in the lowest 500m.
      k       = 2
      kthv    = 1
      ktke    = 1
      tke_max = 0
      thv_min = 9.0e9_r8

      do while (zi(k) <= 500) ! 垂直顺序是从底到顶？
        tke_lim = max(tke(k), 0.0_r8)
        if (tke_max < tke_lim) then
          tke_max = tke_lim
          ktke = k
        end if
        if (thv_min > thetav(k)) then
          thv_min = thetav(k)
          kthv = k
        end if
        k = k + 1
      end do
      tke_eps = min(max(tke_max / 40.0_r8, 0.025_r8), 0.25_r8)

      ! Find thetav-based PBLH (best for daytime).
      if (landsea - 1.5 >= 0) then
        ! Water
        delt_thv = 0.75_r8
      else
        ! Land
        delt_thv = 1.5_r8
      end if

      pblh = 0
      k = kthv + 1
      do while (pblh == 0)
        if (thetav(k) >= thv_min + delt_thv) then
          pblh = zi(k) - dz(k-1) * min((thetav(k) - (thv_min + delt_thv)) / max(thetav(k) - thetav(k-1), 1.0e-6_r8), 1.0_r8)
        end if
        k = k + 1
        if (k == nlev - 1) pblh = zi(2) ! Exit safeguard
      end do

      ! FOR STABLE BOUNDARY LAYERS, USE TKE METHOD TO COMPLEMENT THE
      ! THETAV-BASED DEFINITION (WHEN THE THETA-V BASED PBLH IS BELOW ~0.5 KM).
      ! THE TANH WEIGHTING FUNCTION WILL MAKE THE TKE-BASED DEFINITION NEGLIGIBLE
      ! WHEN THE THETA-V-BASED DEFINITION IS ABOVE ~1 KM.
      ! FIND TKE-BASED PBLH (BEST FOR NOCTURNAL/STABLE CONDITIONS).

      pblh_tke = 0
      k = ktke + 1
      do while (pblh_tke == 0.)
        ! TKE can be negative (if CKmod == 0)... make TKE non-negative.
        tke_lim = max(tke(k  ) / 2.0_r8, 0.0_r8)
        tke_km1 = max(tke(k-1) / 2.0_r8, 0.0_r8)
        if (tke_lim <= tke_eps) then
          pblh_tke = zi(k) - dz(k-1) * min((tke_eps - tke_lim) / max(tke_km1 - tke_lim, 1.0e-6_r8), 1.0_r8)
          ! IN CASE OF NEAR ZERO TKE, SET PBLH = LOWEST LEVEL.
          pblh_tke = max(pblh_tke, zi(2))
        end if
        k = k+1
        if (k == nlev - 1) pblh_tke = zi(2) !exit safeguard
      end do

      ! Blend the two PBLH types.
      wgt = 0.5_r8 * tanh((pblh - sbl_lim) / sbl_damp) + 0.5_r8
      pblh = pblh_tke * (1 - wgt) + pblh * wgt

  end subroutine get_pblh

end module pbl_ysu_mod
