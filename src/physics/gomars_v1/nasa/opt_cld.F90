subroutine opt_cld(q, pl, qxv, qsv, gv, qxi, qsi, gi, qextrefcld, taurefcld, taucld)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_tracers_mod
  use gomars_v1_mp_mod

  implicit none

  real(r8), intent(in ) :: q         (nlev,ntracers)
  real(r8), intent(in ) :: pl        (2*nlev+3)
  real(r8), intent(out) :: qxv       (2*nlev+4,nspectv)
  real(r8), intent(out) :: qsv       (2*nlev+4,nspectv)
  real(r8), intent(out) :: gv        (2*nlev+4,nspectv)
  real(r8), intent(out) :: qxi       (2*nlev+4,nspecti)
  real(r8), intent(out) :: qsi       (2*nlev+4,nspecti)
  real(r8), intent(out) :: gi        (2*nlev+4,nspecti)
  real(r8), intent(out) :: qextrefcld(2*nlev+4)
  real(r8), intent(out) :: taurefcld (2*nlev+4)
  real(r8), intent(out) :: taucld    (2)

  integer n, k, l, i, is, ib
  integer irap
  real(r8) dev2, cst
  real(r8) mantletocore
  real(r8) Mo, No, Rs, Ao
  real(r8) surf(nbin_rt)

  n = 2 * nlev + 3

  do l = 1, n + 1
    qextrefcld(l) = 1
    taurefcld (l) = 0
  end do

  do is = 1, nspectv
    do l = 1, n + 1
      qxv(l,is) = qextrefcld(l)
      qsv(l,is) = qextrefcld(l) * 0.99_r8
      gv (l,is) = 0
    end do
  end do

  do is = 1, nspecti
    do l = 1, n + 1
      qxi(l,is) = qextrefcld(l)
      qsi(l,is) = qextrefcld(l) * 0.99_r8
      gi (l,is) = 0
    end do
  end do

  dev2 = 1.0_r8 / (sqrt2 * dev_ice)

  taucld = 0

  do k = 1, nlev
    if (q(k,iMa_cld) + q(k,iMa_cor) > 1.0e-7 .and. q(k,iNb_cld) > 1) then
      ! Determine the ratio of dust core over that of ice mantle.
      mantletocore = (q(k,iMa_cor) / dpden_dt / (q(k,iMa_cld) / dpden_ice + q(k,iMa_cor) / dpden_dt))**athird
      ! Find the index to which corresponds the optical properties of the core to mantle ratio.
      ! Those properties were determined off-line using the Toon and Ackerman coated spheres code.
      irap = nratio
      do i = 1, nratio
        if (mantletocore < cor_ratio(i) .and. i /= 1) then
          irap = i - 1
          exit
        else if (mantletocore == cor_ratio(i)) then
          irap = i
          exit
        else if (mantletocore < cor_ratio(i) .and. i == 1) then
          irap = i
          exit
        end if
      end do
      ! Calculate the cross-section mean radius (Rs) of the log-normal distribution.
      Mo = q(k,iMa_cld) + q(k,iMa_cor)
      No = q(k,iNb_cld)
      cst = 0.75_r8 / (pi * (q(k,iMa_cld) / dpden_ice + q(k,iMa_cor) / dpden_dt) / Mo)
      Rs = (Mo / No * cst)**athird * exp(-0.5_r8 * dev_ice**2)
      ! Calculate the total cross sectional area (Ao) of water ice particles.
      Ao = No * pi * Rs**2
      ! Define the cross-section weighted distribution, i.e. surface/size bin. Change Rs to Reff.
      Rs = Rs * exp(1.5_r8 * dev_ice**2)
      Rs = 1.0_r8 / min(max(Rs, 1.0e-7_r8), 100.0e-6_r8)
      do ib = 1, nbin_rt
        surf(ib) = 0.5_r8 * (erf(log(radb_rt(ib+1) * Rs) * dev2) - erf(log(radb_rt(ib) * Rs) * dev2))
      end do
      ! Calculate the averaged values of the optical properties for the whole distribution.
      do l = 2 * k + 2, 2 * k + 3
        do is = 1, nspectv
          qxv(l,is) = 0
          qsv(l,is) = 0
          do ib = 1, nbin_rt
            qxv(l,is) = qxv(l,is) + surf(ib) * qextv_cld (irap,ib,is)
            qsv(l,is) = qsv(l,is) + surf(ib) * qscatv_cld(irap,ib,is)
            gv (l,is) = gv (l,is) + surf(ib) * gv_cld    (irap,ib,is)
          end do
          qsv(l,is) = min(qsv(l,is), 0.99999_r8 * qxv(l,is))
        end do
        do is = 1, nspecti
          qxi(l,is) = 0
          qsi(l,is) = 0
          do ib = 1, nbin_rt
            qxi(l,is) = qxi(l,is) + surf(ib) * qexti_cld (irap,ib,is)
            qsi(l,is) = qsi(l,is) + surf(ib) * qscati_cld(irap,ib,is)
            gi (l,is) = gi (l,is) + surf(ib) * gi_cld    (irap,ib,is)
          end do
          qsi(l,is) = min(qsi(l,is), 0.99999_r8 * qxi(l,is))
        end do
        qextrefcld(l) = qxv(l,nrefv)
        taurefcld (l) = Ao * qextrefcld(l) * (pl(l) - pl(l-1)) / g
        ! For diagnostics: cloud opacity at reference wavelengths (vis and ir).
        taucld(1) = taucld(1) + taurefcld(l)
        taucld(2) = taucld(2) + taurefcld(l) / qextrefcld(l) * (qxi(l,4) - qsi(l,4))
      end do
    end if
  end do

end subroutine opt_cld