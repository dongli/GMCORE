subroutine convect(p, p_lev, dp_dry, pk, pt, pt_lev, q, ptop_pbl, ptcon)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  ! This subroutine checks for (superadiabatic) potential temperature
  ! instabilities between adjacent layers. Any instabilities are resolved by the
  ! process of atmospheric convection. The extent and the potential
  ! temperature of the convective zone are also determined. The convective
  ! zone starts at the surface and includes all layers that exchange heat with
  ! the surface through convection. The pressure of PBL top is also determined.

  use gomars_v1_const_mod

  implicit none

  real(r8), intent(in   ) :: p      (nlev)
  real(r8), intent(in   ) :: p_lev  (nlev+1)
  real(r8), intent(in   ) :: dp_dry (nlev)
  real(r8), intent(in   ) :: pk     (nlev)
  real(r8), intent(inout) :: pt     (nlev)
  real(r8), intent(inout) :: pt_lev (nlev+1)
  real(r8), intent(inout) :: q      (nlev,ntracers)
  real(r8), intent(  out) :: ptop_pbl
  real(r8), intent(  out) :: ptcon

  integer k, k2, k3, m
  real(r8) wgt(nlev), wgt_sum, dp_sum
  real(r8) qcon(ntracers)
  real(r8) pcexp

  ! wgt is a weighting factor for convection. wgt * pt is proportional to the
  ! heat energy of the layer.
  do k = 1, nlev
    wgt(k) = dp_dry(k) * pk(k)
  end do

  ! Each time through the following loop adds another layer to the set and then
  ! makes the atmosphere stable for that set.
  ! If levels k and k+1 are stable, the whole set of layers is stable. This is
  ! because if k = 1, there are no other layers, and if k > 4, the previous time
  ! through this loop left pt(k) non-decreasing in altitude for levels 1 to k.
  loop_k: do k = 1, nlev - 1
    if (pt(k) > pt(k+1)) then
      ! If stable, no adjustments are needed.
      if (k + 1 /= nlev) cycle loop_k
      ! Check pt at layer boundary (instead of layer midpoint) at top of
      ! convective zone. Adjust pcon if pt not non-decreasing here.
      if (pt_lev(nlev) > pt(nlev)) then
        ptcon    = (pt_lev(nlev+1) + pt(nlev)) * 0.5_r8
        ptop_pbl = p(nlev)
      else
        ! Shallow convective layer presents.
        ptcon    = pt(nlev)
        pcexp    = 2 * (ptcon - pt_lev(nlev)) / (pt(nlev-1) - pt_lev(nlev))
        ptop_pbl = p_lev(nlev) * (p(nlev-1) / p_lev(nlev))**pcexp
        ptop_pbl = max(ptop_pbl, p_lev(nlev-1))
      end if
      cycle loop_k
    else
      wgt_sum = wgt   (k+1)
      dp_sum  = dp_dry(k+1)
      qcon    = q     (k+1,:)
      ptcon   = pt    (k+1)
      loop_k2: do k2 = k, 1, -1
        if (pt(k2) <= ptcon) then
          ! Found new instabilities.
          ! This 'if' statement presents roundoff.
          if (pt(k2) /= ptcon) ptcon = (pt(k2) * wgt(k2) + ptcon * wgt_sum) / (wgt(k2) + wgt_sum)
          ! Mix tracers.
          do m = 1, ntracers
            if (q(k2,m) /= qcon(m)) qcon(m) = (q(k2,m) * dp_dry(k2) + qcon(m) * dp_sum) / (dp_dry(k2) + dp_sum)
          end do
          ! Adjust the layers from k2 to k.
          loop_k3: do k3 = k2, k + 1
            pt    (k3  ) = ptcon
            pt_lev(k3+1) = ptcon
            q     (k3,:) = qcon
          end do loop_k3
          ptop_pbl = p_lev(k2)
          wgt_sum  = wgt_sum + wgt(k2)
          dp_sum   = dp_sum  + dp_dry(k2)
        else
          ! Found no more instabilities. Only need to check for pcon adjustment.
          ! Then exit the look_k2.
          if (pt_lev(k2+1) >= ptcon) cycle loop_k
          if (k + 1        /= nlev ) cycle loop_k
          ! Check pt at layer boundary (instead of layer midpoint) at top of
          ! convective zone.  Adjust pcon if pt not non-decreasing here.
          pcexp    = 2 * (ptcon - pt_lev(k2+1)) / (pt(k2) - pt_lev(k2+1))
          ptop_pbl = p_lev(k2+1) * (p(k2) / p_lev(k2+1))**pcexp
          ptop_pbl = max(ptop_pbl, p_lev(k2))
          cycle loop_k
        end if
      end do loop_k2
    end if
  end do loop_k

end subroutine convect