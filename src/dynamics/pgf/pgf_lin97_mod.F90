! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module pgf_lin97_mod

  use flogger
  use const_mod
  use namelist_mod
  use latlon_parallel_mod
  use block_mod
  use tracer_mod

  implicit none

contains

  subroutine pgf_lin97_run(block, dstate, dtend)

    type(block_type ), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    type(dtend_type ), intent(inout) :: dtend

    real(r8) dpk24, dpk13, dgz13, dgz42, dpp24, dpp13, dph24, dph13, L, tmp
    integer i, j, k

    !                    o
    !                   /|
    !                  / |
    !                 /  |
    !                /   |
    !   o-----------/------------o
    !   |          /|            |
    !   |         / |            |
    !   |        /  |            |
    !   |       /   |            |
    !   |      o    |            |
    !   o------|    -------------o
    !          |   /
    !          |  /
    !          | /
    !          |/
    !          o
    associate (mesh    => block%mesh          , & ! in
               qm      => tracers(block%id)%qm, & ! in
               pkh_lev => block%aux%pkh_lev   , & ! in
               ph_lev  => dstate%ph_lev       , & ! in
               gz_lev  => dstate%gz_lev       , & ! in
               p_lev   => dstate%p_lev        , & ! in
               dudt    => dtend%dudt          , & ! out
               dvdt    => dtend%dvdt          )   ! out
    do k = mesh%half_kds, mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          pkh_lev%d(i,j,k) = ph_lev%d(i,j,k)**rd_o_cpd
        end do
      end do
    end do
    if (hydrostatic) then
      do k = mesh%full_kds, mesh%full_kde
        !
        !   4             3
        ! i,j,k        i+1,j,k
        !   o-------------o
        !   |             |
        !   |             |
        !   |    i,j,k    |
        !   |             |
        !   |             |
        !   o-------------o
        ! i,j,k+1      i+1,j,k+1  --> east
        !   1             2
        !
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            L = 1 + 0.5_r8 * (qm%d(i,j,k) + qm%d(i+1,j,k))
            dgz13 = gz_lev %d(i  ,j,k+1) - gz_lev %d(i+1,j,k  )
            dgz42 = gz_lev %d(i  ,j,k  ) - gz_lev %d(i+1,j,k+1)
            dpk24 = pkh_lev%d(i+1,j,k+1) - pkh_lev%d(i  ,j,k  )
            dpk13 = pkh_lev%d(i  ,j,k+1) - pkh_lev%d(i+1,j,k  )
            tmp = (dgz13 * dpk24 + dgz42 * dpk13) / (dpk24 + dpk13) / mesh%de_lon(j) / L
            dudt%d(i,j,k) = dudt%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
            dtend%dudt_pgf%d(i,j,k) = tmp
#endif
          end do
        end do
        !
        !   4             3
        ! i,j,k        i,j+1,k
        !   o-------------o
        !   |             |
        !   |             |
        !   |    i,j,k    |
        !   |             |
        !   |             |
        !   o-------------o
        ! i,j,k+1      i,j+1,k+1  --> north
        !   1             2
        !
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            L = 1 + 0.5_r8 * (qm%d(i,j,k) + qm%d(i,j+1,k))
            dgz13 = gz_lev %d(i,j  ,k+1) - gz_lev %d(i,j+1,k  )
            dgz42 = gz_lev %d(i,j  ,k  ) - gz_lev %d(i,j+1,k+1)
            dpk24 = pkh_lev%d(i,j+1,k+1) - pkh_lev%d(i,j  ,k  )
            dpk13 = pkh_lev%d(i,j  ,k+1) - pkh_lev%d(i,j+1,k  )
            tmp = (dgz13 * dpk24 + dgz42 * dpk13) / (dpk24 + dpk13) / mesh%de_lat(j) / L
            dvdt%d(i,j,k) = dvdt%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
            dtend%dvdt_pgf%d(i,j,k) = tmp
#endif
          end do
        end do
      end do
    else if (nonhydrostatic) then
      call wait_halo(p_lev)
      do k = mesh%full_kds, mesh%full_kde
        !
        !   4             3
        ! i,j,k        i+1,j,k
        !   o-------------o
        !   |             |
        !   |             |
        !   |    i,j,k    |
        !   |             |
        !   |             |
        !   o-------------o
        ! i,j,k+1      i+1,j,k+1  --> east
        !   1             2
        !
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            L = 1 + 0.5_r8 * (qm%d(i,j,k) + qm%d(i+1,j,k))
            dgz13 = gz_lev %d(i  ,j,k+1) - gz_lev %d(i+1,j,k  )
            dgz42 = gz_lev %d(i  ,j,k  ) - gz_lev %d(i+1,j,k+1)
            dpk24 = pkh_lev%d(i+1,j,k+1) - pkh_lev%d(i  ,j,k  )
            dpk13 = pkh_lev%d(i  ,j,k+1) - pkh_lev%d(i+1,j,k  )
            dph24 = ph_lev %d(i+1,j,k+1) - ph_lev %d(i  ,j,k  )
            dph13 = ph_lev %d(i  ,j,k+1) - ph_lev %d(i+1,j,k  )
            dpp24 = p_lev  %d(i+1,j,k+1) - p_lev  %d(i  ,j,k  ) - dph24
            dpp13 = p_lev  %d(i  ,j,k+1) - p_lev  %d(i+1,j,k  ) - dph13
            tmp = (                                               &
              (dgz13 * dpk24 + dgz42 * dpk13) / (dpk24 + dpk13) + &
              (dgz13 * dpp24 + dgz42 * dpp13) / (dph24 + dph13)   & ! Nonhydrostatic part
            ) / mesh%de_lon(j) / L
            dudt%d(i,j,k) = dudt%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
            dtend%dudt_pgf%d(i,j,k) = tmp
#endif
          end do
        end do
        !
        !   4             3
        ! i,j,k        i,j+1,k
        !   o-------------o
        !   |             |
        !   |             |
        !   |    i,j,k    |
        !   |             |
        !   |             |
        !   o-------------o
        ! i,j,k+1      i,j+1,k+1  --> north
        !   1             2
        !
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            L = 1 + 0.5_r8 * (qm%d(i,j,k) + qm%d(i,j+1,k))
            dgz13 = gz_lev %d(i,j  ,k+1) - gz_lev %d(i,j+1,k  )
            dgz42 = gz_lev %d(i,j  ,k  ) - gz_lev %d(i,j+1,k+1)
            dpk24 = pkh_lev%d(i,j+1,k+1) - pkh_lev%d(i,j  ,k  )
            dpk13 = pkh_lev%d(i,j  ,k+1) - pkh_lev%d(i,j+1,k  )
            dph24 = ph_lev %d(i,j+1,k+1) - ph_lev %d(i,j  ,k  )
            dph13 = ph_lev %d(i,j  ,k+1) - ph_lev %d(i,j+1,k  )
            dpp24 = p_lev  %d(i,j+1,k+1) - p_lev  %d(i,j  ,k  ) - dph24
            dpp13 = p_lev  %d(i,j  ,k+1) - p_lev  %d(i,j+1,k  ) - dph13
            tmp = (                                               &
              (dgz13 * dpk24 + dgz42 * dpk13) / (dpk24 + dpk13) + &
              (dgz13 * dpp24 + dgz42 * dpp13) / (dph24 + dph13)   & ! Nonhydrostatic part
            ) / mesh%de_lat(j) / L
            dvdt%d(i,j,k) = dvdt%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
            dtend%dvdt_pgf%d(i,j,k) = tmp
#endif
          end do
        end do
      end do
    end if
    end associate

  end subroutine pgf_lin97_run

end module pgf_lin97_mod
