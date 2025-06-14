! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module pgf_swm_mod

  use const_mod
  use block_mod

  implicit none

contains

  subroutine pgf_swm_run(block, dstate, dtend)

    type(block_type ), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    type(dtend_type ), intent(inout) :: dtend

    real(r8) tmp
    integer i, j, k

    associate (mesh => block%mesh, &
               gz   => dstate%gz , & ! in
               dudt => dtend%dudt, & ! out
               dvdt => dtend%dvdt)   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          tmp = -(gz%d(i+1,j,k) - gz%d(i,j,k)) / mesh%de_lon(j)
          dudt%d(i,j,k) = dudt%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
          dtend%dudt_pgf%d(i,j,k) = tmp
#endif
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          tmp = -(gz%d(i,j+1,k) - gz%d(i,j,k)) / mesh%de_lat(j)
          dvdt%d(i,j,k) = dvdt%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
          dtend%dvdt_pgf%d(i,j,k) = tmp
#endif
        end do
      end do
    end do
    end associate

  end subroutine pgf_swm_run

end module pgf_swm_mod
