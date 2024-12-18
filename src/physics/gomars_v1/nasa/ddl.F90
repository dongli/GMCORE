subroutine ddl(ps, ht_pbl, ptop_pbl, dstflx_ddl)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_namelist_mod

  implicit none

  real(r8), intent(in   ) :: ps
  real(r8), intent(in   ) :: ht_pbl
  real(r8), intent(in   ) :: ptop_pbl
  real(r8), intent(  out) :: dstflx_ddl

  real(r8) eta ! Thermodynamic efficiency

  if (ht_pbl > 0) then
    eta = 1 - (ps**(rd_o_cpd + 1) - ptop_pbl**(rd_o_cpd + 1)) / &
              ((ps - ptop_pbl) * (rd_o_cpd + 1) * ps**rd_o_cpd)
    dstflx_ddl = max(0.0_r8, alpha_d * eta * ht_pbl)
  end if

end subroutine ddl
