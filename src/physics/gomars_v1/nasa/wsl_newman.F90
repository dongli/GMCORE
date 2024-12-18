subroutine wsl_newman(rho_sfc, taux, tauy, co2ice_sfc, dstflx_wsl)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod
  use gomars_v1_namelist_mod

  implicit none

  real(r8), intent(in ) :: rho_sfc
  real(r8), intent(in ) :: taux
  real(r8), intent(in ) :: tauy
  real(r8), intent(in ) :: co2ice_sfc
  real(r8), intent(out) :: dstflx_wsl

  real(r8) tau
  real(r8) salt_flux

  tau = sqrt(taux**2 + tauy**2)
  ! No lifting if wind stress is below threshold or CO2 ice is present.
  if (tau <= tau_thresh .or. co2ice_sfc > 0) then
    dstflx_wsl = 0
  else
    ! Calculate the horizontal saltation flux. Refer to Eq. (9) in Newman et al. (2002).
    salt_flux  = max(0.0_r8, 2.61_r8 / g / sqrt(rho_sfc) * tau**1.5_r8 * &
                 (1 - sqrt(tau_thresh / tau)) * &
                 (1 + sqrt(tau_thresh / tau))**2)
    ! Connect the vertical dust lifting flux with the horizontal saltation flux.
    dstflx_wsl = alpha_n * salt_flux
  end if

end subroutine wsl_newman