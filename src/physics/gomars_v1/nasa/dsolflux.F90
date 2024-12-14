subroutine dsolflux(sol, acosz, gweight, fzerov, detau, solar_sfc_dn)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  ! Calculate the direct surface solar flux (solar flux at the surface due to direct solar beam.

  use gomars_v1_const_mod

  implicit none

  real(r8), intent(in ) :: sol(nspectv)
  real(r8), intent(in ) :: acosz
  real(r8), intent(in ) :: gweight(ngauss)
  real(r8), intent(in ) :: fzerov(nspectv)
  real(r8), intent(in ) :: detau(nspectv,ngauss)
  real(r8), intent(out) :: solar_sfc_dn

  integer is, ig
  real(r8) factor

  solar_sfc_dn = 0
  do is = 1, nspectv
    factor = acosz * sol(is)
    do ig = 1, ngauss - 1
      if (detau(is,ig) <= 5) then
        solar_sfc_dn = solar_sfc_dn + factor * exp(-detau(is,ig) / acosz) * gweight(ig) * (1 - fzerov(is))
      end if
    end do
    ig = ngauss
    if (detau(is,ig) <= 5) then
      solar_sfc_dn = solar_sfc_dn + factor * exp(-detau(is,ig) / acosz) * fzerov(ig)
    end if
  end do

end subroutine dsolflux