subroutine dsolflux(sol, acosz, gweight, fzerov, detau, directsol)

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  ! Calculate the direct surface solar flux (solar flux at the surface due to direct solar beam.

  use gomars_v1_const_mod

  implicit none

  real(r8), intent(in) :: sol(l_nspectv)
  real(r8), intent(in) :: acosz
  real(r8), intent(in) :: gweight(l_ngauss)
  real(r8), intent(in) :: fzerov(l_nspectv)
  real(r8), intent(in) :: detau(l_nspectv,l_ngauss)
  real(r8), intent(out) :: directsol

  integer is, ig
  real(r8) factor

  directsol = 0
  do is = 1, l_nspectv
    factor = acosz * sol(is)
    do ig = 1, l_ngauss - 1
      if (detau(is,ig) <= 5) then
        directsol = directsol + factor * exp(-detau(is,ig) / acosz) * gweight(ig) * (1 - fzerov(is))
      end if
    end do
    ig = l_ngauss
    if (detau(is,ig) <= 5) then
      directsol = directsol + factor * exp(-detau(is,ig) / acosz) * fzerov(ig)
    end if
  end do

end subroutine dsolflux