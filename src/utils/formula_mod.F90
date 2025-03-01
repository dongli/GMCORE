! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module formula_mod

  use const_mod

  implicit none

  private

  public wet_mixing_ratio
  public dry_mixing_ratio
  public dry_potential_temperature
  public modified_potential_temperature
  public temperature
  public temperature_from_virtual_temperature_and_wet_mixing_ratio
  public temperature_from_density_and_wet_mixing_ratio
  public virtual_temperature
  public virtual_temperature_from_wet_mixing_ratio
  public virtual_temperature_from_density
  public virtual_temperature_from_modified_potential_temperature
  public virtual_temperature_from_geopotential
  public virtual_potential_temperature
  public dry_air_density
  public moist_air_density
  public buoyancy_frequency
  public local_richardson_number
  public water_vapor_saturation_mixing_ratio_mars
  public dewpoint_temperature_mars

contains

  pure elemental real(r8) function wet_mixing_ratio(qv, qm) result(res)

    real(r8), intent(in) :: qv  ! Dry mixing ratio of water vapor (kg kg-1)
    real(r8), intent(in) :: qm  ! Total dry mixing ratio of water vapor and its condensate (kg kg-1)

    res = qv / (1 + qm)

  end function wet_mixing_ratio

  pure elemental real(r8) function dry_mixing_ratio(qv, qm) result(res)

    real(r8), intent(in) :: qv  ! Specific humidity or wet mixing ratio of water vapor (kg kg-1)
    real(r8), intent(in) :: qm  ! Total dry mixing ratio of water vapor and its condensate (kg kg-1)

    res = qv * (1 + qm)

  end function dry_mixing_ratio

  pure elemental real(r8) function dry_potential_temperature(t, p) result(res)

    real(r8), intent(in) :: t   ! Temperature (K)
    real(r8), intent(in) :: p   ! Full pressure (Pa)

    res = t * (p0 / p)**rd_o_cpd

  end function dry_potential_temperature

  pure elemental real(r8) function potential_temperature(t, p) result(res)

    real(r8), intent(in) :: t   ! Temperature (K)
    real(r8), intent(in) :: p   ! Full pressure (Pa)

    res = t * (p0 / p)**rd_o_cpd

  end function potential_temperature

  pure elemental real(r8) function modified_potential_temperature(t, p, qv) result(res)

    real(r8), intent(in) :: t   ! Temperature (K)
    real(r8), intent(in) :: p   ! Full pressure (Pa)
    real(r8), intent(in) :: qv  ! Dry mixing ratio of water vapor (kg kg-1)

    res = t * (p0 / p)**rd_o_cpd * (1 + rv_o_rd * qv)

  end function modified_potential_temperature

  pure elemental real(r8) function temperature(pt, p, qv) result(res)

    real(r8), intent(in) :: pt  ! Modified potential temperature (K)
    real(r8), intent(in) :: p   ! Full pressure (Pa)
    real(r8), intent(in) :: qv  ! Dry mixing ratio of water vapor (kg kg-1)

    res = pt * (p / p0)**rd_o_cpd / (1 + rv_o_rd * qv)

  end function temperature

  pure elemental real(r8) function temperature_from_virtual_temperature_and_wet_mixing_ratio(tv, qv) result(res)

    real(r8), intent(in) :: tv  ! Virtual temperature (K)
    real(r8), intent(in) :: qv  ! Wet mixing ratio of water vapor (kg kg-1)

    res = tv / (1 + (rv_o_rd - 1) * qv)

  end function temperature_from_virtual_temperature_and_wet_mixing_ratio

  pure elemental real(r8) function temperature_from_density_and_wet_mixing_ratio(rho, p, qv) result(res)

    real(r8), intent(in) :: rho ! Air density (kg m-3)
    real(r8), intent(in) :: p   ! Full pressure (Pa)
    real(r8), intent(in) :: qv  ! Wet mixing ratio of water vapor (kg kg-1)

    res = p / rd / rho / (1 + (rv_o_rd - 1) * qv)

  end function temperature_from_density_and_wet_mixing_ratio

  pure elemental real(r8) function virtual_temperature(t, qv, qm) result(res)

    real(r8), intent(in) :: t   ! Temperature (K)
    real(r8), intent(in) :: qv  ! Dry mixing ratio of water vapor (kg kg-1)
    real(r8), intent(in) :: qm  ! Total dry mixing ratio of water vapor and its condensate (kg kg-1)

    res = t * (1 + rv_o_rd * qv) / (1 + qm) ! See (16) in Lauritzen et al. (2018) for details.

  end function virtual_temperature

  pure elemental real(r8) function virtual_temperature_from_wet_mixing_ratio(t, qv) result(res)

    real(r8), intent(in) :: t  ! Temperature (K)
    real(r8), intent(in) :: qv ! Wet mixing ratio of water vapor (kg kg-1)

    res = t * (1 + (rv_o_rd - 1) * qv)

  end function virtual_temperature_from_wet_mixing_ratio

  pure elemental real(r8) function virtual_temperature_from_density(p, rho) result(res)

    real(r8), intent(in) :: p   ! Full pressure (Pa)
    real(r8), intent(in) :: rho ! Air density (kg m-3)

    res = p / rd / rho

  end function virtual_temperature_from_density

  pure elemental real(r8) function virtual_temperature_from_modified_potential_temperature(pt, pk, qm) result(res)

    real(r8), intent(in) :: pt  ! Modified potential temperature (K)
    real(r8), intent(in) :: pk  ! p**(rd/cpd)
    real(r8), intent(in) :: qm  ! Total dry mixing ratio of water vapor and its condensate (kg kg-1)

    res = pt * pk / pk0 / (1 + qm)

  end function virtual_temperature_from_modified_potential_temperature

  pure elemental real(r8) function virtual_temperature_from_geopotential(p1, p2, gz1, gz2) result(res)

    real(r8), intent(in) :: p1  ! Full pressure at level 1 (Pa)
    real(r8), intent(in) :: p2  ! Full pressure at level 2 (Pa)
    real(r8), intent(in) :: gz1 ! Geopotential at level 1 (m2 s-2)
    real(r8), intent(in) :: gz2 ! Geopotential at level 2 (m2 s-2)

    res = - (gz2 - gz1) / (rd_o_cpd * log(p2 / p1))

  end function virtual_temperature_from_geopotential

  pure elemental real(r8) function virtual_potential_temperature(tv, p) result(res)

    real(r8), intent(in) :: tv  ! Virtual temperature (K)
    real(r8), intent(in) :: p   ! Full pressure (Pa)

    res = tv * (p0 / p)**rd_o_cpd

  end function virtual_potential_temperature

  pure elemental real(r8) function dry_air_density(t, p) result(res)

    real(r8), intent(in) :: t   ! Temperature (K)
    real(r8), intent(in) :: p   ! Full pressure (Pa)

    res = p / rd / t

  end function dry_air_density

  pure elemental real(r8) function moist_air_density(t, p, qv, qm) result(res)

    real(r8), intent(in) :: t   ! Temperature (K)
    real(r8), intent(in) :: p   ! Full pressure (Pa)
    real(r8), intent(in) :: qv  ! Dry mixing ratio of water vapor (kg kg-1)
    real(r8), intent(in) :: qm  ! Total dry mixing ratio of water vapor and its condensate (kg kg-1)

    res = p / rd / virtual_temperature(t, qv, qm)

  end function moist_air_density

  pure elemental real(r8) function buoyancy_frequency(pt1, pt2, z1, z2) result(res)

    real(r8), intent(in) :: pt1
    real(r8), intent(in) :: pt2
    real(r8), intent(in) :: z1
    real(r8), intent(in) :: z2

    res = g * (pt1 - pt2) / (z1 - z2) / (pt1 + pt2) * 2

  end function buoyancy_frequency

  pure elemental real(r8) function local_richardson_number(N2, z1, z2, u1, u2, v1, v2) result(res)

    real(r8), intent(in) :: N2
    real(r8), intent(in) :: z1
    real(r8), intent(in) :: z2
    real(r8), intent(in) :: u1
    real(r8), intent(in) :: u2
    real(r8), intent(in) :: v1
    real(r8), intent(in) :: v2

    real(r8) s2

    s2 = ((u1 - u2)**2 + (v1 - v2)**2) / (z1 - z2)**2

    res = n2 / (s2 + 1.0e-4_r8)

  end function local_richardson_number

  pure elemental real(r8) function water_vapor_saturation_mixing_ratio_mars(t, p) result(res)

    real(r8), intent(in) :: t ! Temperature (K)
    real(r8), intent(in) :: p ! Full pressure (Pa)

    res = 611 * exp(22.5_r8 * (1 - t0_trip / t)) / p * m_h2o / m_co2

  end function water_vapor_saturation_mixing_ratio_mars

  pure elemental real(r8) function dewpoint_temperature_mars(p) result(res)

    real(r8), intent(in) :: p ! Full pressure (Pa)

    res = 3182.48_r8 / (23.3494_r8 - log(p / 100.0_r8))

  end function dewpoint_temperature_mars

end module formula_mod
