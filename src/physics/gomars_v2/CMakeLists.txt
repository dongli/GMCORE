# ==============================================================================
# This file is part of GoMars since 2023.
#
# GoMars is a Martian general circulation model developed in Institute of
# Atmospheric Physics (IAP), Chinese Academy of Sciences (CAS).
#
# GMCORE is a dynamical core for atmospheric model used in GoMars.
#
# GoMars and GMCORE are distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY. You may contact authors for helping or cooperation.
# ==============================================================================

cmake_policy(SET CMP0076 NEW)

target_sources(gmcore PRIVATE
  gomars_v2_const_mod.F90
  gomars_v2_driver_mod.F90
  gomars_v2_lsm_mod.F90
  gomars_v2_mp_mod.F90
  gomars_v2_namelist_mod.F90
  gomars_v2_objects_mod.F90
  gomars_v2_optics_mod.F90
  gomars_v2_orbit_mod.F90
  gomars_v2_output_mod.F90
  gomars_v2_pbl_mod.F90
  gomars_v2_rad_kcoef_mod.F90
  gomars_v2_rad_mod.F90
  gomars_v2_sfc_mod.F90
  gomars_v2_solar_mod.F90
  gomars_v2_spectra_mod.F90
  gomars_v2_tracers_mod.F90
  gomars_v2_types_mod.F90
)