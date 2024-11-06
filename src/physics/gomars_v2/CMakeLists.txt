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