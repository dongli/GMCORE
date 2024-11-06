cmake_policy(SET CMP0076 NEW)

target_sources(gmcore PRIVATE
  gomars_v1_const_mod.F90
  gomars_v1_driver_mod.F90
  gomars_v1_namelist_mod.F90
  gomars_v1_objects_mod.F90
  gomars_v1_types_mod.F90
  nasa/nasa_rad_mod.F90
)