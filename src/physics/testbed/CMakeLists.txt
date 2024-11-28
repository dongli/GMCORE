cmake_policy(SET CMP0076 NEW)

target_sources(gmcore PRIVATE
  testbed_const_mod.F90
  testbed_types_mod.F90
  testbed_objects_mod.F90
  testbed_output_mod.F90
  testbed_pbl_mod.F90
  testbed_driver_mod.F90
)