# ==============================================================================
# This file is part of GMCORE since 2019.
#
# GMCORE is a dynamical core for atmospheric model.
#
# GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY. You may contact authors for helping or cooperation.
# ==============================================================================

cmake_minimum_required(VERSION 3.5)

project(gmcore LANGUAGES Fortran C)

option(R4 "Use single precision" OFF)
option(R16 "Use quadruple precision" OFF)
option(ENSURE_ORDER "Ensure computation order in parallel" OFF)
option(USE_GPTL "Use GPTL to profile performance" OFF)
option(USE_BULK_AERO "Use bulk aerosol model" ON)
option(USE_MODAL_AERO "Use modal aerosol model" OFF)
option(USE_CHEM "Use chemistry" OFF)
option(USE_CAM "Use CAM physics" OFF)
option(USE_CLM "Use CLM land model" OFF)
option(OUTPUT_H1_DTEND "Output each dynamics tendency in h1 file" OFF)

if (R4)
  add_definitions(-DREAL_KIND=4)
  set(CPPFLAGS ${CPPFLAGS} -DREAL_KIND=4)
  message(STATUS "Use single precision.")
elseif (R16)
  add_definitions(-DREAL_KIND=16)
  set(CPPFLAGS ${CPPFLAGS} -DREAL_KIND=16)
  message(STATUS "Use quadruple precision.")
else ()
  add_definitions(-DREAL_KIND=8)
endif ()
if (ENSURE_ORDER)
  add_definitions(-DENSURE_ORDER)
  set(CPPFLAGS ${CPPFLAGS} -DENSURE_ORDER)
  message(STATUS "Ensure computation order in parallel.")
endif ()
if (OUTPUT_H1_DTEND)
  add_definitions(-DOUTPUT_H1_DTEND)
endif ()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_definitions(-DFC_IS_GNU)
  list(APPEND fortran_flags
    -ffree-line-length-none
    -fno-range-check
    -ffpe-summary=none
  )
  if (CMAKE_Fortran_COMPILER_VERSION GREATER 10)
    list(APPEND fortran_flags
      -fallow-argument-mismatch
      -fallow-invalid-boz
    )
  endif ()
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    list(APPEND fortran_flags -O0 -fcheck=all -fno-check-array-temporaries -g -fbacktrace)
  else ()
    list(APPEND fortran_flags -Ofast)
    if (ENSURE_ORDER)
      list(APPEND fortran_flags -fcheck=bounds)
    endif ()
  endif ()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  add_definitions(-DFC_IS_INTEL)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    list(APPEND fortran_flags
      -no-wrap-margin
      -O0
      -check bounds
      -g
      -traceback
      -fp-model precise
    )
  else ()
    list(APPEND fortran_flags
      -no-wrap-margin
      -O3
      -align array64byte
    )
    if (ENSURE_ORDER)
      list(APPEND fortran_flags -fp-model precise)
    endif ()
  endif ()
endif ()

find_package(MPI)

if (DEFINED ENV{NETCDF} AND (NOT DEFINED ENV{NETCDF_ROOT}))
  set(ENV{NETCDF_ROOT} $ENV{NETCDF})
endif ()
if (DEFINED ENV{NETCDF_ROOT})
  include_directories("$ENV{NETCDF_ROOT}/include")
  link_directories("$ENV{NETCDF_ROOT}/lib")
else ()
  find_package(PkgConfig)
  if (PKG_CONFIG_FOUND)
    pkg_search_module(NETCDF REQUIRED netcdf)
    if (NETCDF_FOUND)
      include_directories(${NETCDF_INCLUDEDIR})
      link_directories(${NETCDF_LIBDIR})
    else ()
      message(FATAL_ERROR "Pkg-config could not find netcdf library!")
    endif ()
  else ()
    message(FATAL_ERROR "Unable to find pkg-config library!")
  endif ()
endif ()
set(EXTERNAL_LIBS netcdff netcdf)

if (DEFINED ENV{GPTL} AND (NOT DEFINED ENV{GPTL_ROOT}))
  set(ENV{GPTL_ROOT} $ENV{GPTL})
endif ()
if (USE_GPTL AND DEFINED ENV{GPTL_ROOT})
  message(STATUS "Found GPTL in $ENV{GPTL_ROOT}.")
  add_definitions(-DHAS_GPTL)
  set(HAS_GPTL TRUE)
  include_directories("$ENV{GPTL_ROOT}/include")
  link_directories("$ENV{GPTL_ROOT}/lib")
  list(APPEND EXTERNAL_LIBS gptlf gptl)
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    list(APPEND EXTERNAL_LIBS iomp5)
  endif ()
endif ()

if (DEFINED ENV{PAPI} AND (NOT DEFINED ENV{PAPI_ROOT}))
  set(ENV{PAPI_ROOT} $ENV{PAPI})
endif ()
if (DEFINED ENV{PAPI_ROOT})
  message(STATUS "Found PAPI in $ENV{PAPI_ROOT}.")
  add_definitions(-DHAS_PAPI)
  set(HAS_PAPI TRUE)
  include_directories("$ENV{PAPI_ROOT}/include")
  link_directories("$ENV{PAPI_ROOT}/lib")
  list(APPEND EXTERNAL_LIBS papi)
endif ()

if (DEFINED ENV{PIO} AND (NOT DEFINED ENV{PIO_ROOT}))
  set(ENV{PIO_ROOT} $ENV{PIO})
endif ()
if (USE_CAM AND DEFINED ENV{PIO_ROOT})
  message(STATUS "Found PIO in $ENV{PIO_ROOT}.")
  add_definitions(-DHAS_PIO)
  set(HAS_PIO TRUE)
  include_directories("$ENV{PIO_ROOT}/include")
  link_directories("$ENV{PIO_ROOT}/lib")
  list(APPEND EXTERNAL_LIBS piof pioc)
endif ()

if (DEFINED ENV{LAPACK} AND (NOT DEFINED ENV{LAPACK_ROOT}))
  set(ENV{LAPACK_ROOT} $ENV{LAPACK})
endif ()
if (DEFINED ENV{LAPACK_ROOT})
  message(STATUS "Found LAPACK in $ENV{LAPACK_ROOT}.")
  add_definitions(-DHAS_LAPACK)
  set(HAS_LAPACK TRUE)
  include_directories("$ENV{LAPACK_ROOT}/include")
  if (EXISTS "$ENV{LAPACK_ROOT}/lib")
    link_directories("$ENV{LAPACK_ROOT}/lib")
  endif ()
  if (EXISTS "$ENV{LAPACK_ROOT}/lib64")
    link_directories("$ENV{LAPACK_ROOT}/lib64")
  endif ()
  list(APPEND EXTERNAL_LIBS lapack blas)
endif ()

if (DEFINED ENV{ESMF} AND (NOT DEFINED ENV{ESMF_ROOT}))
  set(ENV{ESMF_ROOT} $ENV{ESMF})
endif ()
if (HAS_LAPACK AND DEFINED ENV{ESMF_ROOT})
  message(STATUS "Found ESMF in $ENV{ESMF_ROOT}.")
  add_definitions(-DHAS_ESMF)
  set(HAS_ESMF TRUE)
  include_directories("$ENV{ESMF_ROOT}/include")
  link_directories("$ENV{ESMF_ROOT}/lib")
  list(APPEND EXTERNAL_LIBS esmf)
endif ()

if (DEFINED ENV{MCT} AND (NOT DEFINED ENV{MCT_ROOT}))
  set(ENV{MCT_ROOT} $ENV{MCT})
endif ()
if (DEFINED ENV{MCT_ROOT})
  message(STATUS "Found MCT in $ENV{MCT_ROOT}.")
  add_definitions(-DHAS_MCT)
  set(HAS_MCT TRUE)
  include_directories("$ENV{MCT_ROOT}/include")
  link_directories("$ENV{MCT_ROOT}/lib")
  list(APPEND EXTERNAL_LIBS mct mpeu)
endif ()

if (DEFINED ENV{ECCODES} AND (NOT DEFINED ENV{ECCODES_ROOT}))
  set(ENV{ECCODES_ROOT} $ENV{ECCODES})
endif ()
if (DEFINED ENV{ECCODES_ROOT})
  message(STATUS "Found EcCodes in $ENV{ECCODES_ROOT}.")
  add_definitions(-DHAS_ECCODES)
  include_directories("$ENV{ECCODES_ROOT}/include")
  if (EXISTS "$ENV{ECCODES_ROOT}/lib")
    link_directories("$ENV{ECCODES_ROOT}/lib")
  endif ()
  if (EXISTS "$ENV{ECCODES_ROOT}/lib64")
    link_directories("$ENV{ECCODES_ROOT}/lib64")
  endif ()
  list(APPEND EXTERNAL_LIBS eccodes_f90)
else ()
  find_package(ECCODES QUIET)
  if (ECCODES_FOUND)
    list(APPEND EXTERNAL_LIBS eccodes_f90)
  else ()
    message(WARNING "Unable to find EcCodes!")
  endif ()
endif ()

if (DEFINED ENV{MKLROOT})
  add_definitions(-DHAS_MKL)
  include_directories("$ENV{MKLROOT}/include")
  link_directories("$ENV{MKLROOT}/lib")
elseif (DEFINED ENV{MKL_INC} AND DEFINED ENV{MKL_LIB})
  add_definitions(-DHAS_MKL)
  include_directories("$ENV{MKL_INC}")
  link_directories("$ENV{MKL_LIB}")
endif ()

if (HAS_MKL)
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    list(APPEND EXTERNAL_LIBS lapack64)
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    if (CMAKE_Fortran_COMPILER_VERSION GREATER 2021.1)
      list(APPEND fortran_flags -qmkl=sequential)
    else ()
      list(APPEND fortran_flags -mkl=sequential)
    endif ()
  endif ()
endif ()

string(REPLACE ";" " " CMAKE_Fortran_FLAGS "${fortran_flags}")

add_subdirectory(lib/container)
add_subdirectory(lib/datetime)
add_subdirectory(lib/string)
add_subdirectory(lib/fiona)
add_subdirectory(lib/flogger)

if (EXISTS ${PROJECT_SOURCE_DIR}/lib/rrtmgp/CMakeLists.txt)
  message(STATUS "Build RRTMGP radiation codes.")
  add_subdirectory(lib/rrtmgp)
  list(APPEND EXTERNAL_LIBS rrtmgp)
  set(HAS_RRTMGP TRUE)
endif ()

if (EXISTS ${PROJECT_SOURCE_DIR}/lib/noahmp/CMakeLists.txt)
  message(STATUS "Build NoahMP land codes.")
  add_subdirectory(lib/noahmp)
  list(APPEND EXTERNAL_LIBS noahmp)
  set(HAS_NOAHMP TRUE)
endif ()

add_library(gmcore)
target_sources(gmcore PRIVATE
  lib/quadpack.f90
  src/dynamics/adv/adv_batch_mod.F90
  src/dynamics/adv/adv_mod.F90
  src/dynamics/adv/ffsl_mod.F90
  src/dynamics/adv/limiter_mod.F90
  src/dynamics/adv/ppm_mod.F90
  src/dynamics/adv/tvd_mod.F90
  src/dynamics/adv/upwind_mod.F90
  src/dynamics/adv/weno_mod.F90
  src/dynamics/damp/damp_mod.F90
  src/dynamics/damp/div_damp_mod.F90
  src/dynamics/damp/vor_damp_mod.F90
  src/dynamics/damp/laplace_damp_mod.F90
  src/dynamics/damp/smag_damp_mod.F90
  src/dynamics/dynamics_types_mod.F90
  src/dynamics/filter_mod.F90
  src/dynamics/filter_types_mod.F90
  src/dynamics/gmcore_mod.F90
  src/dynamics/interp_mod.F90
  src/dynamics/operators_mod.F90
  src/dynamics/nh_mod.F90
  src/dynamics/pgf/pgf_lin97_mod.F90
  src/dynamics/pgf/pgf_ptb_mod.F90
  src/dynamics/pgf/pgf_mod.F90
  src/dynamics/pgf/pgf_swm_mod.F90
  src/dynamics/time_schemes_mod.F90
  src/dynamics/tracer_mod.F90
  src/dynamics/tracer_types_mod.F90
  src/dynamics/vert_coord/hybrid_coord_ecmwf_mod.F90
  src/dynamics/vert_coord/hybrid_coord_mod.F90
  src/dynamics/vert_coord/hybrid_coord_ncep_mod.F90
  src/dynamics/vert_coord/hybrid_coord_wrf_mod.F90
  src/dynamics/vert_coord/hybrid_coord_test_mod.F90
  src/dynamics/vert_coord/mars_vert_coord_mod.F90
  src/dynamics/vert_coord/sigma_coord_mod.F90
  src/dynamics/vert_coord/smooth_coord_mod.F90
  src/dynamics/vert_coord/vert_coord_mod.F90
  src/meshes/latlon/latlon_mesh_mod.F90
  src/meshes/latlon/latlon_halo_mod.F90
  src/meshes/latlon/latlon_field_types_mod.F90
  src/meshes/latlon/latlon_operators_mod.F90
  src/meshes/latlon/latlon_parallel_types_mod.F90
  src/meshes/latlon/latlon_parallel_zonal_mod.F90
  src/meshes/latlon/latlon_parallel_global_mod.F90
  src/meshes/latlon/latlon_parallel_mod.F90
  src/meshes/latlon/latlon_decomp_mod.F90
  src/meshes/latlon/latlon_interp_mod.F90
  src/physics/dp_coupling_mod.F90
  src/physics/physics_mod.F90
  src/physics/physics_mesh_mod.F90
  src/physics/physics_types_mod.F90
  src/physics/simple_physics/simple_physics_v6.f90
  src/physics/simple_physics/kessler.f90
  src/physics/simple_physics/simple_physics_types_mod.F90
  src/physics/simple_physics/simple_physics_objects_mod.F90
  src/physics/simple_physics/simple_physics_output_mod.F90
  src/physics/simple_physics/simple_physics_driver_mod.F90
  src/physics/test_forcing_mod.F90
  src/prepare/latlon/latlon_bkg_mod.F90
  src/prepare/latlon/latlon_topo_mod.F90
  src/prepare/prepare_mod.F90
  src/prepare/readers/era5_reader_mod.F90
  src/prepare/readers/cam_reader_mod.F90
  src/prepare/readers/fnl_reader_mod.F90
  src/prepare/readers/openmars_reader_mod.F90
  src/prepare/readers/topo_reader_mod.F90
  src/prepare/ref_mod.F90
  src/prepare/vert_interp_mod.F90
  src/tests/aquaplanet_test_mod.F90
  src/tests/adv/dcmip12_test_mod.F90
  src/tests/adv/deform_test_mod.F90
  src/tests/adv/moving_vortices_test_mod.F90
  src/tests/adv/solid_rotation_test_mod.F90
  src/tests/baroclinic_wave_test_mod.F90
  src/tests/dcmip31_test_mod.F90
  src/tests/held_suarez_test_mod.F90
  src/tests/ksp15_test_mod.F90
  src/tests/mars_cold_run_mod.F90
  src/tests/mountain_wave_test_mod.F90
  src/tests/rossby_haurwitz_wave_3d_test_mod.F90
  src/tests/steady_state_pgf_test_mod.F90
  src/tests/steady_state_test_mod.F90
  src/tests/swm/cross_pole_flow_test_mod.F90
  src/tests/swm/jet_zonal_flow_test_mod.F90
  src/tests/swm/mountain_zonal_flow_test_mod.F90
  src/tests/swm/rossby_haurwitz_wave_test_mod.F90
  src/tests/swm/shallow_water_waves_test_mod.F90
  src/tests/swm/splash_test_mod.F90
  src/tests/swm/steady_geostrophic_flow_test_mod.F90
  src/tests/swm/vortex_erosion_test_mod.F90
  src/tests/tropical_cyclone_test_mod.F90
  src/tests/colliding_modons_test_mod.F90
  src/tests/supercell_test_mod.F90
  src/utils/accum_mod.F90
  src/utils/regrid_mod.F90
  src/utils/block_mod.F90
  src/utils/const_mod.F90
  src/utils/debug_mod.F90
  src/utils/formula_mod.F90
  src/utils/gas_mod.F90
  src/utils/history_mod.F90
  src/utils/initial_mod.F90
  src/utils/math_mod.F90
  src/utils/namelist_mod.F90
  src/utils/perf_mod.F90
  src/utils/process_mod.F90
  src/utils/restart_mod.F90
  src/utils/sphere_geometry_mod.F90
  src/utils/time_mod.F90
  src/utils/time_mod.F90
  src/utils/tridiag_mkl_mod.F90
  src/utils/tridiag_mod.F90
  src/utils/tridiag_spk_mod.F90
)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_Fortran_COMPILER_VERSION GREATER 8)
    target_sources(gmcore PRIVATE
      src/meshes/cubed_sphere/cubed_sphere_panel_mod.F90
      src/meshes/cubed_sphere/cubed_sphere_mesh_mod.F90
      src/meshes/cubed_sphere/cubed_sphere_array_mod.F90
    )
  endif ()
else ()
  target_sources(gmcore PRIVATE
    src/meshes/cubed_sphere/cubed_sphere_panel_mod.F90
    src/meshes/cubed_sphere/cubed_sphere_mesh_mod.F90
    src/meshes/cubed_sphere/cubed_sphere_array_mod.F90
  )
endif ()

add_subdirectory(src/physics/gomars_v1)
add_subdirectory(src/physics/gomars_v2)
add_subdirectory(src/physics/testbed)

if (HAS_RRTMGP AND HAS_NOAHMP)
  target_sources(gmcore PRIVATE
    src/physics/wrf/wrf_namelist_mod.F90
    src/physics/wrf/wrf_tracers_mod.F90
    src/physics/wrf/wrf_physics_types_mod.F90
    src/physics/wrf/wrf_objects_mod.F90
    src/physics/wrf/lsm_noahmp_types_mod.F90
    src/physics/wrf/lsm_noahmp_driver_mod.F90
    src/physics/wrf/lsm_driver_mod.F90
    src/physics/wrf/pbl_driver_mod.F90
    src/physics/wrf/rad_rrtmgp_types_mod.F90
    src/physics/wrf/rad_rrtmgp_mod.F90
  )
  add_definitions(-DHAS_WRF)
endif ()
target_sources(gmcore PRIVATE
  src/physics/wrf/pbl_ysu_mod.F90
)

target_link_libraries(gmcore fortran_container fortran_datetime fiona flogger ${EXTERNAL_LIBS})
if (USE_CAM AND HAS_PIO AND HAS_LAPACK AND HAS_ESMF AND HAS_MCT)
  message(STATUS "PIO, LAPACK, ESMF, and MCT are found. Building CAM physics and chemistry.")
  include_directories(src/physics/cam/include)
  add_definitions(-DSPMD)
  add_definitions(-DUSE_CSM_SHARE)
  add_definitions(-DPCNST=33)
  add_definitions(-DN_RAD_CNST=15)
  add_definitions(-DHAVE_IEEE_ARITHMETIC)
  add_definitions(-DUSE_CONTIGUOUS=)
  add_definitions(-DCLUBB_SGS)
  add_definitions(-DCLUBB_CAM)
  add_definitions(-DCLUBB_REAL_TYPE=8)
  set(USE_CLUBB_SGS ON)
  add_definitions(-DMODAL_AERO_3MODE)
  add_definitions(-DHAS_CAM)
  add_subdirectory(src/physics/cam)
  add_subdirectory(src/chemistry)
else ()
  message(WARNING "PIO, LAPACK, ESMF, and MCT are not found. Skipping CAM physics and chemistry.")
  message(STATUS "HAS_PIO = ${HAS_PIO}")
  message(STATUS "HAS_LAPACK = ${HAS_LAPACK}")
  message(STATUS "HAS_ESMF = ${HAS_ESMF}")
  message(STATUS "HAS_MCT = ${HAS_MCT}")
endif ()

add_executable(gmcore_prepare.exe src/drivers/gmcore_prepare.F90)
target_link_libraries(gmcore_prepare.exe gmcore)

add_executable(gmcore_driver.exe src/drivers/gmcore_driver.F90)
target_link_libraries(gmcore_driver.exe gmcore)

add_executable(gmcore_adv_driver.exe src/drivers/gmcore_adv_driver.F90)
target_link_libraries(gmcore_adv_driver.exe gmcore)

if (HAS_ESMF AND USE_CLM)
  add_subdirectory(src/components/lnd/clm)

  add_executable(gmcore_esmf_driver.exe
    src/coupler/esmf/comp_wrapper_mod.F90
    src/coupler/esmf/atm_comp_mod.F90
    src/coupler/esmf/lnd_comp_mod.F90
    src/coupler/esmf/ocn_comp_mod.F90
    src/coupler/esmf/ice_comp_mod.F90
    src/coupler/esmf/cpl_comp_mod.F90
    src/drivers/gmcore_esmf_driver.F90
  )
  target_link_libraries(gmcore_esmf_driver.exe gmcore clm)
endif ()
