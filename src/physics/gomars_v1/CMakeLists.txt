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
  gomars_v1_const_mod.F90
  gomars_v1_driver_mod.F90
  gomars_v1_namelist_mod.F90
  gomars_v1_objects_mod.F90
  gomars_v1_tracers_mod.F90
  gomars_v1_types_mod.F90
  gomars_v1_orbit_mod.F90
  gomars_v1_pbl_mod.F90
  gomars_v1_lsm_mod.F90
  gomars_v1_mp_mod.F90
  gomars_v1_damp_mod.F90
  gomars_v1_rad_mod.F90
  nasa/setspi.F90
  nasa/setspv.F90
  nasa/setrad.F90
  nasa/lagrange.F90
  nasa/laginterp.F90
  nasa/ini_optdst.F90
  nasa/ini_optcld.F90
  nasa/dsolflux.F90
  nasa/funcd.F90
  nasa/newtg.F90
  nasa/tempgr.F90
  nasa/solarza.F90
  nasa/potemp1.F90
  nasa/coldair.F90
  nasa/potemp2.F90
  nasa/fillpt.F90
  nasa/opt_dst.F90
  nasa/opt_cld.F90
  nasa/tpindex.F90
  nasa/optcv.F90
  nasa/getdetau.F90
  nasa/gfluxv.F90
  nasa/sfluxv.F90
  nasa/dtridgl.F90
  nasa/dsolver.F90
  nasa/optci.F90
  nasa/gfluxi.F90
  nasa/sfluxi.F90
  nasa/eddycoef.F90
  nasa/newpbl.F90
)