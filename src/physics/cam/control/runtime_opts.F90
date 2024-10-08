module runtime_opts

  !-----------------------------------------------------------------------
  !
  ! Provide driver level routine for making calls to the namelist readers
  ! for the infrastructure and the dycore and physics parameterizations.
  !
  !-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private
  save

  public read_namelist

contains

  subroutine read_namelist(nlfilename)

    use cam_initfiles       , only: cam_initfiles_readnl
    use constituents        , only: cnst_readnl
    use phys_grid           , only: phys_grid_readnl
    use chem_surfvals       , only: chem_surfvals_readnl
    use check_energy        , only: check_energy_readnl
    use radiation           , only: radiation_readnl
    use carma_flags_mod     , only: carma_readnl
    use co2_cycle           , only: co2_cycle_readnl
    use spmd_utils          , only: spmd_utils_readnl
    use cam_history         , only: history_readnl
    use physconst           , only: physconst_readnl
    use air_composition     , only: air_composition_readnl
    use physics_buffer      , only: pbuf_readnl
    use phys_control        , only: phys_ctl_readnl
    use wv_saturation       , only: wv_sat_readnl
    use ref_pres            , only: ref_pres_readnl
    use cam3_aero_data      , only: cam3_aero_data_readnl
    use cam3_ozone_data     , only: cam3_ozone_data_readnl
    use dadadj_cam          , only: dadadj_readnl
    use macrop_driver       , only: macrop_driver_readnl
    use microp_driver       , only: microp_driver_readnl
    use microp_aero         , only: microp_aero_readnl
    use subcol              , only: subcol_readnl
    use cloud_fraction      , only: cldfrc_readnl
    use cldfrc2m            , only: cldfrc2m_readnl
    use rk_stratiform       , only: rk_stratiform_readnl
    use unicon_cam          , only: unicon_cam_readnl
    use zm_conv_intr        , only: zm_conv_readnl
    use hk_conv             , only: hkconv_readnl
    use uwshcu              , only: uwshcu_readnl
    use pkg_cld_sediment    , only: cld_sediment_readnl
    use gw_drag             , only: gw_drag_readnl
    use qbo                 , only: qbo_readnl
    use iondrag             , only: iondrag_readnl
    use waccmx_phys_intr    , only: waccmx_phys_ion_elec_temp_readnl
    use phys_debug_util     , only: phys_debug_readnl
    use conv_water          , only: conv_water_readnl
    use rad_constituents    , only: rad_cnst_readnl
    use radiation_data      , only: rad_data_readnl
    use modal_aer_opt       , only: modal_aer_opt_readnl
    use clubb_intr          , only: clubb_readnl
    use chemistry           , only: chem_readnl
    use prescribed_volcaero , only: prescribed_volcaero_readnl
    use prescribed_strataero, only: prescribed_strataero_readnl
    use aerodep_flx         , only: aerodep_flx_readnl
    use solar_data          , only: solar_data_readnl
    use tropopause          , only: tropopause_readnl
    use aoa_tracers         , only: aoa_tracers_readnl
    use prescribed_ozone    , only: prescribed_ozone_readnl
    use prescribed_aero     , only: prescribed_aero_readnl
    use prescribed_ghg      , only: prescribed_ghg_readnl
    use aircraft_emit       , only: aircraft_emit_readnl
    use cospsimulator_intr  , only: cospsimulator_intr_readnl
    use vertical_diffusion  , only: vd_readnl
    use rayleigh_friction   , only: rayleigh_friction_readnl
    use cam_diagnostics     , only: diag_readnl
    use radheat             , only: radheat_readnl
    use rate_diags          , only: rate_diags_readnl
    use qneg_module         , only: qneg_readnl
    use lunar_tides         , only: lunar_tides_readnl

    character(*), intent(in) :: nlfilename

    character(*), parameter :: subname = "read_namelist"

    ! Call subroutines for modules to read their own namelist.
    ! In some cases namelist default values may depend on settings from
    ! other modules, so there may be an order dependence in the following
    ! calls.
    ! ***N.B.*** In particular, physconst_readnl should be called before
    !            the other readnl methods in case that method is used to set
    !            physical constants, some of which are set at runtime
    !            by the physconst_readnl method.
    ! Modules that read their own namelist are responsible for making sure
    ! all processes receive the values.

    call spmd_utils_readnl               (nlfilename)
    call phys_grid_readnl                (nlfilename)
    call air_composition_readnl          (nlfilename)
    call physconst_readnl                (nlfilename)
    !++bee 13 Oct 2015, need to fix the pbuf_global_allocate functionality, then
    !                   can uncomment the pbuf_readnl line
    ! call pbuf_readnl(nlfilename)
    call cam_initfiles_readnl            (nlfilename)
    call cnst_readnl                     (nlfilename)
    call history_readnl                  (nlfilename)
    call chem_surfvals_readnl            (nlfilename)
    call phys_ctl_readnl                 (nlfilename)
    call wv_sat_readnl                   (nlfilename)
    call ref_pres_readnl                 (nlfilename)
    call cam3_aero_data_readnl           (nlfilename)
    call cam3_ozone_data_readnl          (nlfilename)
    call dadadj_readnl                   (nlfilename)
    call macrop_driver_readnl            (nlfilename)
    call microp_driver_readnl            (nlfilename)
    call microp_aero_readnl              (nlfilename)
    call clubb_readnl                    (nlfilename)
    call subcol_readnl                   (nlfilename)
    call cldfrc_readnl                   (nlfilename)
    call cldfrc2m_readnl                 (nlfilename)
    call unicon_cam_readnl               (nlfilename)
    call zm_conv_readnl                  (nlfilename)
    call rk_stratiform_readnl            (nlfilename)
    call hkconv_readnl                   (nlfilename)
    call uwshcu_readnl                   (nlfilename)
    call cld_sediment_readnl             (nlfilename)
    call gw_drag_readnl                  (nlfilename)
    call qbo_readnl                      (nlfilename)
    call lunar_tides_readnl              (nlfilename)
    call iondrag_readnl                  (nlfilename)
    call waccmx_phys_ion_elec_temp_readnl(nlfilename)
    call phys_debug_readnl               (nlfilename)
    call conv_water_readnl               (nlfilename)
    call radiation_readnl                (nlfilename)
    call rad_cnst_readnl                 (nlfilename)
    call rad_data_readnl                 (nlfilename)
    call modal_aer_opt_readnl            (nlfilename)
    call chem_readnl                     (nlfilename)
    call prescribed_volcaero_readnl      (nlfilename)
    call prescribed_strataero_readnl     (nlfilename)
    call solar_data_readnl               (nlfilename)
    call carma_readnl                    (nlfilename)
    call tropopause_readnl               (nlfilename)
    call aoa_tracers_readnl              (nlfilename)
    call aerodep_flx_readnl              (nlfilename)
    call prescribed_ozone_readnl         (nlfilename)
    call prescribed_aero_readnl          (nlfilename)
    call prescribed_ghg_readnl           (nlfilename)
    call co2_cycle_readnl                (nlfilename)
    call aircraft_emit_readnl            (nlfilename)
    call cospsimulator_intr_readnl       (nlfilename)
    call diag_readnl                     (nlfilename)
    call check_energy_readnl             (nlfilename)
    call radheat_readnl                  (nlfilename)
    call vd_readnl                       (nlfilename)
    call rayleigh_friction_readnl        (nlfilename)
    call rate_diags_readnl               (nlfilename)
    call qneg_readnl                     (nlfilename)

  end subroutine read_namelist

end module runtime_opts
