module camsrfexch

  !-----------------------------------------------------------------------
  ! Module to handle data that is exchanged between the CAM atmosphere
  ! model and the surface models (land, sea-ice, and ocean).
  !-----------------------------------------------------------------------

  use shr_kind_mod,    only: r8 => shr_kind_r8, r4 => shr_kind_r4
  use constituents,    only: pcnst
  use ppgrid,          only: pcols, begchunk, endchunk
  use phys_grid,       only: get_ncols_p, phys_grid_initialized
  use infnan,          only: posinf, assignment(=)
  use cam_abortutils,  only: endrun
  use cam_logfile,     only: iulog
  use srf_field_check, only: active_Sl_ram1, active_Sl_fv, active_Sl_soilw,                &
                             active_Fall_flxdst1, active_Fall_flxvoc, active_Fall_flxfire, &
                             active_Faxa_nhx, active_Faxa_noy

  implicit none
  private

  public atm2hub_alloc              ! Atmosphere to surface data allocation method
  public hub2atm_alloc              ! Merged hub surface to atmosphere data allocation method
  public atm2hub_deallocate
  public hub2atm_deallocate
  public cam_export
  public cam_out_t                  ! Data from atmosphere
  public cam_in_t                   ! Merged surface data

  !---------------------------------------------------------------------------
  ! This is the data that is sent from the atmosphere to the surface models
  !---------------------------------------------------------------------------

  type cam_out_t
    integer lchnk
    integer ncol
    real(r8), allocatable, dimension(:  ) :: tbot                       ! bot level temperature
    real(r8), allocatable, dimension(:  ) :: zbot                       ! bot level height above surface
    real(r8), allocatable, dimension(:  ) :: topo                       ! surface topographic height (m)
    real(r8), allocatable, dimension(:  ) :: ubot                       ! bot level u wind
    real(r8), allocatable, dimension(:  ) :: vbot                       ! bot level v wind
    real(r8), allocatable, dimension(:,:) :: qbot                       ! bot level specific humidity
    real(r8), allocatable, dimension(:  ) :: pbot                       ! bot level pressure
    real(r8), allocatable, dimension(:  ) :: rho                        ! bot level density
    real(r8), allocatable, dimension(:  ) :: netsw                      !
    real(r8), allocatable, dimension(:  ) :: flwds                      !
    real(r8), allocatable, dimension(:  ) :: precsc                     !
    real(r8), allocatable, dimension(:  ) :: precsl                     !
    real(r8), allocatable, dimension(:  ) :: precc                      !
    real(r8), allocatable, dimension(:  ) :: precl                      !
    real(r8), allocatable, dimension(:  ) :: soll                       !
    real(r8), allocatable, dimension(:  ) :: sols                       !
    real(r8), allocatable, dimension(:  ) :: solld                      !
    real(r8), allocatable, dimension(:  ) :: solsd                      !
    real(r8), allocatable, dimension(:  ) :: thbot                      !
    real(r8), allocatable, dimension(:  ) :: co2prog                    ! prognostic co2
    real(r8), allocatable, dimension(:  ) :: co2diag                    ! diagnostic co2
    real(r8), allocatable, dimension(:  ) :: psl     
    real(r8), allocatable, dimension(:  ) :: bcphiwet                   ! wet deposition of hydrophilic black carbon
    real(r8), allocatable, dimension(:  ) :: bcphidry                   ! dry deposition of hydrophilic black carbon
    real(r8), allocatable, dimension(:  ) :: bcphodry                   ! dry deposition of hydrophobic black carbon
    real(r8), allocatable, dimension(:  ) :: ocphiwet                   ! wet deposition of hydrophilic organic carbon
    real(r8), allocatable, dimension(:  ) :: ocphidry                   ! dry deposition of hydrophilic organic carbon
    real(r8), allocatable, dimension(:  ) :: ocphodry                   ! dry deposition of hydrophobic organic carbon
    real(r8), allocatable, dimension(:  ) :: dstwet1                    ! wet deposition of dust (bin1)
    real(r8), allocatable, dimension(:  ) :: dstdry1                    ! dry deposition of dust (bin1)
    real(r8), allocatable, dimension(:  ) :: dstwet2                    ! wet deposition of dust (bin2)
    real(r8), allocatable, dimension(:  ) :: dstdry2                    ! dry deposition of dust (bin2)
    real(r8), allocatable, dimension(:  ) :: dstwet3                    ! wet deposition of dust (bin3)
    real(r8), allocatable, dimension(:  ) :: dstdry3                    ! dry deposition of dust (bin3)
    real(r8), allocatable, dimension(:  ) :: dstwet4                    ! wet deposition of dust (bin4)
    real(r8), allocatable, dimension(:  ) :: dstdry4                    ! dry deposition of dust (bin4)
    real(r8), pointer    , dimension(:  ) :: nhx_nitrogen_flx => null() ! nitrogen deposition fluxes (kgN/m2/s)
    real(r8), pointer    , dimension(:  ) :: noy_nitrogen_flx => null() ! nitrogen deposition fluxes (kgN/m2/s)
  contains
    procedure :: init  => cam_out_init
    procedure :: clear => cam_out_clear
    final cam_out_final
  end type cam_out_t

  !---------------------------------------------------------------------------
  ! This is the merged state of sea-ice, land and ocean surface parameterizations
  !---------------------------------------------------------------------------

  type cam_in_t
    integer lchnk
    integer ncol
    real(r8), allocatable, dimension(:  ) :: asdir              ! albedo: shortwave, direct
    real(r8), allocatable, dimension(:  ) :: asdif              ! albedo: shortwave, diffuse
    real(r8), allocatable, dimension(:  ) :: aldir              ! albedo: longwave, direct
    real(r8), allocatable, dimension(:  ) :: aldif              ! albedo: longwave, diffuse
    real(r8), allocatable, dimension(:  ) :: lwup               ! longwave up radiative flux
    real(r8), allocatable, dimension(:  ) :: lhf                ! latent heat flux
    real(r8), allocatable, dimension(:  ) :: shf                ! sensible heat flux
    real(r8), allocatable, dimension(:  ) :: wsx                ! surface u-stress (N)
    real(r8), allocatable, dimension(:  ) :: wsy                ! surface v-stress (N)
    real(r8), allocatable, dimension(:  ) :: tref               ! ref height surface air temp
    real(r8), allocatable, dimension(:  ) :: qref               ! ref height specific humidity
    real(r8), allocatable, dimension(:  ) :: u10                ! 10m wind speed
    real(r8), allocatable, dimension(:  ) :: ugust_out          ! gustiness added
    real(r8), allocatable, dimension(:  ) :: u10_with_gust      ! 10m wind speed with gusts added
    real(r8), allocatable, dimension(:  ) :: ts                 ! merged surface temp
    real(r8), allocatable, dimension(:  ) :: sst                ! sea surface temp
    real(r8), allocatable, dimension(:  ) :: snowhland          ! snow depth (liquid water equivalent) over land
    real(r8), allocatable, dimension(:  ) :: snowhice           ! snow depth over ice
    real(r8), allocatable, dimension(:  ) :: fco2_lnd           ! co2 flux from lnd
    real(r8), allocatable, dimension(:  ) :: fco2_ocn           ! co2 flux from ocn
    real(r8), allocatable, dimension(:  ) :: fdms               ! dms flux
    real(r8), allocatable, dimension(:  ) :: landfrac           ! land area fraction
    real(r8), allocatable, dimension(:  ) :: icefrac            ! sea-ice areal fraction
    real(r8), allocatable, dimension(:  ) :: ocnfrac            ! ocean areal fraction
    real(r8), allocatable, dimension(:,:) :: cflx               ! constituent flux (emissions)
    real(r8), allocatable, dimension(:  ) :: ustar              ! atm/ocn saved version of ustar
    real(r8), allocatable, dimension(:  ) :: re                 ! atm/ocn saved version of re
    real(r8), allocatable, dimension(:  ) :: ssq                ! atm/ocn saved version of ssq
    real(r8), pointer    , dimension(:  ) :: ram1     => null() ! aerodynamical resistance (s/m) (pcols)
    real(r8), pointer    , dimension(:  ) :: fv       => null() ! friction velocity (m/s) (pcols)
    real(r8), pointer    , dimension(:  ) :: soilw    => null() ! volumetric soil water (m3/m3)
    real(r8), pointer    , dimension(:,:) :: depvel   => null() ! deposition velocities
    real(r8), pointer    , dimension(:,:) :: dstflx   => null() ! dust fluxes
    real(r8), pointer    , dimension(:,:) :: meganflx => null() ! MEGAN fluxes
    real(r8), pointer    , dimension(:,:) :: fireflx  => null() ! wild fire emissions
    real(r8), pointer    , dimension(:  ) :: fireztop => null() ! wild fire emissions vert distribution top
  contains
    procedure :: init  => cam_in_init
    procedure :: clear => cam_in_clear
    final cam_in_final
  end type cam_in_t

contains

  subroutine cam_out_init(this)

    class(cam_out_t), intent(inout) :: this

    call this%clear()

    allocate(this%tbot    (pcols      ))
    allocate(this%zbot    (pcols      ))
    allocate(this%topo    (pcols      ))
    allocate(this%ubot    (pcols      ))
    allocate(this%vbot    (pcols      ))
    allocate(this%qbot    (pcols,PCNST))
    allocate(this%pbot    (pcols      ))
    allocate(this%rho     (pcols      ))
    allocate(this%netsw   (pcols      ))
    allocate(this%flwds   (pcols      ))
    allocate(this%precsc  (pcols      ))
    allocate(this%precsl  (pcols      ))
    allocate(this%precc   (pcols      ))
    allocate(this%precl   (pcols      ))
    allocate(this%soll    (pcols      ))
    allocate(this%sols    (pcols      ))
    allocate(this%solld   (pcols      ))
    allocate(this%solsd   (pcols      ))
    allocate(this%thbot   (pcols      ))
    allocate(this%co2prog (pcols      ))
    allocate(this%co2diag (pcols      ))
    allocate(this%psl     (pcols      ))
    allocate(this%bcphiwet(pcols      ))
    allocate(this%bcphidry(pcols      ))
    allocate(this%bcphodry(pcols      ))
    allocate(this%ocphiwet(pcols      ))
    allocate(this%ocphidry(pcols      ))
    allocate(this%ocphodry(pcols      ))
    allocate(this%dstwet1 (pcols      ))
    allocate(this%dstdry1 (pcols      ))
    allocate(this%dstwet2 (pcols      ))
    allocate(this%dstdry2 (pcols      ))
    allocate(this%dstwet3 (pcols      ))
    allocate(this%dstdry3 (pcols      ))
    allocate(this%dstwet4 (pcols      ))
    allocate(this%dstdry4 (pcols      ))

  end subroutine cam_out_init

  subroutine cam_out_clear(this)

    class(cam_out_t), intent(inout) :: this

    if (allocated(this%tbot    )) deallocate(this%tbot    )
    if (allocated(this%zbot    )) deallocate(this%zbot    )
    if (allocated(this%topo    )) deallocate(this%topo    )
    if (allocated(this%ubot    )) deallocate(this%ubot    )
    if (allocated(this%vbot    )) deallocate(this%vbot    )
    if (allocated(this%qbot    )) deallocate(this%qbot    )
    if (allocated(this%pbot    )) deallocate(this%pbot    )
    if (allocated(this%rho     )) deallocate(this%rho     )
    if (allocated(this%netsw   )) deallocate(this%netsw   )
    if (allocated(this%flwds   )) deallocate(this%flwds   )
    if (allocated(this%precsc  )) deallocate(this%precsc  )
    if (allocated(this%precsl  )) deallocate(this%precsl  )
    if (allocated(this%precc   )) deallocate(this%precc   )
    if (allocated(this%precl   )) deallocate(this%precl   )
    if (allocated(this%soll    )) deallocate(this%soll    )
    if (allocated(this%sols    )) deallocate(this%sols    )
    if (allocated(this%solld   )) deallocate(this%solld   )
    if (allocated(this%solsd   )) deallocate(this%solsd   )
    if (allocated(this%thbot   )) deallocate(this%thbot   )
    if (allocated(this%co2prog )) deallocate(this%co2prog )
    if (allocated(this%co2diag )) deallocate(this%co2diag )
    if (allocated(this%psl     )) deallocate(this%psl     )
    if (allocated(this%bcphiwet)) deallocate(this%bcphiwet)
    if (allocated(this%bcphidry)) deallocate(this%bcphidry)
    if (allocated(this%bcphodry)) deallocate(this%bcphodry)
    if (allocated(this%ocphiwet)) deallocate(this%ocphiwet)
    if (allocated(this%ocphidry)) deallocate(this%ocphidry)
    if (allocated(this%ocphodry)) deallocate(this%ocphodry)
    if (allocated(this%dstwet1 )) deallocate(this%dstwet1 )
    if (allocated(this%dstdry1 )) deallocate(this%dstdry1 )
    if (allocated(this%dstwet2 )) deallocate(this%dstwet2 )
    if (allocated(this%dstdry2 )) deallocate(this%dstdry2 )
    if (allocated(this%dstwet3 )) deallocate(this%dstwet3 )
    if (allocated(this%dstdry3 )) deallocate(this%dstdry3 )
    if (allocated(this%dstwet4 )) deallocate(this%dstwet4 )
    if (allocated(this%dstdry4 )) deallocate(this%dstdry4 )

    if (associated(this%nhx_nitrogen_flx)) deallocate(this%nhx_nitrogen_flx)
    if (associated(this%noy_nitrogen_flx)) deallocate(this%noy_nitrogen_flx)

  end subroutine cam_out_clear

  subroutine cam_out_final(this)

    type(cam_out_t), intent(inout) :: this

    call this%clear()

  end subroutine cam_out_final

  subroutine cam_in_init(this)

    class(cam_in_t), intent(inout) :: this

    call this%clear()

    allocate(this%asdir        (pcols      ))
    allocate(this%asdif        (pcols      ))
    allocate(this%aldir        (pcols      ))
    allocate(this%aldif        (pcols      ))
    allocate(this%lwup         (pcols      ))
    allocate(this%lhf          (pcols      ))
    allocate(this%shf          (pcols      ))
    allocate(this%wsx          (pcols      ))
    allocate(this%wsy          (pcols      ))
    allocate(this%tref         (pcols      ))
    allocate(this%qref         (pcols      ))
    allocate(this%u10          (pcols      ))
    allocate(this%ugust_out    (pcols      ))
    allocate(this%u10_with_gust(pcols      ))
    allocate(this%ts           (pcols      ))
    allocate(this%sst          (pcols      ))
    allocate(this%snowhland    (pcols      ))
    allocate(this%snowhice     (pcols      ))
    allocate(this%fco2_lnd     (pcols      ))
    allocate(this%fco2_ocn     (pcols      ))
    allocate(this%fdms         (pcols      ))
    allocate(this%landfrac     (pcols      ))
    allocate(this%icefrac      (pcols      ))
    allocate(this%ocnfrac      (pcols      ))
    allocate(this%cflx         (pcols,PCNST))
    allocate(this%ustar        (pcols      ))
    allocate(this%re           (pcols      ))
    allocate(this%ssq          (pcols      ))

  end subroutine cam_in_init

  subroutine cam_in_clear(this)

    class(cam_in_t), intent(inout) :: this

    if (allocated(this%asdir        )) deallocate(this%asdir        )
    if (allocated(this%asdif        )) deallocate(this%asdif        )
    if (allocated(this%aldir        )) deallocate(this%aldir        )
    if (allocated(this%aldif        )) deallocate(this%aldif        )
    if (allocated(this%lwup         )) deallocate(this%lwup         )
    if (allocated(this%lhf          )) deallocate(this%lhf          )
    if (allocated(this%shf          )) deallocate(this%shf          )
    if (allocated(this%wsx          )) deallocate(this%wsx          )
    if (allocated(this%wsy          )) deallocate(this%wsy          )
    if (allocated(this%tref         )) deallocate(this%tref         )
    if (allocated(this%qref         )) deallocate(this%qref         )
    if (allocated(this%u10          )) deallocate(this%u10          )
    if (allocated(this%ugust_out    )) deallocate(this%ugust_out    )
    if (allocated(this%u10_with_gust)) deallocate(this%u10_with_gust)
    if (allocated(this%ts           )) deallocate(this%ts           )
    if (allocated(this%sst          )) deallocate(this%sst          )
    if (allocated(this%snowhland    )) deallocate(this%snowhland    )
    if (allocated(this%snowhice     )) deallocate(this%snowhice     )
    if (allocated(this%fco2_lnd     )) deallocate(this%fco2_lnd     )
    if (allocated(this%fco2_ocn     )) deallocate(this%fco2_ocn     )
    if (allocated(this%fdms         )) deallocate(this%fdms         )
    if (allocated(this%landfrac     )) deallocate(this%landfrac     )
    if (allocated(this%icefrac      )) deallocate(this%icefrac      )
    if (allocated(this%ocnfrac      )) deallocate(this%ocnfrac      )
    if (allocated(this%cflx         )) deallocate(this%cflx         )
    if (allocated(this%ustar        )) deallocate(this%ustar        )
    if (allocated(this%re           )) deallocate(this%re           )
    if (allocated(this%ssq          )) deallocate(this%ssq          )

    if (associated(this%ram1        )) deallocate(this%ram1         )
    if (associated(this%fv          )) deallocate(this%fv           )
    if (associated(this%soilw       )) deallocate(this%soilw        )
    if (associated(this%depvel      )) deallocate(this%depvel       )
    if (associated(this%dstflx      )) deallocate(this%dstflx       )
    if (associated(this%meganflx    )) deallocate(this%meganflx     )
    if (associated(this%fireflx     )) deallocate(this%fireflx      )
    if (associated(this%fireztop    )) deallocate(this%fireztop     )

  end subroutine cam_in_clear

  subroutine cam_in_final(this)

    type(cam_in_t), intent(inout) :: this

    call this%clear()

  end subroutine cam_in_final

  subroutine hub2atm_alloc(cam_in)

    ! Allocate space for the surface to atmosphere data type. And initialize the values.

    use seq_drydep_mod   , only: lnd_drydep, n_drydep
    use shr_megan_mod    , only: shr_megan_mechcomps_n
    use shr_fire_emis_mod, only: shr_fire_emis_mechcomps_n

    type(cam_in_t), pointer :: cam_in(:) ! Merged surface state

    integer c, ierror
    character(*), parameter :: sub = 'hub2atm_alloc'

    if (.not. phys_grid_initialized()) call endrun(sub // ': phys_grid not called yet')
    allocate(cam_in(begchunk:endchunk), stat=ierror)
    if (ierror /= 0) then
      write(iulog, *) sub // ': Allocation error: ', ierror
      call endrun(sub // ': allocation error')
    end if

    do c = begchunk, endchunk
      call cam_in(c)%init()
      if (active_Sl_ram1) then
        allocate(cam_in(c)%ram1(pcols), stat=ierror)
        if (ierror /= 0) call endrun(sub // ': allocation error ram1')
      end if
      if (active_Sl_fv) then
        allocate (cam_in(c)%fv(pcols), stat=ierror)
        if (ierror /= 0) call endrun(sub // ': allocation error fv')
      end if
      if (active_Sl_soilw) then
        allocate(cam_in(c)%soilw(pcols), stat=ierror)
        if (ierror /= 0) call endrun(sub // ': allocation error soilw')
      end if
      if (active_Fall_flxdst1) then
        ! Assume 4 bins from surface model ....
        allocate(cam_in(c)%dstflx(pcols,4), stat=ierror)
        if (ierror /= 0) call endrun(sub // ': allocation error dstflx')
      end if
      if (active_Fall_flxvoc .and. shr_megan_mechcomps_n > 0) then
        allocate(cam_in(c)%meganflx(pcols,shr_megan_mechcomps_n), stat=ierror)
        if (ierror /= 0) call endrun(sub // ': allocation error meganflx')
      end if
    end do

    if (lnd_drydep .and. n_drydep > 0) then
      do c = begchunk, endchunk
        allocate (cam_in(c)%depvel(pcols,n_drydep), stat=ierror)
        if (ierror /= 0) call endrun(sub // ': allocation error depvel')
      end do
    end if

    if (active_Fall_flxfire .and. shr_fire_emis_mechcomps_n > 0) then
      do c = begchunk,endchunk
        allocate(cam_in(c)%fireflx(pcols,shr_fire_emis_mechcomps_n), stat=ierror)
        if (ierror /= 0) call endrun(sub // ': allocation error fireflx')
        allocate(cam_in(c)%fireztop(pcols), stat=ierror)
        if (ierror /= 0) call endrun(sub // ': allocation error fireztop')
      end do
    end if

    do c = begchunk, endchunk
      cam_in(c)%lchnk     = c
      cam_in(c)%ncol      = get_ncols_p(c)
      cam_in(c)%asdir     = 0
      cam_in(c)%asdif     = 0
      cam_in(c)%aldir     = 0
      cam_in(c)%aldif     = 0
      cam_in(c)%lwup      = 0
      cam_in(c)%lhf       = 0
      cam_in(c)%shf       = 0
      cam_in(c)%wsx       = 0
      cam_in(c)%wsy       = 0
      cam_in(c)%tref      = 0
      cam_in(c)%qref      = 0
      cam_in(c)%u10       = 0
      cam_in(c)%ts        = 0
      cam_in(c)%sst       = 0
      cam_in(c)%snowhland = 0
      cam_in(c)%snowhice  = 0
      cam_in(c)%fco2_lnd  = 0
      cam_in(c)%fco2_ocn  = 0
      cam_in(c)%fdms      = 0
      cam_in(c)%landfrac  = posinf
      cam_in(c)%icefrac   = posinf
      cam_in(c)%ocnfrac   = posinf

      if (associated(cam_in(c)%ram1    )) cam_in(c)%ram1     = 0.1_r8
      if (associated(cam_in(c)%fv      )) cam_in(c)%fv       = 0.1_r8
      if (associated(cam_in(c)%soilw   )) cam_in(c)%soilw    = 0
      if (associated(cam_in(c)%dstflx  )) cam_in(c)%dstflx   = 0
      if (associated(cam_in(c)%meganflx)) cam_in(c)%meganflx = 0

      cam_in(c)%cflx  = 0
      cam_in(c)%ustar = 0
      cam_in(c)%re    = 0
      cam_in(c)%ssq   = 0
      if (lnd_drydep .and. n_drydep > 0) then
        cam_in(c)%depvel = 0._r8
      end if
      if (active_Fall_flxfire .and. shr_fire_emis_mechcomps_n > 0) then
        cam_in(c)%fireflx  = 0
        cam_in(c)%fireztop = 0
      end if
    end do

  end subroutine hub2atm_alloc

  subroutine atm2hub_alloc(cam_out)

    ! Allocate space for the atmosphere to surface data type. And initialize the values.

    type(cam_out_t), pointer :: cam_out(:) ! Atmosphere to surface input

    integer c, ierror
    character(*), parameter :: sub = 'atm2hub_alloc'

    if (.not. phys_grid_initialized()) call endrun(sub // ': phys_grid not called yet')
    allocate(cam_out(begchunk:endchunk), stat=ierror)
    if (ierror /= 0) then
      write(iulog, *) sub // ': Allocation error: ', ierror
      call endrun(sub//': allocation error: cam_out')
    end if

    do c = begchunk,endchunk
      call cam_out(c)%init()
      cam_out(c)%lchnk    = c
      cam_out(c)%ncol     = get_ncols_p(c)
      cam_out(c)%tbot     = 0
      cam_out(c)%zbot     = 0
      cam_out(c)%topo     = 0
      cam_out(c)%ubot     = 0
      cam_out(c)%vbot     = 0
      cam_out(c)%qbot     = 0
      cam_out(c)%pbot     = 0
      cam_out(c)%rho      = 0
      cam_out(c)%netsw    = 0
      cam_out(c)%flwds    = 0
      cam_out(c)%precsc   = 0
      cam_out(c)%precsl   = 0
      cam_out(c)%precc    = 0
      cam_out(c)%precl    = 0
      cam_out(c)%soll     = 0
      cam_out(c)%sols     = 0
      cam_out(c)%solld    = 0
      cam_out(c)%solsd    = 0
      cam_out(c)%thbot    = 0
      cam_out(c)%co2prog  = 0
      cam_out(c)%co2diag  = 0
      cam_out(c)%psl      = 0
      cam_out(c)%bcphidry = 0
      cam_out(c)%bcphodry = 0
      cam_out(c)%bcphiwet = 0
      cam_out(c)%ocphidry = 0
      cam_out(c)%ocphodry = 0
      cam_out(c)%ocphiwet = 0
      cam_out(c)%dstdry1  = 0
      cam_out(c)%dstwet1  = 0
      cam_out(c)%dstdry2  = 0
      cam_out(c)%dstwet2  = 0
      cam_out(c)%dstdry3  = 0
      cam_out(c)%dstwet3  = 0
      cam_out(c)%dstdry4  = 0
      cam_out(c)%dstwet4  = 0

      if (active_Faxa_nhx) then
        allocate(cam_out(c)%nhx_nitrogen_flx(pcols), stat=ierror)
        if (ierror /= 0) call endrun(sub // ': allocation error nhx_nitrogen_flx')
        cam_out(c)%nhx_nitrogen_flx(:) = 0
      end if
      if (active_Faxa_noy) then
        allocate(cam_out(c)%noy_nitrogen_flx(pcols), stat=ierror)
        if (ierror /= 0) call endrun(sub // ': allocation error noy_nitrogen_flx')
        cam_out(c)%noy_nitrogen_flx(:) = 0
      end if
    end do

  end subroutine atm2hub_alloc

  subroutine atm2hub_deallocate(cam_out)

    type(cam_out_t), pointer :: cam_out(:) ! Atmosphere to surface input

    if (associated(cam_out)) then
      deallocate(cam_out)
    end if
    nullify(cam_out)

  end subroutine atm2hub_deallocate

  subroutine hub2atm_deallocate(cam_in)

    type(cam_in_t), pointer :: cam_in(:) ! Atmosphere to surface input

    integer c

    if (associated(cam_in)) then
      deallocate(cam_in)
    end if
    nullify(cam_in)

  end subroutine hub2atm_deallocate

  subroutine cam_export(state, cam_out, pbuf)

    ! Transfer atmospheric fields into necessary surface data structures

    use physics_types , only: physics_state
    use ppgrid        , only: pver
    use cam_history   , only: outfld
    use chem_surfvals , only: chem_surfvals_get
    use co2_cycle     , only: co2_transport, c_i
    use physconst     , only: rair, mwdry, mwco2, gravit
    use constituents  , only: pcnst
    use physics_buffer, only: pbuf_get_index, pbuf_get_field, physics_buffer_desc

    type(physics_state), intent(in) :: state
    type(cam_out_t), intent(inout) :: cam_out
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer i, m, lchnk, ncol
    integer psl_idx
    integer prec_dp_idx, snow_dp_idx, prec_sh_idx, snow_sh_idx
    integer prec_sed_idx,snow_sed_idx,prec_pcw_idx,snow_pcw_idx

    real(r8), pointer :: psl     (:)
    real(r8), pointer :: prec_dp (:)  ! Total precipitation from ZM convection
    real(r8), pointer :: snow_dp (:)  ! snow from ZM convection
    real(r8), pointer :: prec_sh (:)  ! total precipitation from Hack convection
    real(r8), pointer :: snow_sh (:)  ! snow from Hack convection
    real(r8), pointer :: prec_sed(:)  ! total precipitation from ZM convection
    real(r8), pointer :: snow_sed(:)  ! snow from ZM convection
    real(r8), pointer :: prec_pcw(:)  ! total precipitation from Hack convection
    real(r8), pointer :: snow_pcw(:)  ! snow from Hack convection

    lchnk = state%lchnk
    ncol  = state%ncol

    psl_idx = pbuf_get_index('PSL')
    call pbuf_get_field(pbuf, psl_idx, psl)

    prec_dp_idx  = pbuf_get_index('PREC_DP' , errcode=i)
    snow_dp_idx  = pbuf_get_index('SNOW_DP' , errcode=i)
    prec_sh_idx  = pbuf_get_index('PREC_SH' , errcode=i)
    snow_sh_idx  = pbuf_get_index('SNOW_SH' , errcode=i)
    prec_sed_idx = pbuf_get_index('PREC_SED', errcode=i)
    snow_sed_idx = pbuf_get_index('SNOW_SED', errcode=i)
    prec_pcw_idx = pbuf_get_index('PREC_PCW', errcode=i)
    snow_pcw_idx = pbuf_get_index('SNOW_PCW', errcode=i)

    if (prec_dp_idx  > 0) call pbuf_get_field(pbuf, prec_dp_idx , prec_dp )
    if (snow_dp_idx  > 0) call pbuf_get_field(pbuf, snow_dp_idx , snow_dp )
    if (prec_sh_idx  > 0) call pbuf_get_field(pbuf, prec_sh_idx , prec_sh )
    if (snow_sh_idx  > 0) call pbuf_get_field(pbuf, snow_sh_idx , snow_sh )
    if (prec_sed_idx > 0) call pbuf_get_field(pbuf, prec_sed_idx, prec_sed)
    if (snow_sed_idx > 0) call pbuf_get_field(pbuf, snow_sed_idx, snow_sed)
    if (prec_pcw_idx > 0) call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw)
    if (snow_pcw_idx > 0) call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw)

    do i = 1, ncol
      cam_out%tbot (i) = state%t   (i,pver)
      cam_out%thbot(i) = state%t   (i,pver) * state%exner(i,pver)
      cam_out%zbot (i) = state%zm  (i,pver)
      cam_out%topo (i) = state%phis(i) / gravit
      cam_out%ubot (i) = state%u   (i,pver)
      cam_out%vbot (i) = state%v   (i,pver)
      cam_out%pbot (i) = state%pmid(i,pver)
      cam_out%psl  (i) = psl(i)
      cam_out%rho  (i) = cam_out%pbot(i) / (rair * cam_out%tbot(i))
    end do
    do m = 1, pcnst
      do i = 1, ncol
        cam_out%qbot(i,m) = state%q(i,pver,m)
      end do
    end do

    cam_out%co2diag(:ncol) = chem_surfvals_get('CO2VMR') * 1.0e+6_r8
    if (co2_transport()) then
      do i = 1, ncol
        cam_out%co2prog(i) = state%q(i,pver,c_i(4)) * 1.0e+6_r8 * mwdry / mwco2
      end do
    end if
    !
    ! Precipitation and snow rates from shallow convection, deep convection and stratiform processes.
    ! Compute total convective and stratiform precipitation and snow rates
    !
    do i = 1, ncol
      cam_out%precc (i) = 0
      cam_out%precl (i) = 0
      cam_out%precsc(i) = 0
      cam_out%precsl(i) = 0
      if (prec_dp_idx  > 0) cam_out%precc (i) = cam_out%precc (i) + prec_dp (i)
      if (prec_sh_idx  > 0) cam_out%precc (i) = cam_out%precc (i) + prec_sh (i)
      if (prec_sed_idx > 0) cam_out%precl (i) = cam_out%precl (i) + prec_sed(i)
      if (prec_pcw_idx > 0) cam_out%precl (i) = cam_out%precl (i) + prec_pcw(i)
      if (snow_dp_idx  > 0) cam_out%precsc(i) = cam_out%precsc(i) + snow_dp (i)
      if (snow_sh_idx  > 0) cam_out%precsc(i) = cam_out%precsc(i) + snow_sh (i)
      if (snow_sed_idx > 0) cam_out%precsl(i) = cam_out%precsl(i) + snow_sed(i)
      if (snow_pcw_idx > 0) cam_out%precsl(i) = cam_out%precsl(i) + snow_pcw(i)

      ! NOTE: These checks should not be necessary if they exist in the parameterizations
      if (cam_out%precc (i) < 0) cam_out%precc (i) = 0
      if (cam_out%precl (i) < 0) cam_out%precl (i) = 0
      if (cam_out%precsc(i) < 0) cam_out%precsc(i) = 0
      if (cam_out%precsl(i) < 0) cam_out%precsl(i) = 0
      if (cam_out%precsc(i) > cam_out%precc(i)) cam_out%precsc(i) = cam_out%precc(i)
      if (cam_out%precsl(i) > cam_out%precl(i)) cam_out%precsl(i) = cam_out%precl(i)
    end do

  end subroutine cam_export

end module camsrfexch
