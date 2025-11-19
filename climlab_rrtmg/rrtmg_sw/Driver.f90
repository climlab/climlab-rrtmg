!  CLIMLAB driver for RRTMG_SW radiation
!
!   using the latest version RRTMG_SW_v4.0
!
!  This is a lightweight driver designed for wrapping with f2py
!  and calling from Python
!  Refer to RRTMG_SW source code for more documentation
!
!  Three functions are exposed here:
!   climlab_rrtmg_sw_ini (wrapper for rrtmg_sw_ini)
!   climlab_mcica_subcol_sw (wrapper for mcica_subcol_sw)
!   climlab_rrtmg_sw  (wrapper for rrtmg_sw)
!
!   The call signature for each of these is nearly identical to its
!   equivalent in the RRTMG_SW code.
!
!   See the python module climlab/radiation/rrtm/rrtmg_sw.py
!    to see how these are called from Python
!
!  Brian Rose
!  brose@albany.edu
!

subroutine climlab_rrtmg_sw_ini(cpdair)
    ! Modules
    use parkind, only : rb => kind_rb
    use rrtmg_sw_init, only: rrtmg_sw_ini
    ! Input
    real(kind=rb), intent(in) :: cpdair    ! Specific heat capacity of dry air
                                            ! at constant pressure at 273 K
                                            ! (J kg-1 K-1)
    ! Call the initialization routine
    call rrtmg_sw_ini(cpdair)
end subroutine climlab_rrtmg_sw_ini


subroutine climlab_mcica_subcol_sw &
  (ncol, nlay, icld, permuteseed, irng, play, &
   cldfrac, ciwp, clwp, reic, relq, tauc, ssac, asmc, fsfc, &
   cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, &
   taucmcl, ssacmcl, asmcmcl, fsfcmcl)

! Modules
  use parkind, only : im => kind_im
  use mcica_subcol_gen_sw, only: mcica_subcol_sw
  !use parrrsw, only: nbndsw, ngptsw
  integer(kind=im), parameter :: nbndsw = 14
  integer(kind=im), parameter :: ngptsw = 112

! Input
  integer, parameter :: rb = selected_real_kind(12)
  integer(kind=im), intent(in) :: ncol            ! number of columns
  integer(kind=im), intent(in) :: nlay            ! number of model layers
  integer(kind=im), intent(inout) :: icld         ! Cloud overlap method
                                                  !    0: Clear only
                                                  !    1: Random
                                                  !    2: Maximum/random
                                                  !    3: Maximum
  integer(kind=im), intent(in) :: permuteseed     ! if the cloud generator is called multiple times,
                                                  ! permute the seed between each call.
                                                  ! between calls for LW and SW, recommended
                                                  ! permuteseed differes by 'ngpt'
  integer(kind=im), intent(inout) :: irng         ! flag for random number generator
                                                  !  0 = kissvec
                                                  !  1 = Mersenne Twister
  real(kind=rb), intent(in) :: play(ncol,nlay)    ! Layer pressures (hPa, mb)
  real(kind=rb), intent(in) :: cldfrac(ncol,nlay)        ! layer cloud fraction
  real(kind=rb), intent(in) :: ciwp(ncol,nlay)           ! in-cloud ice water path
  real(kind=rb), intent(in) :: clwp(ncol,nlay)           ! in-cloud liquid water path
  real(kind=rb), intent(in) :: reic(ncol,nlay)           ! cloud ice particle size
  real(kind=rb), intent(in) :: relq(ncol,nlay)           ! cloud liquid particle size
  real(kind=rb), intent(in) :: tauc(nbndsw,ncol,nlay)    ! in-cloud optical depth
  real(kind=rb), intent(in) :: ssac(nbndsw,ncol,nlay)    ! in-cloud single scattering albedo (non-delta scaled)
  real(kind=rb), intent(in) :: asmc(nbndsw,ncol,nlay)    ! in-cloud asymmetry parameter (non-delta scaled)
  real(kind=rb), intent(in) :: fsfc(nbndsw,ncol,nlay)    ! in-cloud forward scattering fraction (non-delta scaled)
  ! Output
  !   These quantities are computed by McICA
  real(kind=rb), intent(out) :: cldfmcl(ngptsw,ncol,nlay)   ! cloud fraction [mcica]
  real(kind=rb), intent(out) :: ciwpmcl(ngptsw,ncol,nlay)   ! in-cloud ice water path [mcica]
  real(kind=rb), intent(out) :: clwpmcl(ngptsw,ncol,nlay)   ! in-cloud liquid water path [mcica]
  real(kind=rb), intent(out) :: reicmcl(ncol,nlay)          ! ice partcle size (microns)
  real(kind=rb), intent(out) :: relqmcl(ncol,nlay)          ! liquid particle size (microns)
  real(kind=rb), intent(out) :: taucmcl(ngptsw,ncol,nlay)   ! in-cloud optical depth [mcica]
  real(kind=rb), intent(out) :: ssacmcl(ngptsw,ncol,nlay)   ! in-cloud single scattering albedo [mcica]
  real(kind=rb), intent(out) :: asmcmcl(ngptsw,ncol,nlay)   ! in-cloud asymmetry parameter [mcica]
  real(kind=rb), intent(out) :: fsfcmcl(ngptsw,ncol,nlay)   ! in-cloud forward scattering fraction [mcica]
!  These are not comments! Necessary directives to f2py to handle array dimensions
!f2py depend(ncol,nlay) play,
!f2py depend(ncol,nlay) cldfrac,ciwp,clwp,reic,relq,tauc
!f2py depend(ncol,nlay) tauc,ssac,asmc,fsfc
!f2py depend(ncol,nlay) reicmcl,relqmcl
!f2py depend(ncol,nlay) cldfmcl,ciwpmcl,clwpmcl,taucmcl,ssacmcl,asmcmcl,fsfcmcl

  ! Call the Monte Carlo Independent Column Approximation
  !   (McICA, Pincus et al., JC, 2003)
  call mcica_subcol_sw(1, ncol, nlay, icld, permuteseed, irng, play, &
                     cldfrac, ciwp, clwp, reic, relq, tauc, ssac, asmc, fsfc, &
                     cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, &
                     taucmcl, ssacmcl, asmcmcl, fsfcmcl)

end subroutine climlab_mcica_subcol_sw


subroutine climlab_rrtmg_sw_expanded &
    (ncol    ,nlay    ,icld     ,ispec    ,iaer    , &
    play    ,plev    ,tlay    ,tlev    ,tsfc   , &
    h2ovmr , o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr ,o2vmr , &
    asdir   ,asdif   ,aldir   ,aldif   , &
    kmodts, coszen  ,adjes   ,dyofyr  ,scon    ,isolvar, &
    inflgsw ,iceflgsw,liqflgsw,cldfmcl , &
    taucmcl ,ssacmcl ,asmcmcl ,fsfcmcl , &
    ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
    tauaer  ,ssaaer  ,asmaer  ,ecaer   , &
    bndsolvar,indsolvar,solcycfrac, &
    swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
    swuflxspec, swdflxspec, swuflxcspec, swdflxcspec, &
    add_aero_layer, r_mu, t_mu, r_bar, t_bar)

! Modules
    use parkind, only : im => kind_im
    use rrtmg_sw_rad, only: rrtmg_sw
    ! use parrrsw, only: nbndsw, ngptsw, naerec
    integer(kind=im), parameter :: nbndsw = 14
    integer(kind=im), parameter :: naerec  = 6
    integer(kind=im), parameter :: ngptsw = 112
! Input
    integer, parameter :: rb = selected_real_kind(12)
    integer(kind=im), intent(in) :: ncol            ! number of columns
    integer(kind=im), intent(in) :: nlay            ! number of model layers
    integer(kind=im), intent(inout) :: icld         ! Cloud overlap method
                                                    !    0: Clear only
                                                    !    1: Random
                                                    !    2: Maximum/random
                                                    !    3: Maximum
    integer(kind=im), intent(inout) :: ispec        ! spectral ASR output flag
    integer(kind=im), intent(inout) :: iaer         ! Aerosol option flag
                                                    !    0: No aerosol
                                                    !    6: ECMWF method
                                                    !    10:Input aerosol optical
                                                    !       properties
    real(kind=rb), intent(in) :: play(ncol,nlay)    ! Layer pressures (hPa, mb)
    real(kind=rb), intent(in) :: plev(ncol,nlay+1)  ! Interface pressures (hPa, mb)
    real(kind=rb), intent(in) :: tlay(ncol,nlay)    ! Layer temperatures (K)
    real(kind=rb), intent(in) :: tlev(ncol,nlay+1)  ! Interface temperatures (K)
    real(kind=rb), intent(in) :: tsfc(ncol)         ! Surface temperature (K)
    real(kind=rb), intent(in) :: h2ovmr(ncol,nlay)  ! H2O volume mixing ratio
    real(kind=rb), intent(in) :: o3vmr(ncol,nlay)   ! O3 volume mixing ratio
    real(kind=rb), intent(in) :: co2vmr(ncol,nlay)  ! CO2 volume mixing ratio
    real(kind=rb), intent(in) :: ch4vmr(ncol,nlay)  ! Methane volume mixing ratio
    real(kind=rb), intent(in) :: n2ovmr(ncol,nlay)  ! Nitrous oxide volume mixing ratio
    real(kind=rb), intent(in) :: o2vmr(ncol,nlay)   ! Oxygen volume mixing ratio
    real(kind=rb), intent(in) :: aldif(ncol)        ! UV/vis surface albedo direct rad
    real(kind=rb), intent(in) :: aldir(ncol)        ! Near-IR surface albedo direct rad
    real(kind=rb), intent(in) :: asdif(ncol)        ! UV/vis surface albedo: diffuse rad
    real(kind=rb), intent(in) :: asdir(ncol)        ! Near-IR surface albedo: diffuse rad
    integer(kind=im), intent(in) :: kmodts
    real(kind=rb), intent(in) :: coszen(ncol)       ! Cosine of solar zenith angle
    real(kind=rb), intent(in) :: adjes(ncol)        ! Flux adjustment (Earth/Sun distance and/or zenith angle compensation)
    integer(kind=im), intent(in) :: dyofyr          ! Day of the year (used to get Earth/Sun
                                                    !  distance if adjflx not provided)
    real(kind=rb), intent(in) :: scon               ! Solar constant (W/m2)
                                                    !    Total solar irradiance averaged
                                                    !    over the solar cycle.
                                                    !    If scon = 0.0, the internal solar
                                                    !    constant, which depends on the
                                                    !    value of isolvar, will be used.
                                                    !    For isolvar=-1, scon=1368.22 Wm-2,
                                                    !    For isolvar=0,1,3, scon=1360.85 Wm-2,
                                                    !    If scon > 0.0, the internal solar
                                                    !    constant will be scaled to the
                                                    !    provided value of scon.
    integer(kind=im), intent(in) :: isolvar         ! Flag for solar variability method
                                                    !   -1 = (when scon .eq. 0.0): No solar variability
                                                    !        and no solar cycle (Kurucz solar irradiance
                                                    !        of 1368.22 Wm-2 only);
                                                    !        (when scon .ne. 0.0): Kurucz solar irradiance
                                                    !        scaled to scon and solar variability defined
                                                    !        (optional) by setting non-zero scale factors
                                                    !        for each band in bndsolvar
                                                    !    0 = (when SCON .eq. 0.0): No solar variability
                                                    !        and no solar cycle (NRLSSI2 solar constant of
                                                    !        1360.85 Wm-2 for the 100-50000 cm-1 spectral
                                                    !        range only), with facular and sunspot effects
                                                    !        fixed to the mean of Solar Cycles 13-24;
                                                    !        (when SCON .ne. 0.0): No solar variability
                                                    !        and no solar cycle (NRLSSI2 solar constant of
                                                    !        1360.85 Wm-2 for the 100-50000 cm-1 spectral
                                                    !        range only), is scaled to SCON
                                                    !    1 = Solar variability (using NRLSSI2  solar
                                                    !        model) with solar cycle contribution
                                                    !        determined by fraction of solar cycle
                                                    !        with facular and sunspot variations
                                                    !        fixed to their mean variations over the
                                                    !        average of Solar Cycles 13-24;
                                                    !        two amplitude scale factors allow
                                                    !        facular and sunspot adjustments from
                                                    !        mean solar cycle as defined by indsolvar
                                                    !    2 = Solar variability (using NRLSSI2 solar
                                                    !        model) over solar cycle determined by
                                                    !        direct specification of Mg (facular)
                                                    !        and SB (sunspot) indices provided
                                                    !        in indsolvar (scon = 0.0 only)
                                                    !    3 = (when scon .eq. 0.0): No solar variability
                                                    !        and no solar cycle (NRLSSI2 solar irradiance
                                                    !        of 1360.85 Wm-2 only);
                                                    !        (when scon .ne. 0.0): NRLSSI2 solar irradiance
                                                    !        scaled to scon and solar variability defined
                                                    !        (optional) by setting non-zero scale factors
                                                    !        for each band in bndsolvar
    real(kind=rb), intent(inout) :: indsolvar(2) ! Facular and sunspot amplitude
                                              ! scale factors (isolvar=1), or
                                              ! Mg and SB indices (isolvar=2)
                                              !    Dimensions: (2)
    real(kind=rb), intent(inout) :: bndsolvar(nbndsw) ! Solar variability scale factors
                                                   ! for each shortwave band
                                                   !    Dimensions: (nbndsw=14)
    real(kind=rb), intent(inout) :: solcycfrac   ! Fraction of averaged solar cycle (0-1)
                                              !    at current time (isolvar=1)

    integer(kind=im), intent(in) :: inflgsw         ! Flag for cloud optical properties
    integer(kind=im), intent(in) :: iceflgsw        ! Flag for ice particle specification
    integer(kind=im), intent(in) :: liqflgsw        ! Flag for liquid droplet specification
    !   These quantities are computed by McICA
    real(kind=rb), intent(in) :: cldfmcl(ngptsw,ncol,nlay)   ! cloud fraction [mcica]
    real(kind=rb), intent(in) :: ciwpmcl(ngptsw,ncol,nlay)   ! in-cloud ice water path [mcica]
    real(kind=rb), intent(in) :: clwpmcl(ngptsw,ncol,nlay)   ! in-cloud liquid water path [mcica]
    real(kind=rb), intent(in) :: reicmcl(ncol,nlay)          ! ice partcle size (microns)
    real(kind=rb), intent(in) :: relqmcl(ncol,nlay)          ! liquid particle size (microns)
    real(kind=rb), intent(in) :: taucmcl(ngptsw,ncol,nlay)   ! in-cloud optical depth [mcica]
    real(kind=rb), intent(in) :: ssacmcl(ngptsw,ncol,nlay)   ! in-cloud single scattering albedo [mcica]
    real(kind=rb), intent(in) :: asmcmcl(ngptsw,ncol,nlay)   ! in-cloud asymmetry parameter [mcica]
    real(kind=rb), intent(in) :: fsfcmcl(ngptsw,ncol,nlay)   ! in-cloud forward scattering fraction [mcica]
    real(kind=rb), intent(in) :: tauaer(ncol,nlay,nbndsw)  ! Aerosol optical depth (iaer=10 only)
    real(kind=rb), intent(in) :: ssaaer(ncol,nlay,nbndsw)  ! Aerosol single scattering albedo (iaer=10 only)
    real(kind=rb), intent(in) :: asmaer(ncol,nlay,nbndsw)  ! Aerosol asymmetry parameter (iaer=10 only)
    real(kind=rb), intent(in) :: ecaer(ncol,nlay,naerec)   ! Aerosol optical depth at 0.55 micron (iaer=6 only)

! Output
    real(kind=rb), intent(out) :: swuflx(ncol,nlay+1)    ! Total sky shortwave upward flux (W/m2)
    real(kind=rb), intent(out) :: swdflx(ncol,nlay+1)    ! Total sky shortwave downward flux (W/m2)
    real(kind=rb), intent(out) :: swhr(ncol,nlay)        ! Total sky shortwave radiative heating rate (K/d)
    real(kind=rb), intent(out) :: swuflxc(ncol,nlay+1)   ! Clear sky shortwave upward flux (W/m2)
    real(kind=rb), intent(out) :: swdflxc(ncol,nlay+1)   ! Clear sky shortwave downward flux (W/m2)
    real(kind=rb), intent(out) :: swhrc(ncol,nlay)       ! Clear sky shortwave radiative heating rate (K/d)
    real(kind=rb), intent(out) :: swuflxspec(ncol,nlay+1,nbndsw)    ! Total sky shortwave upward flux spectrum (W/m2)
    real(kind=rb), intent(out) :: swdflxspec(ncol,nlay+1,nbndsw)    ! Total sky shortwave downward flux spectrum (W/m2)
    real(kind=rb), intent(out) :: swuflxcspec(ncol,nlay+1,nbndsw)   ! Clear sky shortwave upward flux spectrum (W/m2)
    real(kind=rb), intent(out) :: swdflxcspec(ncol,nlay+1,nbndsw)   ! Clear sky shortwave downward flux spectrum (W/m2)
    integer(kind=im), intent(in) :: add_aero_layer
    real(kind=rb), intent(in) :: r_mu(ncol,nlay,nbndsw)    ! Aerosols directional reflection
    real(kind=rb), intent(in) :: t_mu(ncol,nlay,nbndsw)    ! Aerosols directional transmission
    real(kind=rb), intent(in) :: r_bar(ncol,nlay,nbndsw)    ! Aerosols diffusive reflection
    real(kind=rb), intent(in) :: t_bar(ncol,nlay,nbndsw)    ! Aerosols diffusive transmission

!  These are not comments! Necessary directives to f2py to handle array dimensions
!f2py depend(ncol,nlay) play
!f2py depend(ncol,nlay) plev
!f2py depend(ncol,nlay) tlay
!f2py depend(ncol,nlay) tlev
!f2py depend(ncol,nlay) h2ovmr,o3vmr,co2vmr,ch4vmr,n2ovmr,o2vmr
!f2py depend(ncol) asdir,asdif,aldir,aldif,coszen, adjes, tsfc
!f2py depend(ncol,nlay) tauaer,ssaaer,asmaer,ecaer
!f2py depend(ncol,nlay) cldfmcl,taucmcl,ssacmcl,asmcmcl,fsfcmcl
!f2py depend(ncol,nlay) ciwpmcl,clwpmcl
!f2py depend(ncol,nlay) reicmcl,relqmcl
!f2py depend(ncol,nlay) swuflx,swdflx
!f2py depend(ncol,nlay) swhr
!f2py depend(ncol,nlay) swuflxc,swdflxc
!f2py depend(ncol,nlay) swhrc
!f2py depend(ncol,nlay,nbndsw) swuflxspec,swdflxspec,swuflxcspec,swdflxcspec
!f2py depend(ncol,nlay,nbndsw) r_mu, t_mu, r_bar, t_bar

    !  Call the RRTMG_SW driver to compute radiative fluxes
    call rrtmg_sw(ncol    ,nlay    ,icld    ,ispec    ,iaer    , &
              play    ,plev    ,tlay    ,tlev    ,tsfc   , &
              h2ovmr , o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr ,o2vmr , &
              asdir   ,asdif   ,aldir   ,aldif   , &
              kmodts, coszen  ,adjes   ,dyofyr  ,scon    ,isolvar, &
              inflgsw ,iceflgsw,liqflgsw,cldfmcl , &
              taucmcl ,ssacmcl ,asmcmcl ,fsfcmcl , &
              ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
              tauaer  ,ssaaer  ,asmaer  ,ecaer   , &
              swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc, &
              swuflxspec, swdflxspec, swuflxcspec, swdflxcspec, &
              bndsolvar,indsolvar,solcycfrac, &
              add_aero_layer, r_mu, t_mu, r_bar, t_bar)

end subroutine climlab_rrtmg_sw_expanded

subroutine climlab_rrtmg_sw &
    (ncol    ,nlay    ,icld    ,iaer    , &
    play    ,plev    ,tlay    ,tlev    ,tsfc   , &
    h2ovmr , o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr ,o2vmr , &
    asdir   ,asdif   ,aldir   ,aldif   , &
    coszen  ,adjes   ,dyofyr  ,scon    ,isolvar, &
    inflgsw ,iceflgsw,liqflgsw,cldfmcl , &
    taucmcl ,ssacmcl ,asmcmcl ,fsfcmcl , &
    ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
    tauaer  ,ssaaer  ,asmaer  ,ecaer   , &
    bndsolvar,indsolvar,solcycfrac, &
    swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc)

! Modules
    use parkind, only : im => kind_im
    use rrtmg_sw_rad, only: rrtmg_sw
    ! use parrrsw, only: nbndsw, ngptsw, naerec
    integer(kind=im), parameter :: nbndsw = 14
    integer(kind=im), parameter :: naerec  = 6
    integer(kind=im), parameter :: ngptsw = 112
! Input
    integer, parameter :: rb = selected_real_kind(12)
    integer(kind=im), intent(in) :: ncol            ! number of columns
    integer(kind=im), intent(in) :: nlay            ! number of model layers
    integer(kind=im), intent(inout) :: icld         ! Cloud overlap method
                                                    !    0: Clear only
                                                    !    1: Random
                                                    !    2: Maximum/random
                                                    !    3: Maximum
    integer(kind=im), intent(inout) :: iaer         ! Aerosol option flag
                                                    !    0: No aerosol
                                                    !    6: ECMWF method
                                                    !    10:Input aerosol optical
                                                    !       properties
    real(kind=rb), intent(in) :: play(ncol,nlay)    ! Layer pressures (hPa, mb)
    real(kind=rb), intent(in) :: plev(ncol,nlay+1)  ! Interface pressures (hPa, mb)
    real(kind=rb), intent(in) :: tlay(ncol,nlay)    ! Layer temperatures (K)
    real(kind=rb), intent(in) :: tlev(ncol,nlay+1)  ! Interface temperatures (K)
    real(kind=rb), intent(in) :: tsfc(ncol)         ! Surface temperature (K)
    real(kind=rb), intent(in) :: h2ovmr(ncol,nlay)  ! H2O volume mixing ratio
    real(kind=rb), intent(in) :: o3vmr(ncol,nlay)   ! O3 volume mixing ratio
    real(kind=rb), intent(in) :: co2vmr(ncol,nlay)  ! CO2 volume mixing ratio
    real(kind=rb), intent(in) :: ch4vmr(ncol,nlay)  ! Methane volume mixing ratio
    real(kind=rb), intent(in) :: n2ovmr(ncol,nlay)  ! Nitrous oxide volume mixing ratio
    real(kind=rb), intent(in) :: o2vmr(ncol,nlay)   ! Oxygen volume mixing ratio
    real(kind=rb), intent(in) :: aldif(ncol)        ! UV/vis surface albedo direct rad
    real(kind=rb), intent(in) :: aldir(ncol)        ! Near-IR surface albedo direct rad
    real(kind=rb), intent(in) :: asdif(ncol)        ! UV/vis surface albedo: diffuse rad
    real(kind=rb), intent(in) :: asdir(ncol)        ! Near-IR surface albedo: diffuse rad
    real(kind=rb), intent(in) :: coszen(ncol)       ! Cosine of solar zenith angle
    real(kind=rb), intent(in) :: adjes(ncol)        ! Flux adjustment (Earth/Sun distance and/or zenith angle compensation)
    integer(kind=im), intent(in) :: dyofyr          ! Day of the year (used to get Earth/Sun
                                                    !  distance if adjflx not provided)
    real(kind=rb), intent(in) :: scon               ! Solar constant (W/m2)
                                                    !    Total solar irradiance averaged
                                                    !    over the solar cycle.
                                                    !    If scon = 0.0, the internal solar
                                                    !    constant, which depends on the
                                                    !    value of isolvar, will be used.
                                                    !    For isolvar=-1, scon=1368.22 Wm-2,
                                                    !    For isolvar=0,1,3, scon=1360.85 Wm-2,
                                                    !    If scon > 0.0, the internal solar
                                                    !    constant will be scaled to the
                                                    !    provided value of scon.
    integer(kind=im), intent(in) :: isolvar         ! Flag for solar variability method
                                                    !   -1 = (when scon .eq. 0.0): No solar variability
                                                    !        and no solar cycle (Kurucz solar irradiance
                                                    !        of 1368.22 Wm-2 only);
                                                    !        (when scon .ne. 0.0): Kurucz solar irradiance
                                                    !        scaled to scon and solar variability defined
                                                    !        (optional) by setting non-zero scale factors
                                                    !        for each band in bndsolvar
                                                    !    0 = (when SCON .eq. 0.0): No solar variability
                                                    !        and no solar cycle (NRLSSI2 solar constant of
                                                    !        1360.85 Wm-2 for the 100-50000 cm-1 spectral
                                                    !        range only), with facular and sunspot effects
                                                    !        fixed to the mean of Solar Cycles 13-24;
                                                    !        (when SCON .ne. 0.0): No solar variability
                                                    !        and no solar cycle (NRLSSI2 solar constant of
                                                    !        1360.85 Wm-2 for the 100-50000 cm-1 spectral
                                                    !        range only), is scaled to SCON
                                                    !    1 = Solar variability (using NRLSSI2  solar
                                                    !        model) with solar cycle contribution
                                                    !        determined by fraction of solar cycle
                                                    !        with facular and sunspot variations
                                                    !        fixed to their mean variations over the
                                                    !        average of Solar Cycles 13-24;
                                                    !        two amplitude scale factors allow
                                                    !        facular and sunspot adjustments from
                                                    !        mean solar cycle as defined by indsolvar
                                                    !    2 = Solar variability (using NRLSSI2 solar
                                                    !        model) over solar cycle determined by
                                                    !        direct specification of Mg (facular)
                                                    !        and SB (sunspot) indices provided
                                                    !        in indsolvar (scon = 0.0 only)
                                                    !    3 = (when scon .eq. 0.0): No solar variability
                                                    !        and no solar cycle (NRLSSI2 solar irradiance
                                                    !        of 1360.85 Wm-2 only);
                                                    !        (when scon .ne. 0.0): NRLSSI2 solar irradiance
                                                    !        scaled to scon and solar variability defined
                                                    !        (optional) by setting non-zero scale factors
                                                    !        for each band in bndsolvar
    real(kind=rb), intent(inout) :: indsolvar(2) ! Facular and sunspot amplitude
                                              ! scale factors (isolvar=1), or
                                              ! Mg and SB indices (isolvar=2)
                                              !    Dimensions: (2)
    real(kind=rb), intent(inout) :: bndsolvar(nbndsw) ! Solar variability scale factors
                                                   ! for each shortwave band
                                                   !    Dimensions: (nbndsw=14)
    real(kind=rb), intent(inout) :: solcycfrac   ! Fraction of averaged solar cycle (0-1)
                                              !    at current time (isolvar=1)

    integer(kind=im), intent(in) :: inflgsw         ! Flag for cloud optical properties
    integer(kind=im), intent(in) :: iceflgsw        ! Flag for ice particle specification
    integer(kind=im), intent(in) :: liqflgsw        ! Flag for liquid droplet specification
    !   These quantities are computed by McICA
    real(kind=rb), intent(in) :: cldfmcl(ngptsw,ncol,nlay)   ! cloud fraction [mcica]
    real(kind=rb), intent(in) :: ciwpmcl(ngptsw,ncol,nlay)   ! in-cloud ice water path [mcica]
    real(kind=rb), intent(in) :: clwpmcl(ngptsw,ncol,nlay)   ! in-cloud liquid water path [mcica]
    real(kind=rb), intent(in) :: reicmcl(ncol,nlay)          ! ice partcle size (microns)
    real(kind=rb), intent(in) :: relqmcl(ncol,nlay)          ! liquid particle size (microns)
    real(kind=rb), intent(in) :: taucmcl(ngptsw,ncol,nlay)   ! in-cloud optical depth [mcica]
    real(kind=rb), intent(in) :: ssacmcl(ngptsw,ncol,nlay)   ! in-cloud single scattering albedo [mcica]
    real(kind=rb), intent(in) :: asmcmcl(ngptsw,ncol,nlay)   ! in-cloud asymmetry parameter [mcica]
    real(kind=rb), intent(in) :: fsfcmcl(ngptsw,ncol,nlay)   ! in-cloud forward scattering fraction [mcica]
    real(kind=rb), intent(in) :: tauaer(ncol,nlay,nbndsw)  ! Aerosol optical depth (iaer=10 only)
    real(kind=rb), intent(in) :: ssaaer(ncol,nlay,nbndsw)  ! Aerosol single scattering albedo (iaer=10 only)
    real(kind=rb), intent(in) :: asmaer(ncol,nlay,nbndsw)  ! Aerosol asymmetry parameter (iaer=10 only)
    real(kind=rb), intent(in) :: ecaer(ncol,nlay,naerec)   ! Aerosol optical depth at 0.55 micron (iaer=6 only)

! Output
    real(kind=rb), intent(out) :: swuflx(ncol,nlay+1)    ! Total sky shortwave upward flux (W/m2)
    real(kind=rb), intent(out) :: swdflx(ncol,nlay+1)    ! Total sky shortwave downward flux (W/m2)
    real(kind=rb), intent(out) :: swhr(ncol,nlay)        ! Total sky shortwave radiative heating rate (K/d)
    real(kind=rb), intent(out) :: swuflxc(ncol,nlay+1)   ! Clear sky shortwave upward flux (W/m2)
    real(kind=rb), intent(out) :: swdflxc(ncol,nlay+1)   ! Clear sky shortwave downward flux (W/m2)
    real(kind=rb), intent(out) :: swhrc(ncol,nlay)       ! Clear sky shortwave radiative heating rate (K/d)

    integer(kind=im) :: ispec = 0        ! spectral ASR output flag
    integer(kind=im) :: kmodts = 2
    real(kind=rb) :: swuflxspec(ncol,nlay+1,nbndsw)    ! Total sky shortwave upward flux spectrum (W/m2)
    real(kind=rb) :: swdflxspec(ncol,nlay+1,nbndsw)    ! Total sky shortwave downward flux spectrum (W/m2)
    real(kind=rb) :: swuflxcspec(ncol,nlay+1,nbndsw)   ! Clear sky shortwave upward flux spectrum (W/m2)
    real(kind=rb) :: swdflxcspec(ncol,nlay+1,nbndsw)   ! Clear sky shortwave downward flux spectrum (W/m2)
    integer(kind=im) :: add_aero_layer = 0
    real(kind=rb) :: r_mu(ncol,nlay,nbndsw)    ! Aerosols directional reflection
    real(kind=rb) :: t_mu(ncol,nlay,nbndsw)    ! Aerosols directional transmission
    real(kind=rb) :: r_bar(ncol,nlay,nbndsw)    ! Aerosols diffusive reflection
    real(kind=rb) :: t_bar(ncol,nlay,nbndsw)    ! Aerosols diffusive transmission

!  These are not comments! Necessary directives to f2py to handle array dimensions
!f2py depend(ncol,nlay) play
!f2py depend(ncol,nlay) plev
!f2py depend(ncol,nlay) tlay
!f2py depend(ncol,nlay) tlev
!f2py depend(ncol,nlay) h2ovmr,o3vmr,co2vmr,ch4vmr,n2ovmr,o2vmr
!f2py depend(ncol) asdir,asdif,aldir,aldif,coszen, adjes, tsfc
!f2py depend(ncol,nlay) tauaer,ssaaer,asmaer,ecaer
!f2py depend(ncol,nlay) cldfmcl,taucmcl,ssacmcl,asmcmcl,fsfcmcl
!f2py depend(ncol,nlay) ciwpmcl,clwpmcl
!f2py depend(ncol,nlay) reicmcl,relqmcl
!f2py depend(ncol,nlay) swuflx,swdflx
!f2py depend(ncol,nlay) swhr
!f2py depend(ncol,nlay) swuflxc,swdflxc
!f2py depend(ncol,nlay) swhrc

    r_mu = 0.0
    t_mu = 1.0
    r_bar = 0.0
    t_bar = 1.0                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
    !  Call the RRTMG_SW driver to compute radiative fluxes
    call rrtmg_sw(ncol    ,nlay    ,icld    ,ispec    ,iaer    , &
              play    ,plev    ,tlay    ,tlev    ,tsfc   , &
              h2ovmr , o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr ,o2vmr , &
              asdir   ,asdif   ,aldir   ,aldif   , &
              kmodts, coszen  ,adjes   ,dyofyr  ,scon    ,isolvar, &
              inflgsw ,iceflgsw,liqflgsw,cldfmcl , &
              taucmcl ,ssacmcl ,asmcmcl ,fsfcmcl , &
              ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
              tauaer  ,ssaaer  ,asmaer  ,ecaer   , &
              swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc, &
              swuflxspec, swdflxspec, swuflxcspec, swdflxcspec, &
              bndsolvar,indsolvar,solcycfrac, &
              add_aero_layer, r_mu, t_mu, r_bar, t_bar)

end subroutine climlab_rrtmg_sw

subroutine climlab_rrtmg_sw_ensemble &
    (ncol    ,nlay   , &
    permuteseed, irng, n_ens   ,col_by_col   ,do_seed_permutation , &
    icld     ,ispec    ,iaer    , &
    play    ,plev    ,tlay    ,tlev    ,tsfc   , &
    cldfrac, ciwp, clwp, reic, relq, tauc, ssac, asmc, fsfc , &
    h2ovmr , o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr ,o2vmr , &
    asdir   ,asdif   ,aldir   ,aldif   , &
    kmodts, coszen  ,adjes   ,dyofyr  ,scon    ,isolvar, &
    inflgsw ,iceflgsw,liqflgsw, &
    tauaer  ,ssaaer  ,asmaer  ,ecaer   , &
    bndsolvar,indsolvar,solcycfrac, &
    swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
    swuflxspec, swdflxspec, swuflxcspec, swdflxcspec, &
    add_aero_layer, r_mu, t_mu, r_bar, t_bar)
! Modules
    use parkind, only : im => kind_im
    use rrtmg_sw_rad, only: rrtmg_sw
    use mcica_subcol_gen_sw, only: mcica_subcol_sw
    ! use parrrsw, only: nbndsw, naerec
    integer(kind=im), parameter :: nbndsw = 14
    integer(kind=im), parameter :: naerec  = 6
! Input
    integer, parameter :: rb = selected_real_kind(12)
    integer(kind=im), intent(in) :: ncol            ! number of columns
    integer(kind=im), intent(in) :: nlay            ! number of model layers
    integer(kind=im), intent(in) :: n_ens               ! number of ensemble members
    integer(kind=im), intent(in) :: col_by_col          ! whether to run col by col
    integer(kind=im), intent(in) :: do_seed_permutation ! whether to permute the random seed
    integer(kind=im), intent(in) :: permuteseed     ! if the cloud generator is called multiple times,
                                                  ! permute the seed between each call.
                                                  ! between calls for LW and SW, recommended
                                                  ! permuteseed differes by 'ngpt'
    integer(kind=im), intent(inout) :: irng         ! flag for random number generator
                                                  !  0 = kissvec
                                                  !  1 = Mersenne Twister
    integer(kind=im), intent(inout) :: icld         ! Cloud overlap method
                                                    !    0: Clear only
                                                    !    1: Random
                                                    !    2: Maximum/random
                                                    !    3: Maximum
    integer(kind=im), intent(inout) :: ispec        ! spectral ASR output flag
    integer(kind=im), intent(inout) :: iaer         ! Aerosol option flag
                                                    !    0: No aerosol
                                                    !    6: ECMWF method
                                                    !    10:Input aerosol optical
                                                    !       properties
    real(kind=rb), intent(in) :: play(ncol,nlay)    ! Layer pressures (hPa, mb)
    real(kind=rb), intent(in) :: plev(ncol,nlay+1)  ! Interface pressures (hPa, mb)
    real(kind=rb), intent(in) :: tlay(ncol,nlay)    ! Layer temperatures (K)
    real(kind=rb), intent(in) :: tlev(ncol,nlay+1)  ! Interface temperatures (K)
    real(kind=rb), intent(in) :: tsfc(ncol)         ! Surface temperature (K)
    real(kind=rb), intent(in) :: cldfrac(ncol,nlay)        ! layer cloud fraction
    real(kind=rb), intent(in) :: ciwp(ncol,nlay)           ! in-cloud ice water path
    real(kind=rb), intent(in) :: clwp(ncol,nlay)           ! in-cloud liquid water path
    real(kind=rb), intent(in) :: reic(ncol,nlay)           ! cloud ice particle size
    real(kind=rb), intent(in) :: relq(ncol,nlay)           ! cloud liquid particle size
    real(kind=rb), intent(in) :: tauc(nbndsw,ncol,nlay)    ! in-cloud optical depth
    real(kind=rb), intent(in) :: ssac(nbndsw,ncol,nlay)    ! in-cloud single scattering albedo (non-delta scaled)
    real(kind=rb), intent(in) :: asmc(nbndsw,ncol,nlay)    ! in-cloud asymmetry parameter (non-delta scaled)
    real(kind=rb), intent(in) :: fsfc(nbndsw,ncol,nlay)    ! in-cloud forward scattering fraction (non-delta scaled)
    real(kind=rb), intent(in) :: h2ovmr(ncol,nlay)  ! H2O volume mixing ratio
    real(kind=rb), intent(in) :: o3vmr(ncol,nlay)   ! O3 volume mixing ratio
    real(kind=rb), intent(in) :: co2vmr(ncol,nlay)  ! CO2 volume mixing ratio
    real(kind=rb), intent(in) :: ch4vmr(ncol,nlay)  ! Methane volume mixing ratio
    real(kind=rb), intent(in) :: n2ovmr(ncol,nlay)  ! Nitrous oxide volume mixing ratio
    real(kind=rb), intent(in) :: o2vmr(ncol,nlay)   ! Oxygen volume mixing ratio
    real(kind=rb), intent(in) :: aldif(ncol)        ! UV/vis surface albedo direct rad
    real(kind=rb), intent(in) :: aldir(ncol)        ! Near-IR surface albedo direct rad
    real(kind=rb), intent(in) :: asdif(ncol)        ! UV/vis surface albedo: diffuse rad
    real(kind=rb), intent(in) :: asdir(ncol)        ! Near-IR surface albedo: diffuse rad
    integer(kind=im), intent(in) :: kmodts
    real(kind=rb), intent(in) :: coszen(ncol)       ! Cosine of solar zenith angle
    real(kind=rb), intent(in) :: adjes(ncol)        ! Flux adjustment (Earth/Sun distance and/or zenith angle compensation)
    integer(kind=im), intent(in) :: dyofyr          ! Day of the year (used to get Earth/Sun
                                                    !  distance if adjflx not provided)
    real(kind=rb), intent(in) :: scon               ! Solar constant (W/m2)
                                                    !    Total solar irradiance averaged
                                                    !    over the solar cycle.
                                                    !    If scon = 0.0, the internal solar
                                                    !    constant, which depends on the
                                                    !    value of isolvar, will be used.
                                                    !    For isolvar=-1, scon=1368.22 Wm-2,
                                                    !    For isolvar=0,1,3, scon=1360.85 Wm-2,
                                                    !    If scon > 0.0, the internal solar
                                                    !    constant will be scaled to the
                                                    !    provided value of scon.
    integer(kind=im), intent(in) :: isolvar         ! Flag for solar variability method
                                                    !   -1 = (when scon .eq. 0.0): No solar variability
                                                    !        and no solar cycle (Kurucz solar irradiance
                                                    !        of 1368.22 Wm-2 only);
                                                    !        (when scon .ne. 0.0): Kurucz solar irradiance
                                                    !        scaled to scon and solar variability defined
                                                    !        (optional) by setting non-zero scale factors
                                                    !        for each band in bndsolvar
                                                    !    0 = (when SCON .eq. 0.0): No solar variability
                                                    !        and no solar cycle (NRLSSI2 solar constant of
                                                    !        1360.85 Wm-2 for the 100-50000 cm-1 spectral
                                                    !        range only), with facular and sunspot effects
                                                    !        fixed to the mean of Solar Cycles 13-24;
                                                    !        (when SCON .ne. 0.0): No solar variability
                                                    !        and no solar cycle (NRLSSI2 solar constant of
                                                    !        1360.85 Wm-2 for the 100-50000 cm-1 spectral
                                                    !        range only), is scaled to SCON
                                                    !    1 = Solar variability (using NRLSSI2  solar
                                                    !        model) with solar cycle contribution
                                                    !        determined by fraction of solar cycle
                                                    !        with facular and sunspot variations
                                                    !        fixed to their mean variations over the
                                                    !        average of Solar Cycles 13-24;
                                                    !        two amplitude scale factors allow
                                                    !        facular and sunspot adjustments from
                                                    !        mean solar cycle as defined by indsolvar
                                                    !    2 = Solar variability (using NRLSSI2 solar
                                                    !        model) over solar cycle determined by
                                                    !        direct specification of Mg (facular)
                                                    !        and SB (sunspot) indices provided
                                                    !        in indsolvar (scon = 0.0 only)
                                                    !    3 = (when scon .eq. 0.0): No solar variability
                                                    !        and no solar cycle (NRLSSI2 solar irradiance
                                                    !        of 1360.85 Wm-2 only);
                                                    !        (when scon .ne. 0.0): NRLSSI2 solar irradiance
                                                    !        scaled to scon and solar variability defined
                                                    !        (optional) by setting non-zero scale factors
                                                    !        for each band in bndsolvar
    real(kind=rb), intent(inout) :: indsolvar(2) ! Facular and sunspot amplitude
                                              ! scale factors (isolvar=1), or
                                              ! Mg and SB indices (isolvar=2)
                                              !    Dimensions: (2)
    real(kind=rb), intent(inout) :: bndsolvar(nbndsw) ! Solar variability scale factors
                                                   ! for each shortwave band
                                                   !    Dimensions: (nbndsw=14)
    real(kind=rb), intent(inout) :: solcycfrac   ! Fraction of averaged solar cycle (0-1)
                                              !    at current time (isolvar=1)

    integer(kind=im), intent(in) :: inflgsw         ! Flag for cloud optical properties
    integer(kind=im), intent(in) :: iceflgsw        ! Flag for ice particle specification
    integer(kind=im), intent(in) :: liqflgsw        ! Flag for liquid droplet specification
    !   These quantities are computed by McICA
    real(kind=rb), intent(in) :: tauaer(ncol,nlay,nbndsw)  ! Aerosol optical depth (iaer=10 only)
    real(kind=rb), intent(in) :: ssaaer(ncol,nlay,nbndsw)  ! Aerosol single scattering albedo (iaer=10 only)
    real(kind=rb), intent(in) :: asmaer(ncol,nlay,nbndsw)  ! Aerosol asymmetry parameter (iaer=10 only)
    real(kind=rb), intent(in) :: ecaer(ncol,nlay,naerec)   ! Aerosol optical depth at 0.55 micron (iaer=6 only)

! Output
    real(kind=rb), intent(out) :: swuflx(ncol,nlay+1)    ! Total sky shortwave upward flux (W/m2)
    real(kind=rb), intent(out) :: swdflx(ncol,nlay+1)    ! Total sky shortwave downward flux (W/m2)
    real(kind=rb), intent(out) :: swhr(ncol,nlay)        ! Total sky shortwave radiative heating rate (K/d)
    real(kind=rb), intent(out) :: swuflxc(ncol,nlay+1)   ! Clear sky shortwave upward flux (W/m2)
    real(kind=rb), intent(out) :: swdflxc(ncol,nlay+1)   ! Clear sky shortwave downward flux (W/m2)
    real(kind=rb), intent(out) :: swhrc(ncol,nlay)       ! Clear sky shortwave radiative heating rate (K/d)
    real(kind=rb), intent(out) :: swuflxspec(ncol,nlay+1,nbndsw)    ! Total sky shortwave upward flux spectrum (W/m2)
    real(kind=rb), intent(out) :: swdflxspec(ncol,nlay+1,nbndsw)    ! Total sky shortwave downward flux spectrum (W/m2)
    real(kind=rb), intent(out) :: swuflxcspec(ncol,nlay+1,nbndsw)   ! Clear sky shortwave upward flux spectrum (W/m2)
    real(kind=rb), intent(out) :: swdflxcspec(ncol,nlay+1,nbndsw)   ! Clear sky shortwave downward flux spectrum (W/m2)
    integer(kind=im), intent(in) :: add_aero_layer
    real(kind=rb), intent(in) :: r_mu(ncol,nlay,nbndsw)    ! Aerosols directional reflection
    real(kind=rb), intent(in) :: t_mu(ncol,nlay,nbndsw)    ! Aerosols directional transmission
    real(kind=rb), intent(in) :: r_bar(ncol,nlay,nbndsw)    ! Aerosols diffusive reflection
    real(kind=rb), intent(in) :: t_bar(ncol,nlay,nbndsw)    ! Aerosols diffusive transmission

    integer(kind=im) :: ind_ens, n_rrtmg_repeat, icol, ilay

!  These are not comments! Necessary directives to f2py to handle array dimensions
!f2py depend(ncol,nlay) play
!f2py depend(ncol,nlay) plev
!f2py depend(ncol,nlay) tlay
!f2py depend(ncol,nlay) tlev
!f2py depend(ncol,nlay) cldfrac,ciwp,clwp,reic,relq
!f2py depend(ncol,nlay) tauc,ssac,asmc,fsfc
!f2py depend(ncol,nlay) h2ovmr,o3vmr,co2vmr,ch4vmr,n2ovmr,o2vmr
!f2py depend(ncol) asdir,asdif,aldir,aldif,coszen, adjes, tsfc
!f2py depend(ncol,nlay) tauaer,ssaaer,asmaer,ecaer
!f2py depend(ncol,nlay) swuflx,swdflx
!f2py depend(ncol,nlay) swhr
!f2py depend(ncol,nlay) swuflxc,swdflxc
!f2py depend(ncol,nlay) swhrc
!f2py depend(ncol,nlay,nbndsw) swuflxspec,swdflxspec,swuflxcspec,swdflxcspec
!f2py depend(ncol,nlay,nbndsw) r_mu, t_mu, r_bar, t_bar
    if ((icld > 0).and.(do_seed_permutation.eq.1)) then
        n_rrtmg_repeat = n_ens
    else
        n_rrtmg_repeat = 1
    end if
    ! initialize output vars
    swuflx = 0._rb
    swdflx = 0._rb
    swhr = 0._rb
    swuflxc = 0._rb
    swdflxc = 0._rb
    swhrc = 0._rb
    swuflxspec = 0._rb
    swdflxspec = 0._rb
    swuflxcspec = 0._rb
    swdflxcspec = 0._rb
    ! These are not comments, they are using for multi-processing
    !$omp parallel do collapse(2) default(shared) &
    !$omp private(i, ilay) schedule(static)
    do icol=1, ncol
        do ind_ens=1, n_rrtmg_repeat
            call climlab_rrtmg_sw_single_col &
            (icol, ind_ens, ncol    ,nlay   , &
            permuteseed, irng, n_ens   ,col_by_col   ,do_seed_permutation , &
            icld     ,ispec    ,iaer    , &
            play    ,plev    ,tlay    ,tlev    ,tsfc   , &
            cldfrac, ciwp, clwp, reic, relq, tauc, ssac, asmc, fsfc , &
            h2ovmr , o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr ,o2vmr , &
            asdir   ,asdif   ,aldir   ,aldif   , &
            kmodts, coszen  ,adjes   ,dyofyr  ,scon    ,isolvar, &
            inflgsw ,iceflgsw,liqflgsw, &
            tauaer  ,ssaaer  ,asmaer  ,ecaer   , &
            bndsolvar,indsolvar,solcycfrac, &
            swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
            swuflxspec, swdflxspec, swuflxcspec, swdflxcspec, &
            add_aero_layer, r_mu, t_mu, r_bar, t_bar)
        end do
    end do
    !$omp end parallel do
    ! These are not comments, they are using for multi-processing
    swuflx = swuflx / n_rrtmg_repeat
    swdflx = swdflx / n_rrtmg_repeat
    swhr = swhr / n_rrtmg_repeat
    swuflxc = swuflxc / n_rrtmg_repeat
    swdflxc = swdflxc / n_rrtmg_repeat
    swhrc = swhrc / n_rrtmg_repeat
    swuflxspec = swuflxspec / n_rrtmg_repeat
    swdflxspec = swdflxspec / n_rrtmg_repeat
    swuflxcspec = swuflxcspec / n_rrtmg_repeat
    swdflxcspec = swdflxcspec / n_rrtmg_repeat
end subroutine climlab_rrtmg_sw_ensemble

subroutine climlab_rrtmg_sw_single_col&
    (icol, ind_ens, ncol    ,nlay   , &
    permuteseed, irng, n_ens   ,col_by_col   ,do_seed_permutation , &
    icld     ,ispec    ,iaer    , &
    play    ,plev    ,tlay    ,tlev    ,tsfc   , &
    cldfrac, ciwp, clwp, reic, relq, tauc, ssac, asmc, fsfc , &
    h2ovmr , o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr ,o2vmr , &
    asdir   ,asdif   ,aldir   ,aldif   , &
    kmodts, coszen  ,adjes   ,dyofyr  ,scon    ,isolvar, &
    inflgsw ,iceflgsw,liqflgsw, &
    tauaer  ,ssaaer  ,asmaer  ,ecaer   , &
    bndsolvar,indsolvar,solcycfrac, &
    swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
    swuflxspec, swdflxspec, swuflxcspec, swdflxcspec, &
    add_aero_layer, r_mu, t_mu, r_bar, t_bar)
! Modules
    use parkind, only : im => kind_im
    use rrtmg_sw_rad, only: rrtmg_sw
    use mcica_subcol_gen_sw, only: mcica_subcol_sw
    ! use parrrsw, only: nbndsw, ngptsw, naerec
    integer(kind=im), parameter :: nbndsw = 14
    integer(kind=im), parameter :: naerec  = 6
    integer(kind=im), parameter :: ngptsw = 112
! Input
    integer, parameter :: rb = selected_real_kind(12)
    integer(kind=im), intent(in) :: icol            ! index for calculation
    integer(kind=im), intent(in) :: ind_ens         ! index for ensamble sampling
    integer(kind=im), intent(in) :: ncol            ! number of columns
    integer(kind=im), intent(in) :: nlay            ! number of model layers
    integer(kind=im), intent(in) :: n_ens               ! number of ensemble members
    integer(kind=im), intent(in) :: col_by_col          ! whether to run col by col
    integer(kind=im), intent(in) :: do_seed_permutation ! whether to permute the random seed
    integer(kind=im), intent(in) :: permuteseed     ! if the cloud generator is called multiple times,
                                                  ! permute the seed between each call.
                                                  ! between calls for LW and SW, recommended
                                                  ! permuteseed differes by 'ngpt'
    integer(kind=im), intent(inout) :: irng         ! flag for random number generator
                                                  !  0 = kissvec
                                                  !  1 = Mersenne Twister
    integer(kind=im), intent(inout) :: icld         ! Cloud overlap method
                                                    !    0: Clear only
                                                    !    1: Random
                                                    !    2: Maximum/random
                                                    !    3: Maximum
    integer(kind=im), intent(inout) :: ispec        ! spectral ASR output flag
    integer(kind=im), intent(inout) :: iaer         ! Aerosol option flag
                                                    !    0: No aerosol
                                                    !    6: ECMWF method
                                                    !    10:Input aerosol optical
                                                    !       properties
    real(kind=rb), intent(in) :: play(ncol,nlay)    ! Layer pressures (hPa, mb)
    real(kind=rb), intent(in) :: plev(ncol,nlay+1)  ! Interface pressures (hPa, mb)
    real(kind=rb), intent(in) :: tlay(ncol,nlay)    ! Layer temperatures (K)
    real(kind=rb), intent(in) :: tlev(ncol,nlay+1)  ! Interface temperatures (K)
    real(kind=rb), intent(in) :: tsfc(ncol)         ! Surface temperature (K)
    real(kind=rb), intent(in) :: cldfrac(ncol,nlay)        ! layer cloud fraction
    real(kind=rb), intent(in) :: ciwp(ncol,nlay)           ! in-cloud ice water path
    real(kind=rb), intent(in) :: clwp(ncol,nlay)           ! in-cloud liquid water path
    real(kind=rb), intent(in) :: reic(ncol,nlay)           ! cloud ice particle size
    real(kind=rb), intent(in) :: relq(ncol,nlay)           ! cloud liquid particle size
    real(kind=rb), intent(in) :: tauc(nbndsw,ncol,nlay)    ! in-cloud optical depth
    real(kind=rb), intent(in) :: ssac(nbndsw,ncol,nlay)    ! in-cloud single scattering albedo (non-delta scaled)
    real(kind=rb), intent(in) :: asmc(nbndsw,ncol,nlay)    ! in-cloud asymmetry parameter (non-delta scaled)
    real(kind=rb), intent(in) :: fsfc(nbndsw,ncol,nlay)    ! in-cloud forward scattering fraction (non-delta scaled)
    real(kind=rb), intent(in) :: h2ovmr(ncol,nlay)  ! H2O volume mixing ratio
    real(kind=rb), intent(in) :: o3vmr(ncol,nlay)   ! O3 volume mixing ratio
    real(kind=rb), intent(in) :: co2vmr(ncol,nlay)  ! CO2 volume mixing ratio
    real(kind=rb), intent(in) :: ch4vmr(ncol,nlay)  ! Methane volume mixing ratio
    real(kind=rb), intent(in) :: n2ovmr(ncol,nlay)  ! Nitrous oxide volume mixing ratio
    real(kind=rb), intent(in) :: o2vmr(ncol,nlay)   ! Oxygen volume mixing ratio
    real(kind=rb), intent(in) :: aldif(ncol)        ! UV/vis surface albedo direct rad
    real(kind=rb), intent(in) :: aldir(ncol)        ! Near-IR surface albedo direct rad
    real(kind=rb), intent(in) :: asdif(ncol)        ! UV/vis surface albedo: diffuse rad
    real(kind=rb), intent(in) :: asdir(ncol)        ! Near-IR surface albedo: diffuse rad
    integer(kind=im), intent(in) :: kmodts
    real(kind=rb), intent(in) :: coszen(ncol)       ! Cosine of solar zenith angle
    real(kind=rb), intent(in) :: adjes(ncol)        ! Flux adjustment (Earth/Sun distance and/or zenith angle compensation)
    integer(kind=im), intent(in) :: dyofyr          ! Day of the year (used to get Earth/Sun
                                                    !  distance if adjflx not provided)
    real(kind=rb), intent(in) :: scon               ! Solar constant (W/m2)
                                                    !    Total solar irradiance averaged
                                                    !    over the solar cycle.
                                                    !    If scon = 0.0, the internal solar
                                                    !    constant, which depends on the
                                                    !    value of isolvar, will be used.
                                                    !    For isolvar=-1, scon=1368.22 Wm-2,
                                                    !    For isolvar=0,1,3, scon=1360.85 Wm-2,
                                                    !    If scon > 0.0, the internal solar
                                                    !    constant will be scaled to the
                                                    !    provided value of scon.
    integer(kind=im), intent(in) :: isolvar         ! Flag for solar variability method
                                                    !   -1 = (when scon .eq. 0.0): No solar variability
                                                    !        and no solar cycle (Kurucz solar irradiance
                                                    !        of 1368.22 Wm-2 only);
                                                    !        (when scon .ne. 0.0): Kurucz solar irradiance
                                                    !        scaled to scon and solar variability defined
                                                    !        (optional) by setting non-zero scale factors
                                                    !        for each band in bndsolvar
                                                    !    0 = (when SCON .eq. 0.0): No solar variability
                                                    !        and no solar cycle (NRLSSI2 solar constant of
                                                    !        1360.85 Wm-2 for the 100-50000 cm-1 spectral
                                                    !        range only), with facular and sunspot effects
                                                    !        fixed to the mean of Solar Cycles 13-24;
                                                    !        (when SCON .ne. 0.0): No solar variability
                                                    !        and no solar cycle (NRLSSI2 solar constant of
                                                    !        1360.85 Wm-2 for the 100-50000 cm-1 spectral
                                                    !        range only), is scaled to SCON
                                                    !    1 = Solar variability (using NRLSSI2  solar
                                                    !        model) with solar cycle contribution
                                                    !        determined by fraction of solar cycle
                                                    !        with facular and sunspot variations
                                                    !        fixed to their mean variations over the
                                                    !        average of Solar Cycles 13-24;
                                                    !        two amplitude scale factors allow
                                                    !        facular and sunspot adjustments from
                                                    !        mean solar cycle as defined by indsolvar
                                                    !    2 = Solar variability (using NRLSSI2 solar
                                                    !        model) over solar cycle determined by
                                                    !        direct specification of Mg (facular)
                                                    !        and SB (sunspot) indices provided
                                                    !        in indsolvar (scon = 0.0 only)
                                                    !    3 = (when scon .eq. 0.0): No solar variability
                                                    !        and no solar cycle (NRLSSI2 solar irradiance
                                                    !        of 1360.85 Wm-2 only);
                                                    !        (when scon .ne. 0.0): NRLSSI2 solar irradiance
                                                    !        scaled to scon and solar variability defined
                                                    !        (optional) by setting non-zero scale factors
                                                    !        for each band in bndsolvar
    real(kind=rb), intent(inout) :: indsolvar(2) ! Facular and sunspot amplitude
                                              ! scale factors (isolvar=1), or
                                              ! Mg and SB indices (isolvar=2)
                                              !    Dimensions: (2)
    real(kind=rb), intent(inout) :: bndsolvar(nbndsw) ! Solar variability scale factors
                                                   ! for each shortwave band
                                                   !    Dimensions: (nbndsw=14)
    real(kind=rb), intent(inout) :: solcycfrac   ! Fraction of averaged solar cycle (0-1)
                                              !    at current time (isolvar=1)

    integer(kind=im), intent(in) :: inflgsw         ! Flag for cloud optical properties
    integer(kind=im), intent(in) :: iceflgsw        ! Flag for ice particle specification
    integer(kind=im), intent(in) :: liqflgsw        ! Flag for liquid droplet specification
    !   These quantities are computed by McICA
    real(kind=rb), intent(in) :: tauaer(ncol,nlay,nbndsw)  ! Aerosol optical depth (iaer=10 only)
    real(kind=rb), intent(in) :: ssaaer(ncol,nlay,nbndsw)  ! Aerosol single scattering albedo (iaer=10 only)
    real(kind=rb), intent(in) :: asmaer(ncol,nlay,nbndsw)  ! Aerosol asymmetry parameter (iaer=10 only)
    real(kind=rb), intent(in) :: ecaer(ncol,nlay,naerec)   ! Aerosol optical depth at 0.55 micron (iaer=6 only)

! Output
    real(kind=rb), intent(inout) :: swuflx(ncol,nlay+1)    ! Total sky shortwave upward flux (W/m2)
    real(kind=rb), intent(inout) :: swdflx(ncol,nlay+1)    ! Total sky shortwave downward flux (W/m2)
    real(kind=rb), intent(inout) :: swhr(ncol,nlay)        ! Total sky shortwave radiative heating rate (K/d)
    real(kind=rb), intent(inout) :: swuflxc(ncol,nlay+1)   ! Clear sky shortwave upward flux (W/m2)
    real(kind=rb), intent(inout) :: swdflxc(ncol,nlay+1)   ! Clear sky shortwave downward flux (W/m2)
    real(kind=rb), intent(inout) :: swhrc(ncol,nlay)       ! Clear sky shortwave radiative heating rate (K/d)
    real(kind=rb), intent(inout) :: swuflxspec(ncol,nlay+1,nbndsw)    ! Total sky shortwave upward flux spectrum (W/m2)
    real(kind=rb), intent(inout) :: swdflxspec(ncol,nlay+1,nbndsw)    ! Total sky shortwave downward flux spectrum (W/m2)
    real(kind=rb), intent(inout) :: swuflxcspec(ncol,nlay+1,nbndsw)   ! Clear sky shortwave upward flux spectrum (W/m2)
    real(kind=rb), intent(inout) :: swdflxcspec(ncol,nlay+1,nbndsw)   ! Clear sky shortwave downward flux spectrum (W/m2)
    integer(kind=im), intent(in) :: add_aero_layer
    real(kind=rb), intent(in) :: r_mu(ncol,nlay,nbndsw)    ! Aerosols directional reflection
    real(kind=rb), intent(in) :: t_mu(ncol,nlay,nbndsw)    ! Aerosols directional transmission
    real(kind=rb), intent(in) :: r_bar(ncol,nlay,nbndsw)    ! Aerosols diffusive reflection
    real(kind=rb), intent(in) :: t_bar(ncol,nlay,nbndsw)    ! Aerosols diffusive transmission

    ! Column internal variables
    real(kind=rb) :: play1(1,nlay)    ! Layer pressures (hPa, mb)
    real(kind=rb) :: plev1(1,nlay+1)  ! Interface pressures (hPa, mb)
    real(kind=rb) :: tlay1(1,nlay)    ! Layer temperatures (K)
    real(kind=rb) :: tlev1(1,nlay+1)  ! Interface temperatures (K)
    real(kind=rb) :: tsfc1(1)         ! Surface temperature (K)
    real(kind=rb) :: cldfrac1(1,nlay)        ! layer cloud fraction
    real(kind=rb) :: ciwp1(1,nlay)           ! in-cloud ice water path
    real(kind=rb) :: clwp1(1,nlay)           ! in-cloud liquid water path
    real(kind=rb) :: reic1(1,nlay)           ! cloud ice particle size
    real(kind=rb) :: relq1(1,nlay)           ! cloud liquid particle size
    real(kind=rb) :: tauc1(nbndsw,1,nlay)    ! in-cloud optical depth
    real(kind=rb) :: ssac1(nbndsw,1,nlay)    ! in-cloud single scattering albedo (non-delta scaled)
    real(kind=rb) :: asmc1(nbndsw,1,nlay)    ! in-cloud asymmetry parameter (non-delta scaled)
    real(kind=rb) :: fsfc1(nbndsw,1,nlay)    ! in-cloud forward scattering fraction (non-delta scaled)
    real(kind=rb) :: h2ovmr1(1,nlay)  ! H2O volume mixing ratio
    real(kind=rb) :: o3vmr1(1,nlay)   ! O3 volume mixing ratio
    real(kind=rb) :: co2vmr1(1,nlay)  ! CO2 volume mixing ratio
    real(kind=rb) :: ch4vmr1(1,nlay)  ! Methane volume mixing ratio
    real(kind=rb) :: n2ovmr1(1,nlay)  ! Nitrous oxide volume mixing ratio
    real(kind=rb) :: o2vmr1(1,nlay)   ! Oxygen volume mixing ratio
    real(kind=rb) :: aldif1(1)        ! UV/vis surface albedo direct rad
    real(kind=rb) :: aldir1(1)        ! Near-IR surface albedo direct rad
    real(kind=rb) :: asdif1(1)        ! UV/vis surface albedo: diffuse rad
    real(kind=rb) :: asdir1(1)        ! Near-IR surface albedo: diffuse rad
    real(kind=rb) :: coszen1(1)       ! Cosine of solar zenith angle
    real(kind=rb) :: adjes1(1)        ! Flux adjustment (Earth/Sun distance and/or zenith angle compensation)
    real(kind=rb) :: tauaer1(1,nlay,nbndsw)  ! Aerosol optical depth (iaer=10 only)
    real(kind=rb) :: ssaaer1(1,nlay,nbndsw)  ! Aerosol single scattering albedo (iaer=10 only)
    real(kind=rb) :: asmaer1(1,nlay,nbndsw)  ! Aerosol asymmetry parameter (iaer=10 only)
    real(kind=rb) :: ecaer1(1,nlay,naerec)   ! Aerosol optical depth at 0.55 micron (iaer=6 only)
    real(kind=rb) :: r_mu1(1,nlay,nbndsw)    ! Aerosols directional reflection
    real(kind=rb) :: t_mu1(1,nlay,nbndsw)    ! Aerosols directional transmission
    real(kind=rb) :: r_bar1(1,nlay,nbndsw)    ! Aerosols diffusive reflection
    real(kind=rb) :: t_bar1(1,nlay,nbndsw)    ! Aerosols diffusive transmission

    real(kind=rb) :: swuflx1(1,nlay+1)    ! Total sky shortwave upward flux (W/m2)
    real(kind=rb) :: swdflx1(1,nlay+1)    ! Total sky shortwave downward flux (W/m2)
    real(kind=rb) :: swhr1(1,nlay)        ! Total sky shortwave radiative heating rate (K/d)
    real(kind=rb) :: swuflxc1(1,nlay+1)   ! Clear sky shortwave upward flux (W/m2)
    real(kind=rb) :: swdflxc1(1,nlay+1)   ! Clear sky shortwave downward flux (W/m2)
    real(kind=rb) :: swhrc1(1,nlay)       ! Clear sky shortwave radiative heating rate (K/d)
    real(kind=rb) :: swuflxspec1(1,nlay+1,nbndsw)    ! Total sky shortwave upward flux spectrum (W/m2)
    real(kind=rb) :: swdflxspec1(1,nlay+1,nbndsw)    ! Total sky shortwave downward flux spectrum (W/m2)
    real(kind=rb) :: swuflxcspec1(1,nlay+1,nbndsw)   ! Clear sky shortwave upward flux spectrum (W/m2)
    real(kind=rb) :: swdflxcspec1(1,nlay+1,nbndsw)   ! Clear sky shortwave downward flux spectrum (W/m2)

    ! These quantities are computed by McICA
    real(kind=rb) :: cldfmcl1(ngptsw,1,nlay)   ! cloud fraction [mcica]
    real(kind=rb) :: ciwpmcl1(ngptsw,1,nlay)   ! in-cloud ice water path [mcica]
    real(kind=rb) :: clwpmcl1(ngptsw,1,nlay)   ! in-cloud liquid water path [mcica]
    real(kind=rb) :: reicmcl1(1,nlay)          ! ice partcle size (microns)
    real(kind=rb) :: relqmcl1(1,nlay)          ! liquid particle size (microns)
    real(kind=rb) :: taucmcl1(ngptsw,1,nlay)   ! in-cloud optical depth [mcica]
    real(kind=rb) :: ssacmcl1(ngptsw,1,nlay)   ! in-cloud single scattering albedo [mcica]
    real(kind=rb) :: asmcmcl1(ngptsw,1,nlay)   ! in-cloud asymmetry parameter [mcica]
    real(kind=rb) :: fsfcmcl1(ngptsw,1,nlay)   ! in-cloud forward scattering fraction [mcica]

    integer(kind=im) :: seed, ilay, shift

    ! Arrays initialization
    swuflx1 = 0._rb
    swdflx1 = 0._rb
    swhr1 = 0._rb
    swuflxc1 = 0._rb
    swdflxc1 = 0._rb
    swhrc1 = 0._rb
    swuflxspec1 = 0._rb
    swdflxspec1 = 0._rb
    swuflxcspec1 = 0._rb
    swdflxcspec1 = 0._rb

    cldfmcl1 = 0._rb
    ciwpmcl1 = 0._rb
    clwpmcl1 = 0._rb
    reicmcl1 = 0._rb
    relqmcl1 = 0._rb
    taucmcl1 = 0._rb
    ssacmcl1 = 0._rb
    asmcmcl1 = 0._rb
    fsfcmcl1 = 0._rb
    !  Call the Monte Carlo Independent Column Approximation (McICA, Pincus et al., JC, 2003)
    if (do_seed_permutation.eq.1) then
        if (col_by_col.eq.0) then
            shift = ind_ens * ncol + icol - 1
        else 
            shift = ind_ens
        end if
        seed = permuteseed + shift * ngptsw
    else
        seed = permuteseed
    end if
    play1(1,:) = play(icol,:)
    plev1(1,:) = plev(icol,:)
    tlay1(1,:) = tlay(icol,:)
    tlev1(1,:) = tlev(icol,:)
    tsfc1(1) = tsfc(icol)
    cldfrac1(1,:) = cldfrac(icol,:)
    ciwp1(1,:) = ciwp(1,:)
    clwp1(1,:) = clwp(1,:)
    reic1(1,:) = reic(1,:)
    relq1(1,:) = relq(1,:)
    h2ovmr1(1,:) = h2ovmr(icol,:)
    o3vmr1(1,:) = o3vmr(icol,:)
    co2vmr1(1,:) = co2vmr(icol,:)
    ch4vmr1(1,:) = ch4vmr(icol,:)
    n2ovmr1(1,:) = n2ovmr(icol,:)
    o2vmr1(1,:) = o2vmr(icol,:)
    aldif1(1) = aldif(icol)
    aldir1(1) = aldir(icol)
    asdif1(1) = asdif(icol)
    asdir1(1) = asdir(icol)
    coszen1(1) = coszen(icol)
    adjes1(1) = adjes(icol)
    do ilay=1, nlay
        tauc1(:,1,ilay) = tauc(:,icol,ilay)
        ssac1(:,1,ilay) = ssac(:,icol,ilay)
        asmc1(:,1,ilay) = asmc(:,icol,ilay)
        fsfc1(:,1,ilay) = fsfc(:,icol,ilay)
        tauaer1(1,ilay,:) = tauaer(icol,ilay,:)
        ssaaer1(1,ilay,:) = ssaaer(icol,ilay,:)
        asmaer1(1,ilay,:) = asmaer(icol,ilay,:)
        ecaer1(1,ilay,:) = ecaer(icol,ilay,:)
        r_mu1(1,ilay,:) = r_mu(icol,ilay,:)
        t_mu1(1,ilay,:) = t_mu(icol,ilay,:)
        r_bar1(1,ilay,:) = r_bar(icol,ilay,:)
        t_bar1(1,ilay,:) = t_bar(icol,ilay,:)
    end do
    call mcica_subcol_sw(1, 1, nlay, icld, seed, irng, play1, &
        cldfrac1, ciwp1, clwp1, reic1, relq1, &
        tauc1, ssac1, asmc1, fsfc1, &
        cldfmcl1, ciwpmcl1, clwpmcl1, reicmcl1, relqmcl1, &
        taucmcl1, ssacmcl1, asmcmcl1, fsfcmcl1)

    !  Call the RRTMG_SW driver to compute radiative fluxes
    call rrtmg_sw(1    ,nlay    ,icld    ,ispec    ,iaer    , &
            play1    ,plev1    ,tlay1    ,tlev1    ,tsfc1   , &
            h2ovmr1 , o3vmr1   ,co2vmr1  ,ch4vmr1  ,n2ovmr1 ,o2vmr1 , &
            asdir1   ,asdif1   ,aldir1   ,aldif1   , &
            kmodts, coszen1  ,adjes1   ,dyofyr  ,scon    ,isolvar, &
            inflgsw ,iceflgsw,liqflgsw,cldfmcl1 , &
            taucmcl1 ,ssacmcl1 ,asmcmcl1 ,fsfcmcl1 , &
            ciwpmcl1 ,clwpmcl1 ,reicmcl1 ,relqmcl1 , &
            tauaer1  ,ssaaer1  ,asmaer1  ,ecaer1   , &
            swuflx1  ,swdflx1  ,swhr1    ,swuflxc1 ,swdflxc1 ,swhrc1, &
            swuflxspec1, swdflxspec1, swuflxcspec1, swdflxcspec1, &
            bndsolvar,indsolvar,solcycfrac, &
            add_aero_layer, r_mu1, t_mu1, r_bar1, t_bar1)
    swuflx(icol,:) = swuflx(icol,:) + swuflx1(1,:)
    swdflx(icol,:) = swdflx(icol,:) + swdflx1(1,:)
    swhr(icol,:) = swhr(icol,:) + swhr1(1,:)
    swuflxc(icol,:) = swuflxc(icol,:) + swuflxc1(1,:)
    swdflxc(icol,:) = swdflxc(icol,:) + swdflxc1(1,:)
    swhrc(icol,:) = swhrc(icol,:) + swhrc1(1,:)
    swuflxspec(icol,:,:) = swuflxspec(icol,:,:) + swuflxspec1(1,:,:)
    swdflxspec(icol,:,:) = swdflxspec(icol,:,:) + swdflxspec1(1,:,:)
    swuflxcspec(icol,:,:) = swuflxcspec(icol,:,:) + swuflxcspec1(1,:,:)
    swdflxcspec(icol,:,:) = swdflxcspec(icol,:,:) + swdflxcspec1(1,:,:)
end subroutine climlab_rrtmg_sw_single_col
