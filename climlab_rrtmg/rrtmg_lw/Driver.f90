!  CLIMLAB driver for RRTMG_LW radiation
!
!   using the latest version RRTMG_LW v4.85
!
!  This is a lightweight driver designed for wrapping with f2py
!  and calling from Python
!  Refer to RRTMG_LW source code for more documentation
!
!  Three functions are exposed here:
!   climlab_rrtmg_lw_ini (wrapper for rrtmg_lw_ini)
!   climlab_mcica_subcol_lw (wrapper for mcica_subcol_lw)
!   climlab_rrtmg_lw  (wrapper for rrtmg_lw)
!
!   The call signature for each of these is nearly identical to its
!   equivalent in the RRTMG_LW code.
!
!   See the python module climlab/radiation/rrtm/rrtmg_lw.py
!    to see how these are called from Python
!
!  Brian Rose
!  brose@albany.edu
!

subroutine climlab_rrtmg_lw_ini(cpdair)
    ! Modules
    use rrtmg_lw_init, only: rrtmg_lw_ini
    ! Input
    integer, parameter :: rb = selected_real_kind(12)
    real(kind=rb), intent(in) :: cpdair    ! Specific heat capacity of dry air
                                            ! at constant pressure at 273 K
                                            ! (J kg-1 K-1)
    ! Call the initialization routine
    call rrtmg_lw_ini(cpdair)
end subroutine climlab_rrtmg_lw_ini


subroutine climlab_mcica_subcol_lw &
    (ncol, nlay, icld, permuteseed, irng, play, &
    cldfrac, ciwp, clwp, reic, relq, tauc, cldfmcl, &
    ciwpmcl, clwpmcl, reicmcl, relqmcl, taucmcl)

    ! Modules
    use parkind, only : im => kind_im
    use mcica_subcol_gen_lw, only: mcica_subcol_lw
    !use parrrtm, only: nbndlw, ngptlw
    integer(kind=im), parameter :: nbndlw = 16
    integer(kind=im), parameter :: ngptlw = 140  

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
    real(kind=rb), intent(in) :: cldfrac(ncol,nlay)      ! layer cloud fraction
    real(kind=rb), intent(in) :: ciwp(ncol,nlay)         ! in-cloud ice water path
    real(kind=rb), intent(in) :: clwp(ncol,nlay)         ! in-cloud liquid water path
    real(kind=rb), intent(in) :: reic(ncol,nlay)         ! cloud ice particle size
    real(kind=rb), intent(in) :: relq(ncol,nlay)         ! cloud liquid particle size
    real(kind=rb), intent(in) :: tauc(nbndlw,ncol,nlay)  ! in-cloud optical depth
! Output
    !   These quantities are computed by McICA
    real(kind=rb), intent(out) :: cldfmcl(ngptlw,ncol,nlay)    ! cloud fraction [mcica]
    real(kind=rb), intent(out) :: ciwpmcl(ngptlw,ncol,nlay)    ! in-cloud ice water path [mcica]
    real(kind=rb), intent(out) :: clwpmcl(ngptlw,ncol,nlay)    ! in-cloud liquid water path [mcica]
    real(kind=rb), intent(out) :: reicmcl(ncol,nlay)           ! ice partcle size (microns)  [mcica]
    real(kind=rb), intent(out) :: relqmcl(ncol,nlay)           ! liquid particle size (microns) [mcica]
    real(kind=rb), intent(out) :: taucmcl(ngptlw,ncol,nlay)    ! in-cloud optical depth [mcica]
!  These are not comments! Necessary directives to f2py to handle array dimensions
!f2py depend(ncol,nlay) play,
!f2py depend(ncol,nlay) cldfrac,ciwp,clwp,reic,relq
!f2py depend(ncol,nlay) tauc
!f2py depend(ncol,nlay) cldfmcl,ciwpmcl,clwpmcl,taucmcl
!f2py depend(ncol,nlay) reicmcl,relqmcl

    ! Call the Monte Carlo Independent Column Approximation
    !   (McICA, Pincus et al., JC, 2003)
    call mcica_subcol_lw(1, ncol, nlay, icld, permuteseed, irng, play, &
                       cldfrac, ciwp, clwp, reic, relq, tauc, cldfmcl, &
                       ciwpmcl, clwpmcl, reicmcl, relqmcl, taucmcl)

end subroutine climlab_mcica_subcol_lw


subroutine climlab_rrtmg_lw &
    (ncol    ,nlay    ,icld     , ispec   , idrv    , &
    play    , plev    , tlay    , tlev    , tsfc    , &
    h2ovmr  , o3vmr   , co2vmr  , ch4vmr  , n2ovmr  , o2vmr , &
    cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr , emis    , &
    inflglw , iceflglw, liqflglw, cldfmcl , &
    taucmcl , ciwpmcl , clwpmcl , reicmcl , relqmcl , tauaer  , &
    olr_sr  , uflx    , dflx    , hr      , uflxc   , dflxc,  hrc, &
    duflx_dt,duflxc_dt)

! Modules
    use parkind, only : im => kind_im
    use rrtmg_lw_rad, only: rrtmg_lw
    !use parrrtm, only: nbndlw, ngptlw
    integer(kind=im), parameter :: nbndlw = 16
    integer(kind=im), parameter :: ngptlw = 140  

! Input
    integer, parameter :: rb = selected_real_kind(12)
    integer(kind=im), intent(in) :: ncol            ! number of columns
    integer(kind=im), intent(in) :: nlay            ! number of model layers
    integer(kind=im), intent(inout) :: icld         ! Cloud overlap method
    integer(kind=im), intent(inout) :: ispec        ! spectral OLR output flag
    integer(kind=im), intent(in) :: idrv            ! Flag for calculation of dFdT, the change
                                                    !    in upward flux as a function of
                                                    !    surface temperature [0=off, 1=on]
                                                    !    0: Normal forward calculation
                                                    !    1: Normal forward calculation with
                                                    !       duflx_dt and duflxc_dt output
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
    real(kind=rb), intent(in) :: cfc11vmr(ncol,nlay)  ! CFC11 volume mixing ratio
    real(kind=rb), intent(in) :: cfc12vmr(ncol,nlay)  ! CFC12 volume mixing ratio
    real(kind=rb), intent(in) :: cfc22vmr(ncol,nlay)  ! CFC22 volume mixing ratio
    real(kind=rb), intent(in) :: ccl4vmr(ncol,nlay)   ! CCL4 volume mixing ratio
    real(kind=rb), intent(in) :: emis(ncol,nbndlw)  ! Surface emissivity
    integer(kind=im), intent(in) :: inflglw         ! Flag for cloud optical properties
    integer(kind=im), intent(in) :: iceflglw       ! Flag for ice particle specification
    integer(kind=im), intent(in) :: liqflglw        ! Flag for liquid droplet specification
    real(kind=rb), intent(in) :: tauaer(ncol,nlay,nbndlw) ! aerosol optical depth at mid-point of LW spectral bands
    !   These quantities are computed by McICA
    real(kind=rb), intent(in) :: cldfmcl(ngptlw,ncol,nlay)    ! cloud fraction [mcica]
    real(kind=rb), intent(in) :: ciwpmcl(ngptlw,ncol,nlay)    ! in-cloud ice water path [mcica]
    real(kind=rb), intent(in) :: clwpmcl(ngptlw,ncol,nlay)    ! in-cloud liquid water path [mcica]
    real(kind=rb), intent(in) :: reicmcl(ncol,nlay)           ! ice partcle size (microns)  [mcica]
    real(kind=rb), intent(in) :: relqmcl(ncol,nlay)           ! liquid particle size (microns) [mcica]
    real(kind=rb), intent(in) :: taucmcl(ngptlw,ncol,nlay)    ! in-cloud optical depth [mcica]

! Output
    real(kind=rb), intent(out) :: olr_sr(ncol,nbndlw)    ! Spectrally-decomposed OLR (W/m2)
    real(kind=rb), intent(out) :: uflx(ncol,nlay+1)      ! Total sky longwave upward flux (W/m2)
    real(kind=rb), intent(out) :: dflx(ncol,nlay+1)      ! Total sky longwave downward flux (W/m2)
    real(kind=rb), intent(out) :: hr(ncol,nlay)          ! Total sky longwave radiative heating rate (K/d)
    real(kind=rb), intent(out) :: uflxc(ncol,nlay+1)     ! Clear sky longwave upward flux (W/m2)
    real(kind=rb), intent(out) :: dflxc(ncol,nlay+1)     ! Clear sky longwave downward flux (W/m2)
    real(kind=rb), intent(out) :: hrc(ncol,nlay)         ! Clear sky longwave radiative heating rate (K/d)

    real(kind=rb), intent(out) :: duflx_dt(ncol,nlay+1)  ! change in upward longwave flux (w/m2/K)
                                                         ! with respect to surface temperature
    real(kind=rb), intent(out) :: duflxc_dt(ncol,nlay+1) ! change in clear sky upward longwave flux (w/m2/K)
                                                         ! with respect to surface temperature
    real(kind=rb) :: uflxspec(ncol,nlay+1,nbndlw) ! Total sky longwave upward flux spectrum (W/m2)
    real(kind=rb) :: dflxspec(ncol,nlay+1,nbndlw) ! Total sky longwave downward flux spectrum (W/m2)
    real(kind=rb) :: uflxcspec(ncol,nlay+1,nbndlw) ! Clear sky longwave upward flux spectrum (W/m2)
    real(kind=rb) :: dflxcspec(ncol,nlay+1,nbndlw) ! Clear sky longwave downward flux spectrum (W/m2)

!  These are not comments! Necessary directives to f2py to handle array dimensions
!f2py depend(ncol,nlay) play, plev, tlay, tlev
!f2py depend(ncol,nlay) h2ovmr,o3vmr,co2vmr,ch4vmr,n2ovmr,o2vmr
!f2py depend(ncol,nlay) cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr
!f2py depend(ncol) tsfc, emis
!f2py depend(ncol) olr_sr
!f2py depend(ncol,nlay) tauaer
!f2py depend(ncol,nlay) cldfmcl,ciwpmcl,clwpmcl,taucmcl
!f2py depend(ncol,nlay) reicmcl,relqmcl
!f2py depend(ncol,nlay) uflx,dflx,hr,uflxc,dflxc,hrc,duflx_dt,duflxc_dt

    !  Call the RRTMG_LW driver to compute radiative fluxes
    call rrtmg_lw(ncol    ,nlay    ,icld    ,ispec           ,idrv  , &
             play    , plev    , tlay    , tlev    , tsfc    , &
             h2ovmr  , o3vmr   , co2vmr  , ch4vmr  , n2ovmr  , o2vmr , &
             cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr , emis    , &
             inflglw , iceflglw, liqflglw, cldfmcl , &
             taucmcl , ciwpmcl , clwpmcl , reicmcl , relqmcl , tauaer  , &
             olr_sr  , uflx, dflx, hr, uflxc, dflxc,  hrc, &
             duflx_dt,duflxc_dt, uflxspec, dflxspec, uflxcspec, dflxcspec)

end subroutine climlab_rrtmg_lw

subroutine climlab_rrtmg_lw_expanded &
    (ncol    ,nlay    ,icld     , ispec   , idrv    , &
    play    , plev    , tlay    , tlev    , tsfc    , &
    h2ovmr  , o3vmr   , co2vmr  , ch4vmr  , n2ovmr  , o2vmr , &
    cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr , emis    , &
    inflglw , iceflglw, liqflglw, cldfmcl , &
    taucmcl , ciwpmcl , clwpmcl , reicmcl , relqmcl , tauaer  , &
    olr_sr  , uflx    , dflx    , hr      , uflxc   , dflxc,  hrc, &
    duflx_dt,duflxc_dt, uflxspec, dflxspec, uflxcspec, dflxcspec)

! Modules
    use parkind, only : im => kind_im
    use rrtmg_lw_rad, only: rrtmg_lw
    !use parrrtm, only: nbndlw, ngptlw
    integer(kind=im), parameter :: nbndlw = 16
    integer(kind=im), parameter :: ngptlw = 140  

! Input
    integer, parameter :: rb = selected_real_kind(12)
    integer(kind=im), intent(in) :: ncol            ! number of columns
    integer(kind=im), intent(in) :: nlay            ! number of model layers
    integer(kind=im), intent(inout) :: icld         ! Cloud overlap method
    integer(kind=im), intent(inout) :: ispec        ! spectral OLR output flag
    integer(kind=im), intent(in) :: idrv            ! Flag for calculation of dFdT, the change
                                                    !    in upward flux as a function of
                                                    !    surface temperature [0=off, 1=on]
                                                    !    0: Normal forward calculation
                                                    !    1: Normal forward calculation with
                                                    !       duflx_dt and duflxc_dt output
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
    real(kind=rb), intent(in) :: cfc11vmr(ncol,nlay)  ! CFC11 volume mixing ratio
    real(kind=rb), intent(in) :: cfc12vmr(ncol,nlay)  ! CFC12 volume mixing ratio
    real(kind=rb), intent(in) :: cfc22vmr(ncol,nlay)  ! CFC22 volume mixing ratio
    real(kind=rb), intent(in) :: ccl4vmr(ncol,nlay)   ! CCL4 volume mixing ratio
    real(kind=rb), intent(in) :: emis(ncol,nbndlw)  ! Surface emissivity
    integer(kind=im), intent(in) :: inflglw         ! Flag for cloud optical properties
    integer(kind=im), intent(in) :: iceflglw       ! Flag for ice particle specification
    integer(kind=im), intent(in) :: liqflglw        ! Flag for liquid droplet specification
    real(kind=rb), intent(in) :: tauaer(ncol,nlay,nbndlw) ! aerosol optical depth at mid-point of LW spectral bands
    !   These quantities are computed by McICA
    real(kind=rb), intent(in) :: cldfmcl(ngptlw,ncol,nlay)    ! cloud fraction [mcica]
    real(kind=rb), intent(in) :: ciwpmcl(ngptlw,ncol,nlay)    ! in-cloud ice water path [mcica]
    real(kind=rb), intent(in) :: clwpmcl(ngptlw,ncol,nlay)    ! in-cloud liquid water path [mcica]
    real(kind=rb), intent(in) :: reicmcl(ncol,nlay)           ! ice partcle size (microns)  [mcica]
    real(kind=rb), intent(in) :: relqmcl(ncol,nlay)           ! liquid particle size (microns) [mcica]
    real(kind=rb), intent(in) :: taucmcl(ngptlw,ncol,nlay)    ! in-cloud optical depth [mcica]

! Output
    real(kind=rb), intent(out) :: olr_sr(ncol,nbndlw)    ! Spectrally-decomposed OLR (W/m2)
    real(kind=rb), intent(out) :: uflx(ncol,nlay+1)      ! Total sky longwave upward flux (W/m2)
    real(kind=rb), intent(out) :: dflx(ncol,nlay+1)      ! Total sky longwave downward flux (W/m2)
    real(kind=rb), intent(out) :: hr(ncol,nlay)          ! Total sky longwave radiative heating rate (K/d)
    real(kind=rb), intent(out) :: uflxc(ncol,nlay+1)     ! Clear sky longwave upward flux (W/m2)
    real(kind=rb), intent(out) :: dflxc(ncol,nlay+1)     ! Clear sky longwave downward flux (W/m2)
    real(kind=rb), intent(out) :: hrc(ncol,nlay)         ! Clear sky longwave radiative heating rate (K/d)

    real(kind=rb), intent(out) :: duflx_dt(ncol,nlay+1)  ! change in upward longwave flux (w/m2/K)
                                                         ! with respect to surface temperature
    real(kind=rb), intent(out) :: duflxc_dt(ncol,nlay+1) ! change in clear sky upward longwave flux (w/m2/K)
                                                         ! with respect to surface temperature
    real(kind=rb), intent(out) :: uflxspec(ncol,nlay+1,nbndlw) ! Total sky longwave upward flux spectrum (W/m2)
    real(kind=rb), intent(out) :: dflxspec(ncol,nlay+1,nbndlw) ! Total sky longwave downward flux spectrum (W/m2)
    real(kind=rb), intent(out) :: uflxcspec(ncol,nlay+1,nbndlw) ! Clear sky longwave upward flux spectrum (W/m2)
    real(kind=rb), intent(out) :: dflxcspec(ncol,nlay+1,nbndlw) ! Clear sky longwave downward flux spectrum (W/m2)

!  These are not comments! Necessary directives to f2py to handle array dimensions
!f2py depend(ncol,nlay) play, plev, tlay, tlev
!f2py depend(ncol,nlay) h2ovmr,o3vmr,co2vmr,ch4vmr,n2ovmr,o2vmr
!f2py depend(ncol,nlay) cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr
!f2py depend(ncol) tsfc, emis
!f2py depend(ncol) olr_sr
!f2py depend(ncol,nlay) tauaer
!f2py depend(ncol,nlay) cldfmcl,ciwpmcl,clwpmcl,taucmcl
!f2py depend(ncol,nlay) reicmcl,relqmcl
!f2py depend(ncol,nlay) uflx,dflx,hr,uflxc,dflxc,hrc,duflx_dt,duflxc_dt
!f2py depend(ncol,nlay,nbndlw) uflxspec, dflxspec, uflxcspec, dflxcspec

    !  Call the RRTMG_LW driver to compute radiative fluxes
    call rrtmg_lw(ncol    ,nlay    ,icld    ,ispec           ,idrv  , &
             play    , plev    , tlay    , tlev    , tsfc    , &
             h2ovmr  , o3vmr   , co2vmr  , ch4vmr  , n2ovmr  , o2vmr , &
             cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr , emis    , &
             inflglw , iceflglw, liqflglw, cldfmcl , &
             taucmcl , ciwpmcl , clwpmcl , reicmcl , relqmcl , tauaer  , &
             olr_sr  , uflx, dflx, hr, uflxc, dflxc,  hrc, &
             duflx_dt,duflxc_dt, uflxspec, dflxspec, uflxcspec, dflxcspec)

end subroutine climlab_rrtmg_lw_expanded

subroutine climlab_rrtmg_lw_ensemble &
    (ncol    ,nlay    , &
    permuteseed, irng, n_ens   ,col_by_col   ,do_seed_permutation , &
    icld     , ispec   , idrv    , &
    play    , plev    , tlay    , tlev    , tsfc    , &
    cldfrac, ciwp, clwp, reic, relq, tauc, &
    h2ovmr  , o3vmr   , co2vmr  , ch4vmr  , n2ovmr  , o2vmr , &
    cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr , emis    , &
    inflglw , iceflglw, liqflglw, tauaer  , &
    olr_sr  , uflx    , dflx    , hr      , uflxc   , dflxc,  hrc, &
    duflx_dt,duflxc_dt, uflxspec, dflxspec, uflxcspec, dflxcspec)

! Modules
    use parkind, only : im => kind_im
    !use parrrtm, only: nbndlw
    integer(kind=im), parameter :: nbndlw = 16

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
    integer(kind=im), intent(inout) :: ispec        ! spectral OLR output flag
    integer(kind=im), intent(in) :: idrv            ! Flag for calculation of dFdT, the change
                                                    !    in upward flux as a function of
                                                    !    surface temperature [0=off, 1=on]
                                                    !    0: Normal forward calculation
                                                    !    1: Normal forward calculation with
                                                    !       duflx_dt and duflxc_dt output
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
    real(kind=rb), intent(in) :: tauc(nbndlw,ncol,nlay)    ! in-cloud optical depth
    real(kind=rb), intent(in) :: h2ovmr(ncol,nlay)  ! H2O volume mixing ratio
    real(kind=rb), intent(in) :: o3vmr(ncol,nlay)   ! O3 volume mixing ratio
    real(kind=rb), intent(in) :: co2vmr(ncol,nlay)  ! CO2 volume mixing ratio
    real(kind=rb), intent(in) :: ch4vmr(ncol,nlay)  ! Methane volume mixing ratio
    real(kind=rb), intent(in) :: n2ovmr(ncol,nlay)  ! Nitrous oxide volume mixing ratio
    real(kind=rb), intent(in) :: o2vmr(ncol,nlay)   ! Oxygen volume mixing ratio
    real(kind=rb), intent(in) :: cfc11vmr(ncol,nlay)  ! CFC11 volume mixing ratio
    real(kind=rb), intent(in) :: cfc12vmr(ncol,nlay)  ! CFC12 volume mixing ratio
    real(kind=rb), intent(in) :: cfc22vmr(ncol,nlay)  ! CFC22 volume mixing ratio
    real(kind=rb), intent(in) :: ccl4vmr(ncol,nlay)   ! CCL4 volume mixing ratio
    real(kind=rb), intent(in) :: emis(ncol,nbndlw)  ! Surface emissivity
    integer(kind=im), intent(in) :: inflglw         ! Flag for cloud optical properties
    integer(kind=im), intent(in) :: iceflglw       ! Flag for ice particle specification
    integer(kind=im), intent(in) :: liqflglw        ! Flag for liquid droplet specification
    real(kind=rb), intent(in) :: tauaer(ncol,nlay,nbndlw) ! aerosol optical depth at mid-point of LW spectral bands

! Output
    real(kind=rb), intent(out) :: olr_sr(ncol,nbndlw)    ! Spectrally-decomposed OLR (W/m2)
    real(kind=rb), intent(out) :: uflx(ncol,nlay+1)      ! Total sky longwave upward flux (W/m2)
    real(kind=rb), intent(out) :: dflx(ncol,nlay+1)      ! Total sky longwave downward flux (W/m2)
    real(kind=rb), intent(out) :: hr(ncol,nlay)          ! Total sky longwave radiative heating rate (K/d)
    real(kind=rb), intent(out) :: uflxc(ncol,nlay+1)     ! Clear sky longwave upward flux (W/m2)
    real(kind=rb), intent(out) :: dflxc(ncol,nlay+1)     ! Clear sky longwave downward flux (W/m2)
    real(kind=rb), intent(out) :: hrc(ncol,nlay)         ! Clear sky longwave radiative heating rate (K/d)

    real(kind=rb), intent(out) :: duflx_dt(ncol,nlay+1)  ! change in upward longwave flux (w/m2/K)
                                                         ! with respect to surface temperature
    real(kind=rb), intent(out) :: duflxc_dt(ncol,nlay+1) ! change in clear sky upward longwave flux (w/m2/K)
                                                         ! with respect to surface temperature
    real(kind=rb), intent(out) :: uflxspec(ncol,nlay+1,nbndlw) ! Total sky longwave upward flux spectrum (W/m2)
    real(kind=rb), intent(out) :: dflxspec(ncol,nlay+1,nbndlw) ! Total sky longwave downward flux spectrum (W/m2)
    real(kind=rb), intent(out) :: uflxcspec(ncol,nlay+1,nbndlw) ! Clear sky longwave upward flux spectrum (W/m2)
    real(kind=rb), intent(out) :: dflxcspec(ncol,nlay+1,nbndlw) ! Clear sky longwave downward flux spectrum (W/m2)
    integer(kind=im) :: ind_ens, n_rrtmg_repeat, icol, ilay

!  These are not comments! Necessary directives to f2py to handle array dimensions
!f2py depend(ncol,nlay) play, plev, tlay, tlev
!f2py depend(ncol,nlay) cldfrac,ciwp,clwp,reic,relq
!f2py depend(ncol,nlay) tauc
!f2py depend(ncol,nlay) h2ovmr,o3vmr,co2vmr,ch4vmr,n2ovmr,o2vmr
!f2py depend(ncol,nlay) cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr
!f2py depend(ncol) tsfc, emis
!f2py depend(ncol) olr_sr
!f2py depend(ncol,nlay) tauaer
!f2py depend(ncol,nlay) uflx,dflx,hr,uflxc,dflxc,hrc,duflx_dt,duflxc_dt
!f2py depend(ncol,nlay,nbndlw) uflxspec, dflxspec, uflxcspec, dflxcspec

    if ((icld > 0).and.(do_seed_permutation.eq.1)) then
        n_rrtmg_repeat = n_ens
    else
        n_rrtmg_repeat = 1
    end if
    ! initialize output vars
    olr_sr = 0._rb
    uflx = 0._rb
    dflx = 0._rb
    hr = 0._rb
    uflxc = 0._rb
    dflxc = 0._rb
    hrc = 0._rb
    duflx_dt = 0._rb
    duflxc_dt = 0._rb
    uflxspec = 0._rb
    dflxspec = 0._rb
    uflxcspec = 0._rb
    dflxcspec = 0._rb
    ! These are not comments, they are using for multi-processing
    !$omp parallel do collapse(2) default(shared) &
    !$omp private(i, ilay) schedule(static)
    do icol=1, ncol
        do ind_ens=1, n_rrtmg_repeat
            call climlab_rrtmg_lw_single_col &
                (icol, ind_ens, ncol, nlay, &
                permuteseed, irng, col_by_col   ,do_seed_permutation , &
                icld     , ispec   , idrv    , &
                play    , plev    , tlay    , tlev    , tsfc    , &
                cldfrac, ciwp, clwp, reic, relq, tauc, &
                h2ovmr  , o3vmr   , co2vmr  , ch4vmr  , n2ovmr  , o2vmr , &
                cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr , emis    , &
                inflglw , iceflglw, liqflglw, tauaer  , &
                olr_sr  , uflx    , dflx    , hr      , uflxc   , dflxc,  hrc, &
                duflx_dt,duflxc_dt, uflxspec, dflxspec, uflxcspec, dflxcspec)
        end do
    end do
    !$omp end parallel do
    ! These are not comments, they are using for multi-processing

    olr_sr = olr_sr / n_rrtmg_repeat
    uflx = uflx / n_rrtmg_repeat
    dflx = dflx / n_rrtmg_repeat
    hr = hr / n_rrtmg_repeat
    uflxc = uflxc / n_rrtmg_repeat
    dflxc = dflxc / n_rrtmg_repeat
    hrc = hrc / n_rrtmg_repeat
    duflx_dt = duflx_dt / n_rrtmg_repeat
    duflxc_dt = duflxc_dt / n_rrtmg_repeat
    uflxspec = uflxspec / n_rrtmg_repeat
    dflxspec = dflxspec / n_rrtmg_repeat
    uflxcspec = uflxcspec / n_rrtmg_repeat
    dflxcspec = dflxcspec / n_rrtmg_repeat
end subroutine climlab_rrtmg_lw_ensemble

subroutine climlab_rrtmg_lw_single_col &
    (icol, ind_ens, ncol,    nlay    , &
    permuteseed, irng, col_by_col   ,do_seed_permutation , &
    icld     , ispec   , idrv    , &
    play    , plev    , tlay    , tlev    , tsfc    , &
    cldfrac, ciwp, clwp, reic, relq, tauc, &
    h2ovmr  , o3vmr   , co2vmr  , ch4vmr  , n2ovmr  , o2vmr , &
    cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr , emis    , &
    inflglw , iceflglw, liqflglw, tauaer  , &
    olr_sr  , uflx    , dflx    , hr      , uflxc   , dflxc,  hrc, &
    duflx_dt,duflxc_dt, uflxspec, dflxspec, uflxcspec, dflxcspec)

! Modules
    use parkind, only : im => kind_im
    use rrtmg_lw_rad, only: rrtmg_lw
    use mcica_subcol_gen_lw, only: mcica_subcol_lw
    !use parrrtm, only: nbndlw, ngptlw
    integer(kind=im), parameter :: nbndlw = 16
    integer(kind=im), parameter :: ngptlw = 140  

! Input
    integer, parameter :: rb = selected_real_kind(12)
    integer(kind=im), intent(in) :: icol            ! index for calculation
    integer(kind=im), intent(in) :: ind_ens         ! index for ensamble sampling
    integer(kind=im), intent(in) :: ncol            ! number of columns
    integer(kind=im), intent(in) :: nlay            ! number of model layers
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
    integer(kind=im), intent(inout) :: ispec        ! spectral OLR output flag
    integer(kind=im), intent(in) :: idrv            ! Flag for calculation of dFdT, the change
                                                    !    in upward flux as a function of
                                                    !    surface temperature [0=off, 1=on]
                                                    !    0: Normal forward calculation
                                                    !    1: Normal forward calculation with
                                                    !       duflx_dt and duflxc_dt output
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
    real(kind=rb), intent(in) :: tauc(nbndlw,ncol,nlay)    ! in-cloud optical depth
    real(kind=rb), intent(in) :: h2ovmr(ncol,nlay)  ! H2O volume mixing ratio
    real(kind=rb), intent(in) :: o3vmr(ncol,nlay)   ! O3 volume mixing ratio
    real(kind=rb), intent(in) :: co2vmr(ncol,nlay)  ! CO2 volume mixing ratio
    real(kind=rb), intent(in) :: ch4vmr(ncol,nlay)  ! Methane volume mixing ratio
    real(kind=rb), intent(in) :: n2ovmr(ncol,nlay)  ! Nitrous oxide volume mixing ratio
    real(kind=rb), intent(in) :: o2vmr(ncol,nlay)   ! Oxygen volume mixing ratio
    real(kind=rb), intent(in) :: cfc11vmr(ncol,nlay)  ! CFC11 volume mixing ratio
    real(kind=rb), intent(in) :: cfc12vmr(ncol,nlay)  ! CFC12 volume mixing ratio
    real(kind=rb), intent(in) :: cfc22vmr(ncol,nlay)  ! CFC22 volume mixing ratio
    real(kind=rb), intent(in) :: ccl4vmr(ncol,nlay)   ! CCL4 volume mixing ratio
    real(kind=rb), intent(in) :: emis(ncol,nbndlw)  ! Surface emissivity
    integer(kind=im), intent(in) :: inflglw         ! Flag for cloud optical properties
    integer(kind=im), intent(in) :: iceflglw       ! Flag for ice particle specification
    integer(kind=im), intent(in) :: liqflglw        ! Flag for liquid droplet specification
    real(kind=rb), intent(in) :: tauaer(ncol,nlay,nbndlw) ! aerosol optical depth at mid-point of LW spectral bands

! Output
    real(kind=rb), intent(inout) :: olr_sr(ncol,nbndlw)    ! Spectrally-decomposed OLR (W/m2)
    real(kind=rb), intent(inout) :: uflx(ncol,nlay+1)      ! Total sky longwave upward flux (W/m2)
    real(kind=rb), intent(inout) :: dflx(ncol,nlay+1)      ! Total sky longwave downward flux (W/m2)
    real(kind=rb), intent(inout) :: hr(ncol,nlay)          ! Total sky longwave radiative heating rate (K/d)
    real(kind=rb), intent(inout) :: uflxc(ncol,nlay+1)     ! Clear sky longwave upward flux (W/m2)
    real(kind=rb), intent(inout) :: dflxc(ncol,nlay+1)     ! Clear sky longwave downward flux (W/m2)
    real(kind=rb), intent(inout) :: hrc(ncol,nlay)         ! Clear sky longwave radiative heating rate (K/d)

    real(kind=rb), intent(inout) :: duflx_dt(ncol,nlay+1)  ! change in upward longwave flux (w/m2/K)
                                                         ! with respect to surface temperature
    real(kind=rb), intent(inout) :: duflxc_dt(ncol,nlay+1) ! change in clear sky upward longwave flux (w/m2/K)
                                                         ! with respect to surface temperature
    real(kind=rb), intent(inout) :: uflxspec(ncol,nlay+1,nbndlw) ! Total sky longwave upward flux spectrum (W/m2)
    real(kind=rb), intent(inout) :: dflxspec(ncol,nlay+1,nbndlw) ! Total sky longwave downward flux spectrum (W/m2)
    real(kind=rb), intent(inout) :: uflxcspec(ncol,nlay+1,nbndlw) ! Clear sky longwave upward flux spectrum (W/m2)
    real(kind=rb), intent(inout) :: dflxcspec(ncol,nlay+1,nbndlw) ! Clear sky longwave downward flux spectrum (W/m2)

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
    real(kind=rb) :: tauc1(nbndlw,1,nlay)    ! in-cloud optical depth
    real(kind=rb) :: h2ovmr1(1,nlay)  ! H2O volume mixing ratio
    real(kind=rb) :: o3vmr1(1,nlay)   ! O3 volume mixing ratio
    real(kind=rb) :: co2vmr1(1,nlay)  ! CO2 volume mixing ratio
    real(kind=rb) :: ch4vmr1(1,nlay)  ! Methane volume mixing ratio
    real(kind=rb) :: n2ovmr1(1,nlay)  ! Nitrous oxide volume mixing ratio
    real(kind=rb) :: o2vmr1(1,nlay)   ! Oxygen volume mixing ratio
    real(kind=rb) :: cfc11vmr1(1,nlay)  ! CFC11 volume mixing ratio
    real(kind=rb) :: cfc12vmr1(1,nlay)  ! CFC12 volume mixing ratio
    real(kind=rb) :: cfc22vmr1(1,nlay)  ! CFC22 volume mixing ratio
    real(kind=rb) :: ccl4vmr1(1,nlay)   ! CCL4 volume mixing ratio
    real(kind=rb) :: emis1(1,nbndlw)  ! Surface emissivity
    real(kind=rb) :: tauaer1(1,nlay,nbndlw) ! aerosol optical depth at mid-point of LW spectral bands

    !   These quantities are computed by McICA
    real(kind=rb) :: cldfmcl1(ngptlw,1,nlay)   ! cloud fraction [mcica]
    real(kind=rb) :: ciwpmcl1(ngptlw,1,nlay)   ! in-cloud ice water path [mcica]
    real(kind=rb) :: clwpmcl1(ngptlw,1,nlay)   ! in-cloud liquid water path [mcica]
    real(kind=rb) :: reicmcl1(1,nlay)          ! ice partcle size (microns)
    real(kind=rb) :: relqmcl1(1,nlay)          ! liquid particle size (microns)
    real(kind=rb) :: taucmcl1(ngptlw,1,nlay)   ! in-cloud optical depth [mcica]

    real(kind=rb) :: olr_sr1(1,nbndlw)    ! Spectrally-decomposed OLR (W/m2)
    real(kind=rb) :: uflx1(1,nlay+1)      ! Total sky longwave upward flux (W/m2)
    real(kind=rb) :: dflx1(1,nlay+1)      ! Total sky longwave downward flux (W/m2)
    real(kind=rb) :: hr1(1,nlay)          ! Total sky longwave radiative heating rate (K/d)
    real(kind=rb) :: uflxc1(1,nlay+1)     ! Clear sky longwave upward flux (W/m2)
    real(kind=rb) :: dflxc1(1,nlay+1)     ! Clear sky longwave downward flux (W/m2)
    real(kind=rb) :: hrc1(1,nlay)         ! Clear sky longwave radiative heating rate (K/d)
    real(kind=rb) :: duflx_dt1(1,nlay+1)  ! change in upward longwave flux (w/m2/K)
                                          ! with respect to surface temperature
    real(kind=rb) :: duflxc_dt1(1,nlay+1) ! change in clear sky upward longwave flux (w/m2/K)
                                                         ! with respect to surface temperature
    real(kind=rb) :: uflxspec1(1,nlay+1,nbndlw) ! Total sky longwave upward flux spectrum (W/m2)
    real(kind=rb) :: dflxspec1(1,nlay+1,nbndlw) ! Total sky longwave downward flux spectrum (W/m2)
    real(kind=rb) :: uflxcspec1(1,nlay+1,nbndlw) ! Clear sky longwave upward flux spectrum (W/m2)
    real(kind=rb) :: dflxcspec1(1,nlay+1,nbndlw) ! Clear sky longwave downward flux spectrum (W/m2)
    integer(kind=im) :: seed, ilay, shift

    ! Arrays initialization
    olr_sr1 = 0._rb
    uflx1 = 0._rb
    dflx1 = 0._rb
    hr1 = 0._rb
    uflxc1 = 0._rb
    dflxc1 = 0._rb
    hrc1 = 0._rb
    duflx_dt1 = 0._rb
    duflxc_dt1 = 0._rb
    uflxspec1 = 0._rb
    dflxspec1 = 0._rb
    uflxcspec1 = 0._rb
    dflxcspec1 = 0._rb

    cldfmcl1 = 0._rb
    ciwpmcl1 = 0._rb
    clwpmcl1 = 0._rb
    reicmcl1 = 0._rb
    relqmcl1 = 0._rb
    taucmcl1 = 0._rb

    !  Call the Monte Carlo Independent Column Approximation (McICA, Pincus et al., JC, 2003)
    if (do_seed_permutation.eq.1) then
        if (col_by_col.eq.0) then
            shift = ind_ens * ncol + icol - 1
        else 
            shift = ind_ens
        end if
        seed = permuteseed + shift * ngptlw
    else
        seed = permuteseed
    end if
    play1(1,:) = play(icol,:)
    plev1(1,:) = plev(icol,:)
    tlay1(1,:) = tlay(icol,:)
    tlev1(1,:) = tlev(icol,:)
    cldfrac1(1,:) = cldfrac(icol,:)
    ciwp1(1,:) = ciwp(icol,:)
    clwp1(1,:) = clwp(icol,:)
    reic1(1,:) = reic(icol,:)
    relq1(1,:) = relq(icol,:)
    tsfc1(1) = tsfc(icol)
    h2ovmr1(1,:) = h2ovmr(icol,:)
    o3vmr1(1,:) = o3vmr(icol,:)
    co2vmr1(1,:) = co2vmr(icol,:)
    ch4vmr1(1,:) = ch4vmr(icol,:)
    n2ovmr1(1,:) = n2ovmr(icol,:)
    o2vmr1(1,:) = o2vmr(icol,:)
    cfc11vmr1(1,:) = cfc11vmr(icol,:)
    cfc12vmr1(1,:) = cfc12vmr(icol,:)
    cfc22vmr1(1,:) = cfc22vmr(icol,:)
    ccl4vmr1(1,:) = ccl4vmr(icol,:)
    emis1(1,:) = emis(icol,:)
    do ilay=1, nlay
        tauc1(:,1,ilay) = tauc(:,icol,ilay)
        tauaer1(1,ilay,:) = tauaer(icol,ilay,:) 
    end do
    call mcica_subcol_lw(1, 1, nlay, icld, seed, irng, play1, &
            cldfrac1, ciwp1, clwp1, reic1, relq1, tauc1, cldfmcl1, &
            ciwpmcl1, clwpmcl1, reicmcl1, relqmcl1, taucmcl1)
    !  Call the RRTMG_LW driver to compute radiative fluxes
    call rrtmg_lw(1    ,nlay    ,icld    ,ispec           ,idrv  , &
            play1    , plev1    , tlay1    , tlev1    , tsfc1    , &
            h2ovmr1  , o3vmr1   , co2vmr1 , ch4vmr1  , n2ovmr1  , o2vmr1 , &
            cfc11vmr1, cfc12vmr1, cfc22vmr1, ccl4vmr1 , emis1    , &
            inflglw , iceflglw, liqflglw, cldfmcl1 , &
            taucmcl1 , ciwpmcl1 , clwpmcl1 , reicmcl1 , relqmcl1 , tauaer1  , &
            olr_sr1  , uflx1, dflx1, hr1, uflxc1, dflxc1,  hrc1, &
            duflx_dt1, duflxc_dt1, uflxspec1, dflxspec1, uflxcspec1, dflxcspec1)
    olr_sr(icol,:) = olr_sr(icol,:) + olr_sr1(1,:)
    uflx(icol,:) = uflx(icol,:) + uflx1(1,:)
    dflx(icol,:) = dflx(icol,:) + dflx1(1,:)
    hr(icol,:) = hr(icol,:) + hr1(1,:)
    uflxc(icol,:) = uflxc(icol,:) + uflxc1(1,:)
    dflxc(icol,:) = dflxc(icol,:) + dflxc1(1,:)
    hrc(icol,:) = hrc(icol,:) + hrc1(1,:)
    duflx_dt(icol,:) = duflx_dt(icol,:) + duflx_dt1(1,:)
    duflxc_dt(icol,:) = duflxc_dt(icol,:) + duflxc_dt1(1,:)
    uflxspec(icol,:,:) = uflxspec(icol,:,:) + uflxspec1(1,:,:)
    dflxspec(icol,:,:) = dflxspec(icol,:,:) + dflxspec1(1,:,:)
    uflxcspec(icol,:,:) = uflxcspec(icol,:,:) + uflxcspec1(1,:,:)
    dflxcspec(icol,:,:) = dflxcspec(icol,:,:) + dflxcspec1(1,:,:)
end subroutine climlab_rrtmg_lw_single_col