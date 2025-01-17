import pytest
import numpy as np
from climlab_rrtmg import rrtmg_lw, rrtmg_sw

# Specific heat at constant pressure
cp = 1004.

# Set up the pressure domain
ps = 1000. #  Surface pressure in hPa
# RRTM code expects arrays with (ncol, nlay)
# and with pressure decreasing from surface at element 0
ncol = 1
nlay = 30
deltap = ps / nlay  # pressure interval
plev = np.linspace(ps, 0., nlay+1)  # pressure bounds
plev = plev[np.newaxis, ...]
play = np.linspace(ps-deltap/2., deltap/2., nlay)
play = play[np.newaxis, ...]
# Set the temperatures
#  Using a linearly decreasing temperature from surface to TOA
tsfc = 288.
tlay = np.linspace(288.-10., 200., nlay)
tlay = tlay[np.newaxis, ...]
tlev = np.linspace(288., 200., nlay+1)
tlev = tlev[np.newaxis, ...]

# atmospheric composition
# specific humidity profile computed from those temperatures using Manabe's
# fixed relative humidity profile
specific_humidity = np.array([4.13141097e-03, 3.41509495e-03, 2.81099479e-03, 2.30359570e-03,
       1.87920573e-03, 1.52578624e-03, 1.23279279e-03, 9.91026303e-04,
       7.92494475e-04, 6.30283118e-04, 4.98437246e-04, 3.91851633e-04,
       3.06170488e-04, 2.37695932e-04, 1.83304857e-04, 1.40373783e-04,
       1.06711275e-04, 8.04974602e-05, 6.02302082e-05, 4.46774859e-05,
       3.28354282e-05, 2.38916443e-05, 1.71932832e-05, 1.22193649e-05,
       8.55682965e-06, 5.87957411e-06, 5.00000000e-06, 5.00000000e-06,
       5.00000000e-06, 5.00000000e-06])
# Convert to volume mixing ratio from mass mixing ratio
#  just multiplying by ratio of molecular weights of dry air and H2O
h2ovmr = specific_humidity * 28.97 / 18.01528
h2ovmr = h2ovmr[np.newaxis, ...]
#  A global-mean ozone climatology
o3vmr = np.array([2.25573888e-08, 2.38730436e-08, 2.52586476e-08, 2.66442517e-08,
       2.80298557e-08, 2.97254145e-08, 3.14254923e-08, 3.31238355e-08,
       3.46767916e-08, 3.62297478e-08, 3.76122833e-08, 3.86410454e-08,
       3.96698075e-08, 4.08899052e-08, 4.21303310e-08, 4.39781220e-08,
       4.60528063e-08, 4.87636254e-08, 5.16974065e-08, 5.57122567e-08,
       6.17914190e-08, 7.15771368e-08, 9.29020109e-08, 1.29109217e-07,
       1.75914529e-07, 2.45552383e-07, 3.92764464e-07, 7.61726407e-07,
       2.25137178e-06, 7.27500161e-06])
o3vmr = o3vmr[np.newaxis, ...]
# Other values taken from the AquaPlanet Experiment protocols,
# except for O2 which is set the realistic value 0.21
co2vmr = 348. / 1E6 * np.ones_like(play)
ch4vmr = 1650. / 1E9 * np.ones_like(play)
n2ovmr = 306. / 1E9 * np.ones_like(play)
o2vmr = 0.21 * np.ones_like(play)
cfc11vmr = 0. * np.ones_like(play)
cfc12vmr = 0. * np.ones_like(play)
cfc22vmr = 0. * np.ones_like(play)
ccl4vmr = 0. * np.ones_like(play)

#  Cloud parameters
cloud_level_index = 8
cldfrac = 0.5*np.exp(-(play-play[0,cloud_level_index])**2/(2*25.)**2)  # Layer cloud fraction: a Gaussian centered on a pressure level
clwp = 60. * np.ones_like(play)  # in-cloud liquid water path (g/m2)
ciwp = 0. * np.ones_like(play)   # in-cloud ice water path (g/m2)
relq = 14. * np.ones_like(play) # Cloud water drop effective radius (microns)
reic = 0. * np.ones_like(play)  # Cloud ice particle effective size (microns)


def test_rrtmg_lw_clearsky():
    nbndlw = int(rrtmg_lw.parrrtm.nbndlw)
    ngptlw = int(rrtmg_lw.parrrtm.ngptlw)
    #  Initialize absorption data
    rrtmg_lw.climlab_rrtmg_lw_ini(cp)

    # Lots of RRTMG parameters
    icld = 0    # Cloud overlap method, 0: Clear only, 1: Random, 2,  Maximum/random] 3: Maximum
    irng = 1  # more monte carlo stuff
    idrv = 0  # whether to also calculate the derivative of flux with respect to surface temp
    permuteseed = 300
    inflglw  = 2
    iceflglw = 1
    liqflglw = 1

    #  These arrays have an extra dimension for number of bands
    # in-cloud optical depth [nbndlw,ncol,nlay]
    tauc = 0. * np.ones_like(play)
    #  broadcast to get [nbndlw,ncol,nlay]
    tauc = tauc * np.ones([nbndlw,ncol,nlay])
    # Aerosol optical depth at mid-point of LW spectral bands [ncol,nlay,nbndlw]
    tauaer = 0. * np.ones_like(play)
    #  broadcast and transpose to get [ncol,nlay,nbndlw]
    tauaer = np.transpose(tauaer * np.ones([nbndlw,ncol,nlay]), (1,2,0))

    # surface emissivity
    emis = 1. * np.ones((ncol,nbndlw))

    # Clear-sky only
    cldfmcl = np.zeros((ngptlw,ncol,nlay))
    ciwpmcl = np.zeros((ngptlw,ncol,nlay))
    clwpmcl = np.zeros((ngptlw,ncol,nlay))
    reicmcl = np.zeros((ncol,nlay))
    relqmcl = np.zeros((ncol,nlay))
    taucmcl = np.zeros((ngptlw,ncol,nlay))

    for ispec in [0, 1]: # Spectral OLR output flag, 0: only calculate total fluxes, 1: also return spectral OLR
        # Call the RRTMG_LW driver
        (olr_sr, uflx, dflx, hr, uflxc, dflxc, hrc, duflx_dt, duflxc_dt) = \
                rrtmg_lw.climlab_rrtmg_lw(ncol, nlay, icld, ispec, idrv,
                     play, plev, tlay, tlev, tsfc,
                     h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
                     cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, emis,
                     inflglw, iceflglw, liqflglw, cldfmcl,
                     taucmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl,
                     tauaer)

def test_rrtmg_lw_mcica():
    nbndlw = int(rrtmg_lw.parrrtm.nbndlw)
    ngptlw = int(rrtmg_lw.parrrtm.ngptlw)
    #  Initialize absorption data
    rrtmg_lw.climlab_rrtmg_lw_ini(cp)

    # Lots of RRTMG parameters
    icld = 1    # Cloud overlap method, 0: Clear only, 1: Random, 2,  Maximum/random] 3: Maximum
    irng = 1  # more monte carlo stuff
    idrv = 0  # whether to also calculate the derivative of flux with respect to surface temp
    permuteseed = 300
    inflglw  = 2
    iceflglw = 1
    liqflglw = 1

    #  These arrays have an extra dimension for number of bands
    # in-cloud optical depth [nbndlw,ncol,nlay]
    tauc = 0. * np.ones_like(play)
    #  broadcast to get [nbndlw,ncol,nlay]
    tauc = tauc * np.ones([nbndlw,ncol,nlay])
    # Aerosol optical depth at mid-point of LW spectral bands [ncol,nlay,nbndlw]
    tauaer = 0. * np.ones_like(play)
    #  broadcast and transpose to get [ncol,nlay,nbndlw]
    tauaer = np.transpose(tauaer * np.ones([nbndlw,ncol,nlay]), (1,2,0))

    # surface emissivity
    emis = 1. * np.ones((ncol,nbndlw))

    #  Call the Monte Carlo Independent Column Approximation (McICA, Pincus et al., JC, 2003)
    (cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, taucmcl) = \
        rrtmg_lw.climlab_mcica_subcol_lw(
                    ncol, nlay, icld,
                    permuteseed, irng, play,
                    cldfrac, ciwp, clwp, reic, relq, tauc)

    for ispec in [0, 1]: # Spectral OLR output flag, 0: only calculate total fluxes, 1: also return spectral OLR
        # Call the RRTMG_LW driver
        (olr_sr, uflx, dflx, hr, uflxc, dflxc, hrc, duflx_dt, duflxc_dt) = \
                rrtmg_lw.climlab_rrtmg_lw(ncol, nlay, icld, ispec, idrv,
                     play, plev, tlay, tlev, tsfc,
                     h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
                     cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, emis,
                     inflglw, iceflglw, liqflglw, cldfmcl,
                     taucmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl,
                     tauaer)


def test_rrtmg_sw_clearsky():
    nbndsw = int(rrtmg_sw.parrrsw.nbndsw)
    naerec = int(rrtmg_sw.parrrsw.naerec)
    ngptsw = int(rrtmg_sw.parrrsw.ngptsw)
    #  Initialize absorption data
    rrtmg_sw.climlab_rrtmg_sw_ini(cp)
    #  Lots of RRTMG parameters
    icld = 0    # Cloud overlap method, 0: Clear only, 1: Random, 2,  Maximum/random] 3: Maximum
    irng = 1  # more monte carlo stuff
    permuteseed = 150
    dyofyr = 0       # day of the year used to get Earth/Sun distance (if not adjes)
    inflgsw  = 2
    iceflgsw = 1
    liqflgsw = 1
    # AEROSOLS
    iaer = 0   #! Aerosol option flag
                #!    0: No aerosol
                #!    6: ECMWF method: use six ECMWF aerosol types input aerosol optical depth at 0.55 microns for each aerosol type (ecaer)
                #!    10:Input aerosol optical properties: input total aerosol optical depth, single scattering albedo and asymmetry parameter (tauaer, ssaaer, asmaer) directly
    tauaer = 0. * np.ones_like(play)   # Aerosol optical depth (iaer=10 only), Dimensions,  (ncol,nlay,nbndsw)] #  (non-delta scaled)
    #  broadcast and transpose to get [ncol,nlay,nbndsw]
    tauaer = np.transpose(tauaer * np.ones([nbndsw,ncol,nlay]), (1,2,0))
    ssaaer = 0. * np.ones_like(play)   # Aerosol single scattering albedo (iaer=10 only), Dimensions,  (ncol,nlay,nbndsw)] #  (non-delta scaled)
    #  broadcast and transpose to get [ncol,nlay,nbndsw]
    ssaaer = np.transpose(ssaaer * np.ones([nbndsw,ncol,nlay]), (1,2,0))
    asmaer = 0. * np.ones_like(play)   # Aerosol asymmetry parameter (iaer=10 only), Dimensions,  (ncol,nlay,nbndsw)] #  (non-delta scaled)
    #  broadcast and transpose to get [ncol,nlay,nbndsw]
    asmaer = np.transpose(asmaer * np.ones([nbndsw,ncol,nlay]), (1,2,0))
    ecaer  = 0. * np.ones_like(play)   # Aerosol optical depth at 0.55 micron (iaer=6 only), Dimensions,  (ncol,nlay,naerec)] #  (non-delta scaled)
    #  broadcast and transpose to get [ncol,nlay,naerec]
    ecaer = np.transpose(ecaer * np.ones([naerec,ncol,nlay]), (1,2,0))
    # insolation
    scon = 1365.2  # solar constant
    coszen = 1/4  # cosine of zenith angle
    adjes = 1.  # instantaneous irradiance = scon * eccentricity_factor
    dyofyr = 0       # day of the year used to get Earth/Sun distance (if not adjes)
    # new arguments for RRTMG_SW version 4.0
    isolvar = -1    # ! Flag for solar variability method
                    # !   -1 = (when scon .eq. 0.0): No solar variability
                    # !        and no solar cycle (Kurucz solar irradiance
                    # !        of 1368.22 Wm-2 only);
                    # !        (when scon .ne. 0.0): Kurucz solar irradiance
                    # !        scaled to scon and solar variability defined
                    # !        (optional) by setting non-zero scale factors
                    # !        for each band in bndsolvar
                    # !    0 = (when SCON .eq. 0.0): No solar variability
                    # !        and no solar cycle (NRLSSI2 solar constant of
                    # !        1360.85 Wm-2 for the 100-50000 cm-1 spectral
                    # !        range only), with facular and sunspot effects
                    # !        fixed to the mean of Solar Cycles 13-24;
                    # !        (when SCON .ne. 0.0): No solar variability
                    # !        and no solar cycle (NRLSSI2 solar constant of
                    # !        1360.85 Wm-2 for the 100-50000 cm-1 spectral
                    # !        range only), is scaled to SCON
                    # !    1 = Solar variability (using NRLSSI2  solar
                    # !        model) with solar cycle contribution
                    # !        determined by fraction of solar cycle
                    # !        with facular and sunspot variations
                    # !        fixed to their mean variations over the
                    # !        average of Solar Cycles 13-24;
                    # !        two amplitude scale factors allow
                    # !        facular and sunspot adjustments from
                    # !        mean solar cycle as defined by indsolvar
                    # !    2 = Solar variability (using NRLSSI2 solar
                    # !        model) over solar cycle determined by
                    # !        direct specification of Mg (facular)
                    # !        and SB (sunspot) indices provided
                    # !        in indsolvar (scon = 0.0 only)
                    # !    3 = (when scon .eq. 0.0): No solar variability
                    # !        and no solar cycle (NRLSSI2 solar irradiance
                    # !        of 1360.85 Wm-2 only);
                    # !        (when scon .ne. 0.0): NRLSSI2 solar irradiance
                    # !        scaled to scon and solar variability defined
                    # !        (optional) by setting non-zero scale factors
                    # !        for each band in bndsolvar
    indsolvar = np.ones(2)      # Facular and sunspot amplitude scale factors (isolvar=1),
                                 # or Mg and SB indices (isolvar=2)
    bndsolvar = np.ones(nbndsw) # Solar variability scale factors for each shortwave band
    solcycfrac = 1.              # Fraction of averaged solar cycle (0-1) at current time (isolvar=1)

    # surface albedo
    aldif = 0.3
    aldir = 0.3
    asdif = 0.3
    asdir = 0.3

    # Clear-sky only
    cldfmcl = np.zeros((ngptsw,ncol,nlay))
    ciwpmcl = np.zeros((ngptsw,ncol,nlay))
    clwpmcl = np.zeros((ngptsw,ncol,nlay))
    reicmcl = np.zeros((ncol,nlay))
    relqmcl = np.zeros((ncol,nlay))
    taucmcl = np.zeros((ngptsw,ncol,nlay))
    ssacmcl = np.zeros((ngptsw,ncol,nlay))
    asmcmcl = np.zeros((ngptsw,ncol,nlay))
    fsfcmcl = np.zeros((ngptsw,ncol,nlay))

    (swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc) = \
            rrtmg_sw.climlab_rrtmg_sw(ncol, nlay, icld, iaer,
                play, plev, tlay, tlev, tsfc,
                h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
                asdir, asdif, aldir, aldif,
                coszen, adjes, dyofyr, scon, isolvar,
                inflgsw, iceflgsw, liqflgsw, cldfmcl,
                taucmcl, ssacmcl, asmcmcl, fsfcmcl,
                ciwpmcl, clwpmcl, reicmcl, relqmcl,
                tauaer, ssaaer, asmaer, ecaer,
                bndsolvar, indsolvar, solcycfrac)

def test_rrtmg_sw_mcica():
    nbndsw = int(rrtmg_sw.parrrsw.nbndsw)
    naerec = int(rrtmg_sw.parrrsw.naerec)
    ngptsw = int(rrtmg_sw.parrrsw.ngptsw)
    #  Initialize absorption data
    rrtmg_sw.climlab_rrtmg_sw_ini(cp)
    #  Lots of RRTMG parameters
    icld = 1    # Cloud overlap method, 0: Clear only, 1: Random, 2,  Maximum/random] 3: Maximum
    irng = 1  # more monte carlo stuff
    permuteseed = 150
    dyofyr = 0       # day of the year used to get Earth/Sun distance (if not adjes)
    inflgsw  = 2
    iceflgsw = 1
    liqflgsw = 1
    tauc = 0.
    ssac = 0.  # In-cloud single scattering albedo
    asmc = 0.  # In-cloud asymmetry parameter
    fsfc = 0.  # In-cloud forward scattering fraction (delta function pointing forward "forward peaked scattering")
    # AEROSOLS
    iaer = 0   #! Aerosol option flag
                #!    0: No aerosol
                #!    6: ECMWF method: use six ECMWF aerosol types input aerosol optical depth at 0.55 microns for each aerosol type (ecaer)
                #!    10:Input aerosol optical properties: input total aerosol optical depth, single scattering albedo and asymmetry parameter (tauaer, ssaaer, asmaer) directly
    tauaer = 0. * np.ones_like(play)   # Aerosol optical depth (iaer=10 only), Dimensions,  (ncol,nlay,nbndsw)] #  (non-delta scaled)
    #  broadcast and transpose to get [ncol,nlay,nbndsw]
    tauaer = np.transpose(tauaer * np.ones([nbndsw,ncol,nlay]), (1,2,0))
    ssaaer = 0. * np.ones_like(play)   # Aerosol single scattering albedo (iaer=10 only), Dimensions,  (ncol,nlay,nbndsw)] #  (non-delta scaled)
    #  broadcast and transpose to get [ncol,nlay,nbndsw]
    ssaaer = np.transpose(ssaaer * np.ones([nbndsw,ncol,nlay]), (1,2,0))
    asmaer = 0. * np.ones_like(play)   # Aerosol asymmetry parameter (iaer=10 only), Dimensions,  (ncol,nlay,nbndsw)] #  (non-delta scaled)
    #  broadcast and transpose to get [ncol,nlay,nbndsw]
    asmaer = np.transpose(asmaer * np.ones([nbndsw,ncol,nlay]), (1,2,0))
    ecaer  = 0. * np.ones_like(play)   # Aerosol optical depth at 0.55 micron (iaer=6 only), Dimensions,  (ncol,nlay,naerec)] #  (non-delta scaled)
    #  broadcast and transpose to get [ncol,nlay,naerec]
    ecaer = np.transpose(ecaer * np.ones([naerec,ncol,nlay]), (1,2,0))
    #  These cloud arrays have an extra dimension for number of bands
    # in-cloud optical depth [nbndsw,ncol,nlay]
    tauc = 0. * np.ones_like(play)
    #  broadcast to get [nbndsw,ncol,nlay]
    tauc = tauc * np.ones([nbndsw,ncol,nlay])
    # In-cloud single scattering albedo, same operation
    ssac = 0. * np.ones_like(play) * np.ones([nbndsw,ncol,nlay])
    # In-cloud asymmetry parameter
    asmc = 0. * np.ones_like(play) * np.ones([nbndsw,ncol,nlay])
    # In-cloud forward scattering fraction (delta function pointing forward "forward peaked scattering")
    fsfc = 0. * np.ones_like(play) * np.ones([nbndsw,ncol,nlay])
    # insolation
    scon = 1365.2  # solar constant
    coszen = 1/4  # cosine of zenith angle
    adjes = 1.  # instantaneous irradiance = scon * eccentricity_factor
    dyofyr = 0       # day of the year used to get Earth/Sun distance (if not adjes)
    # new arguments for RRTMG_SW version 4.0
    isolvar = -1    # ! Flag for solar variability method
    indsolvar = np.ones(2)      # Facular and sunspot amplitude scale factors (isolvar=1),
                                 # or Mg and SB indices (isolvar=2)
    bndsolvar = np.ones(nbndsw) # Solar variability scale factors for each shortwave band
    solcycfrac = 1.              # Fraction of averaged solar cycle (0-1) at current time (isolvar=1)

    # surface albedo
    aldif = 0.3
    aldir = 0.3
    asdif = 0.3
    asdir = 0.3

    #  Call the Monte Carlo Independent Column Approximation (McICA, Pincus et al., JC, 2003)
    (cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, taucmcl,
    ssacmcl, asmcmcl, fsfcmcl) = rrtmg_sw.climlab_mcica_subcol_sw(
                    ncol, nlay, icld, permuteseed, irng, play,
                    cldfrac, ciwp, clwp, reic, relq, tauc, ssac, asmc, fsfc)

    (swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc) = \
            rrtmg_sw.climlab_rrtmg_sw(ncol, nlay, icld, iaer,
                play, plev, tlay, tlev, tsfc,
                h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
                asdir, asdif, aldir, aldif,
                coszen, adjes, dyofyr, scon, isolvar,
                inflgsw, iceflgsw, liqflgsw, cldfmcl,
                taucmcl, ssacmcl, asmcmcl, fsfcmcl,
                ciwpmcl, clwpmcl, reicmcl, relqmcl,
                tauaer, ssaaer, asmaer, ecaer,
                bndsolvar, indsolvar, solcycfrac)

def test_rrtmg_sw_multicol():
    ncol = 2
    plev_2d = np.tile(plev, [ncol, 1])
    play_2d = np.tile(play, [ncol, 1])
    tlay_2d = np.tile(tlay, [ncol, 1])
    tlev_2d = np.tile(tlev, [ncol, 1])
    tsfc_2d = np.tile(tsfc, [ncol])
    h2ovmr_2d = np.tile(h2ovmr, [ncol, 1])
    o3vmr_2d = np.tile(o3vmr, [ncol, 1])
    co2vmr_2d = np.tile(co2vmr, [ncol, 1])
    ch4vmr_2d = np.tile(ch4vmr, [ncol, 1])
    n2ovmr_2d = np.tile(n2ovmr, [ncol, 1])
    o2vmr_2d = np.tile(o2vmr, [ncol, 1])

    nbndsw = int(rrtmg_sw.parrrsw.nbndsw)
    naerec = int(rrtmg_sw.parrrsw.naerec)
    ngptsw = int(rrtmg_sw.parrrsw.ngptsw)
    #  Initialize absorption data
    rrtmg_sw.climlab_rrtmg_sw_ini(cp)
    #  Lots of RRTMG parameters
    icld = 0    # Cloud overlap method, 0: Clear only, 1: Random, 2,  Maximum/random] 3: Maximum
    dyofyr = 0       # day of the year used to get Earth/Sun distance (if not adjes)
    inflgsw  = 2
    iceflgsw = 1
    liqflgsw = 1
    # AEROSOLS
    iaer = 0   #! Aerosol option flag
                #!    0: No aerosol
                #!    6: ECMWF method: use six ECMWF aerosol types input aerosol optical depth at 0.55 microns for each aerosol type (ecaer)
                #!    10:Input aerosol optical properties: input total aerosol optical depth, single scattering albedo and asymmetry parameter (tauaer, ssaaer, asmaer) directly
    tauaer_2d = 0. * np.ones_like(play_2d)   # Aerosol optical depth (iaer=10 only), Dimensions,  (ncol,nlay,nbndsw)] #  (non-delta scaled)
    #  broadcast and transpose to get [ncol,nlay,nbndsw]
    tauaer_2d = np.transpose(tauaer_2d * np.ones([nbndsw,ncol,nlay]), (1,2,0))
    ssaaer_2d = 0. * np.ones_like(play)   # Aerosol single scattering albedo (iaer=10 only), Dimensions,  (ncol,nlay,nbndsw)] #  (non-delta scaled)
    #  broadcast and transpose to get [ncol,nlay,nbndsw]
    ssaaer_2d = np.transpose(ssaaer_2d * np.ones([nbndsw,ncol,nlay]), (1,2,0))
    asmaer_2d = 0. * np.ones_like(play_2d)   # Aerosol asymmetry parameter (iaer=10 only), Dimensions,  (ncol,nlay,nbndsw)] #  (non-delta scaled)
    #  broadcast and transpose to get [ncol,nlay,nbndsw]
    asmaer_2d = np.transpose(asmaer_2d * np.ones([nbndsw,ncol,nlay]), (1,2,0))
    ecaer_2d  = 0. * np.ones_like(play_2d)   # Aerosol optical depth at 0.55 micron (iaer=6 only), Dimensions,  (ncol,nlay,naerec)] #  (non-delta scaled)
    #  broadcast and transpose to get [ncol,nlay,naerec]
    ecaer_2d = np.transpose(ecaer_2d * np.ones([naerec,ncol,nlay]), (1,2,0))
    # insolation
    scon = 1365.2  # solar constant
    coszen_2d = np.tile(1/4, [ncol])  # cosine of zenith angle
    adjes = 1.  # instantaneous irradiance = scon * eccentricity_factor
    dyofyr = 0       # day of the year used to get Earth/Sun distance (if not adjes)
    # new arguments for RRTMG_SW version 4.0
    isolvar = -1    # ! Flag for solar variability method
    indsolvar = np.ones(2)      # Facular and sunspot amplitude scale factors (isolvar=1),
                                 # or Mg and SB indices (isolvar=2)
    bndsolvar = np.ones(nbndsw) # Solar variability scale factors for each shortwave band
    solcycfrac = 1.              # Fraction of averaged solar cycle (0-1) at current time (isolvar=1)

    # surface albedo
    aldif_2d = 0.3 * np.ones_like(tsfc_2d)
    aldir_2d = 0.3 * np.ones_like(tsfc_2d)
    asdif_2d = 0.3 * np.ones_like(tsfc_2d)
    asdir_2d = 0.3 * np.ones_like(tsfc_2d)

    # Clear-sky only
    cldfmcl_2d = np.zeros((ngptsw,ncol,nlay))
    ciwpmcl_2d = np.zeros((ngptsw,ncol,nlay))
    clwpmcl_2d = np.zeros((ngptsw,ncol,nlay))
    reicmcl_2d = np.zeros((ncol,nlay))
    relqmcl_2d = np.zeros((ncol,nlay))
    taucmcl_2d = np.zeros((ngptsw,ncol,nlay))
    ssacmcl_2d = np.zeros((ngptsw,ncol,nlay))
    asmcmcl_2d = np.zeros((ngptsw,ncol,nlay))
    fsfcmcl_2d = np.zeros((ngptsw,ncol,nlay))

    (swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc) = \
            rrtmg_sw.climlab_rrtmg_sw(ncol, nlay, icld, iaer,
                play_2d, plev_2d, tlay_2d, tlev_2d, tsfc_2d,
                h2ovmr_2d, o3vmr_2d, co2vmr_2d, ch4vmr_2d, n2ovmr_2d, o2vmr_2d,
                asdir_2d, asdif_2d, aldir_2d, aldif_2d,
                coszen_2d, adjes, dyofyr, scon, isolvar,
                inflgsw, iceflgsw, liqflgsw, cldfmcl_2d,
                taucmcl_2d, ssacmcl_2d, asmcmcl_2d, fsfcmcl_2d,
                ciwpmcl_2d, clwpmcl_2d, reicmcl_2d, relqmcl_2d,
                tauaer_2d, ssaaer_2d, asmaer_2d, ecaer_2d,
                bndsolvar, indsolvar, solcycfrac)
