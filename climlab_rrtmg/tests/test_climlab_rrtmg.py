import pytest
from climlab_rrtmg import rrtmg_lw
from climlab_rrtmg import rrtmg_sw

cp = 1004.

def test_rrtmg_lw():
    nbndlw = int(rrtmg_lw.parrrtm.nbndlw)
    ngptlw = int(rrtmg_lw.parrrtm.ngptlw)
    #  Initialize absorption data
    rrtmg_lw.climlab_rrtmg_lw_ini(cp)


def test_rrtmg_sw():
    nbndsw = int(rrtmg_sw.parrrsw.nbndsw)
    naerec = int(rrtmg_sw.parrrsw.naerec)
    ngptsw = int(rrtmg_sw.parrrsw.ngptsw)
    #  Initialize absorption data
    rrtmg_sw.climlab_rrtmg_sw_ini(cp)
