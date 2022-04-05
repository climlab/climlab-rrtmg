import pytest
from climlab_rrtmg._rrtmg_lw import _rrtmg_lw
from climlab_rrtmg._rrtmg_sw import _rrtmg_sw

cp = 1004.

def test_rrtmg_lw():
    nbndlw = int(_rrtmg_lw.parrrtm.nbndlw)
    ngptlw = int(_rrtmg_lw.parrrtm.ngptlw)
    #  Initialize absorption data
    _rrtmg_lw.climlab_rrtmg_lw_ini(cp)


def test_rrtmg_sw():
    nbndsw = int(_rrtmg_sw.parrrsw.nbndsw)
    naerec = int(_rrtmg_sw.parrrsw.naerec)
    ngptsw = int(_rrtmg_sw.parrrsw.ngptsw)
    #  Initialize absorption data
    _rrtmg_sw.climlab_rrtmg_sw_ini(cp)
