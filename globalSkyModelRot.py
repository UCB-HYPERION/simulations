import matplotlib.pylab as pl
import numpy as np
import aipy as a
import uvtools
import os
import globalSignalVis

import sys
sys.path.append('/home/kara/capo/kmk/scripts')
from GlobalSkyModel import GlobalSkyModel

def makeGSM(nside, freq, path='/home/kara/capo/kmk/gsm/gsm_raw/'):
        gsmMap = a.healpix.HealpixMap(nside=nside)
        g = a.healpix.HealpixMap(nside=512)
        g.map = GlobalSkyModel(freq=1000*freq, GSMlocation=path, GSMNSIDE=512).map # galactic coordinates
        gsmMap.from_hpm(g)
        return gsmMap

def hourAngle(sidereal_time, ra):
    return sidereal_time - ra

calfile = 'hyperion_deployment_aug2017'
freqs = np.array([0.100]) # GHz

# set array and system parameters
aa = a.cal.get_aa(calfile, freqs)
aa.set_jultime()
jd_init = aa.get_jultime()
print jd_init
jds = np.linspace(jd_init, jd_init + 0.5, 12)

# generate galactic sky
hpm = makeGSM(nside=64, freq=freqs[0])
sky = a.healpix.HealpixMap(nside=64)

for i in np.arange(len(jds)):
    aa.set_jultime(jds[i])
    gxyz = gx,gy,gz, = hpm.px2crd(np.arange(hpm.npix())) # galactic
    ga2eq = a.coord.convert_m(isys='eq', osys='ga', oepoch=aa.epoch)
    # conversion matrix so the isys and osys are reversed
    exyz = ex,ey,ez = np.dot(ga2eq, hpm.px2crd(np.arange(hpm.npix()))) # equatorial
    txyz = tx,ty,tz = np.dot(a.coord.eq2top_m(aa.sidereal_time(), aa.lat), exyz)
    az,alt = a.coord.top2azalt(txyz)
    sky.map = hpm[az,alt]
    uvtools.plot.plot_hmap_ortho(sky, res=1)
    pl.show()
