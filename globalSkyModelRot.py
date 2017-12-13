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
    filename = "gsm_%dMHz.fits" % int(1000*freq)
    if os.path.isfile(filename)==True:
        gsmMap.from_fits(filename)
    else:
        g = a.healpix.HealpixMap(nside=512)
        g.map = GlobalSkyModel(freq=1000*freq, GSMlocation=path, GSMNSIDE=512).map # galactic coordinates
        gsmMap.from_hpm(g)
        gsmMap.to_fits(filename)
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
jds = np.array([jd_init])
#jds = np.linspace(jd_init, jd_init + 0.5, 12)

# generate galactic sky
hpm = makeGSM(nside=64, freq=freqs[0])
sky = a.healpix.HealpixMap(nside=64)

for i in np.arange(len(jds)):
    aa.set_jultime(jds[i])
    txyz = tx,ty,tz, = sky.px2crd(np.arange(sky.npix())) # topocentric
    eq2ga = a.coord.convert_m(isys='eq', osys='ga', oepoch=aa.epoch)
    top2eq = np.linalg.inv(aa.eq2top_m)
    top2ga = np.dot(eq2ga, top2eq)
    gxyz = gx,gy,gz = np.dot(top2ga, txyz) # galactic
    sky.map = hpm[gx,gy,gz]
    uvtools.plot.plot_hmap_ortho(sky, res=1)
    pl.show()
