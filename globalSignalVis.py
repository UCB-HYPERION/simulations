#! /usr/bin/env python
import aipy as a
#import capo
import numpy as np
import os
#from GlobalSkyModel import GlobalSkyModel

class AntennaArray(a.pol.AntennaArray):
    def geometric_dly(self, i, j):
        """ converts unique baseline numbers into nanosecond delays """
        bx, by, bz = self.get_baseline(i,j)
        return np.sqrt(bx**2 + by**2 + bz**2)        

def calcFreq(array, ij, ref_ij, ref_freq, min_freq, max_freq):
    """ calculates nuCal frequencies (in GHz) based on an initial frequency, reference baseline, 
        and input baseline. must include min and max frequencies for array in order to
        ensure that calculated frequencies are legal. """

    bl = bl2delay(array, ij)
    ref_bl = bl2delay(array, ref_ij)
    f = ref_freq * (ref_bl / bl)
    assert(max_freq >= f >= min_freq)

    i, j = a.miriad.bl2ij(bl)
    bl = aa.get_baseline(i, j)
    bl = np.sqrt(np.dot(bl, bl))
    f = ref_freq * (ref_bl / baseline)
    assert(max_freq >= f >= min_freq) # check for invalid freq

    return f

def makeFlatMap(nside, freq, Tsky=1.0):
    """ fill sky map of given size and frequency with flat temperature across whole sky,
        returns map in Janskys. """
    hpm = a.healpix.HealpixMap(nside=nside)
    hpm.map = Tsky*np.ones(shape = hpm.map.shape)
    print "flat map size = " + str(hpm.map.shape)
    return hpm_TtoJy(hpm, freq)

def makeSynchMap(nside, freq=0.150):
    """ fill sky map of given size and frequency with flat temperature across whole sky based on synchrotron spectrum,
        returns map in Janskys. """
    hpm = a.healpix.HealpixMap(nside=nside)
    Tsky = 237 * (freq/0.150)**-2.5
    hpm.map = Tsky*np.ones(shape = hpm.map.shape)
    print "flat map size = " + str(hpm.map.shape)
    return hpm_TtoJy(hpm, freq)

def hpm_TtoJy(hpm, freq):
    """ converts given Healpix map from brightness temperature to Janskys, provided map
        frequency """
    wvlen = a.const.c / (freq*1e9)
    hpm.map *= 1e23*(4*np.pi/hpm.npix())*2*a.const.k/wvlen**2 # convert to Jy to make summable
    return hpm

def absorber(theta, sky, abs, nside, dB, theta_cutoff = 0.785):
    """ modify sky with absorptive baffle structure
    takes as arguments:
        theta = spherical coordinate theta
        sky = Healpix map of sky
        abs = Healpix map of flat temperature absorber (Tabs = 300 K)
          NOTE: sky and abs must have same size!
        nside = size of HEALPIX map, must be same as sky and abs
        dB = decibels of attenuation from absorber
        theta_cutoff = radians from zenith let into antenna"""
    obs = a.healpix.HealpixMap(nside=nside)
    absorber = 10**(-1.*dB/10)
    # NOTE: is this going to work???
    obs.map = np.where(theta > theta_cutoff, sky.map * absorber + abs.map * (1 - absorber), sky.map)
    return obs

#def makeGSM(path, nside, freq):
#        gsmMap = a.healpix.HealpixMap(nside=nside)
#        g = a.healpix.HealpixMap(nside=512)
#        g.map = GlobalSkyModel(freq=1000*freq, GSMlocation=path, 
#        GSMNSIDE=512).map # galactic
#        g = hpm_TtoJy(g, freq=freq)
#        gsmMap.from_hpm(g)
#        return gsmMap

#def makeGSMMap(array, nside, filename, freq, path='/home/kara/capo/kmk/gsm/gsm_raw/'):
#    """ create a Healpix map of a given size filled with a simulated global 
#    sky model
#        at a given frequency """
#    makeGSM(path=path, filename=filename, sfreq=freq, sdf=10, num=1)
#    g = a.healpix.HealpixMap(nside=512)
#    gsm = a.healpix.HealpixMap(nside=nside)
#    print "GSM map size = " + str(gsm.map.shape)
#    d = np.loadtxt(path + filename + str(1001) + '.dat')
#    g.map = d
#    g = hpm_TtoJy(g, freq)
#    gsm.from_hpm(g) # hack to lower resolution to prevent memory overload
#    # convert to topocentric coordinates
#    i, j, k = ijk = makeTop(gsm) #topocentric
#    return gsm

def makeTop(hpm):
    # returns topocentric coordinates for a given healpix map
    gxyz = gx,gy,gz, = hpm.px2crd(np.arange(hpm.npix())) # galactic
    ga2eq = a.coord.convert_m(isys='eq', osys='ga', oepoch=aa.epoch) 
    # conversion matrix so the isys and osys are reversed
    exyz = ex,ey,ez = np.dot(ga2eq, hpm.px2crd(np.arange(hpm.npix()))) # equatorial
    txyz = tx,ty,tz = np.dot(a.coord.eq2top_m(aa.sidereal_time(), aa.lat), exyz) # topocentric
    return txyz

def calcVis(aa, sky, txyz, bxyz, freq):
    # TODO: import array + GSMMap, calculate topocentric coordinates on the 
    # fly, generate PB on the fly, include time
    """ simulate sky visibilities for a given baseline and primary beam, 
        provided with a sky map at a known frequency and its coordinate system """
    # NOTE: I don't think these are in the same units and I don't know what to 
    # do about it
    tx,ty,tz = txyz # topocentric
    bx,by,bz = bxyz
    # generate proper PB
    beam = aa[0].bm_response(txyz, pol='x')**2 # topocentric
    # attenuate sky signal and visibility by primary beam
    obs_sky = beam * sky.map
    phs = np.exp(np.complex128(-2j*np.pi*freq*np.dot(bxyz, txyz)))
    vis = np.sum(np.where(tz>0, obs_sky*phs, 0))
    return vis

if __name__ == '__main__':

    import matplotlib.pylab as pl
    import optparse, sys

    import optparse, sys
    import matplotlib.pylab as pl

    o = optparse.OptionParser()
    o.add_option('--calfile', default='hyperion_deployment_aug2017')
    o.add_option('--sim_dir', default='/home/kara/capo/kmk/scripts/')
    o.add_option('--gsm_dir', default='/home/kara/capo/kmk/gsm/gsm_raw/')

    o.add_option('--fileout', default='sim_results.uv')

    opts,args = o.parse_args(sys.argv[1:])
    calfile = opts.calfile
    sim_dir = opts.sim_dir
    gsm_dir = opts.gsm_dir

    fileout = opts.fileout

    # select 150 MHz and 160 MHz for u-mode calibration test
    freqs = np.array([0.100]) 
    #freqs = np.array([0.070, 0.080, 0.090, 0.100, 0.110, 0.120])
    aa = a.cal.get_aa(calfile, freqs)

    N = 64

    # make absorber brightness at room temperature
    # make iterable over many frequencies
    I_abs = makeFlatMap(nside=N, Tsky=300.0, freq=freqs[0])
    top = I_abs.px2crd(np.arange(I_abs.npix())) #topocentric
    sph_crd = theta,phi = I_abs.px2crd(np.arange(I_abs.npix()),ncrd=2) #topocentric

    # make model sky signal (synchrotron or GSM)
    # make iterable over many frequencies
    I_sky = makeSynchMap(nside=N, freq=freqs[0])
    # modify observed sky to account for presence of absorber
    I_obs = absorber(theta=theta, sky=I_sky, abs = I_abs, nside=N, dB=15)

    # create array of baselines (in ns)
    # !!!!!!!!!!!!!!!!!!!!!!!!!
    ant0 = aa.ants[0].get_params()
    x0,y0,z0 = xyz0 = np.array([ant0['x'],ant0['y'],ant0['z']])
    ant1 = aa.ants[1].get_params()
    x1,y1,z1 = xyz1 = np.array([ant1['x'],ant1['y'],ant1['z']])
    # convert from eq to top, call get_baseline with opt 'z'
    bx,by,bz = bxyz = xyz1 - xyz0

    # number of baselines in sim, in array form
    #bl = 1
    bl = np.arange(1)

    sim_data = []
    for i in xrange(len(bl)):
        for j in xrange(len(freqs)):
            obs_vis = calcVis(aa=aa, sky=I_obs, txyz=top, bxyz=bxyz, freq=freqs[j])
            # turn into dictionary
            vis_data = [bl[i], freqs[j], obs_vis]
            sim_data.append(vis_data)

    np.array(sim_data)
    np.savez(sim_dir+'sim_output',sim_data)
    print sim_data

