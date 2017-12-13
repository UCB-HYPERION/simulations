#! /usr/bin/env python
import aipy as a
import uvtools
import numpy as np
import os
import hypersim.absorber as ab
import matplotlib.pylab as pl
import sys
sys.path.append('/home/kara/capo/kmk/scripts')
from GlobalSkyModel import GlobalSkyModel

def makeFlatMap(nside, freq, Tsky=1.0):
    """ fill sky map of given size and frequency with flat temperature across whole sky,
        returns map in Janskys. """
    hpm = a.healpix.HealpixMap(nside=nside)
    hpm.map = Tsky*np.ones(shape = hpm.map.shape)
    print "flat map size = " + str(hpm.map.shape)
    return hpm

def makeSynchMap(nside, freq=0.150):
    """ fill sky map of given size and frequency with flat temperature across whole sky based on synchrotron spectrum,
        returns map in Janskys. """
    hpm = a.healpix.HealpixMap(nside=nside)
    Tsky = 237 * (freq/0.150)**-2.5
    hpm.map = Tsky*np.ones(shape = hpm.map.shape)
    print "synchrotron map size = " + str(hpm.map.shape)
    return hpm

def makePointSourceMap(nside, xyz, freq=0.150, Tsource=1.):
    """ make sky map with point source with brightness temperature Tsource at zenith, for testing purposes"""
    hpm = a.healpix.HealpixMap(nside=nside)
    x,y,z = xyz
    hpm[x,y,z] = Tsource
    print "point source map size = " + str(hpm.map.shape)
    return hpm

def makeGSM(nside, freq, path='/home/kara/capo/kmk/gsm/gsm_raw/'):
    """ make sky map with Global Sky Model, in Kelvin"""
    # I'm p sure GlobalSkyModel returns maps in brightness temp from the source 
    # code, but not 100% sure??
    gsmMap = a.healpix.HealpixMap(nside=nside)
    filename = "gsm_%dMHz.fits" % int(1000*freq)
    if os.path.isfile(filename)==True:
        gsmMap.from_fits(filename)
    else:
        g = a.healpix.HealpixMap(nside=512)
        g.map = GlobalSkyModel(freq=1000*freq, GSMlocation=path, GSMNSIDE=512).m
        gsmMap.from_hpm(g)
        gsmMap.to_fits(filename)
    return gsmMap

def calcVis(aa, sky, txyz, nside, bl, freq, smooth, theta_cutoff, abs_file, flat = False, make_plot = True):
    """ simulate sky visibilities for a given baseline and primary beam, 
        provided with a sky map at a known frequency and its coordinate system """
    tx,ty,tz = txyz
    bxyz = aa.get_baseline(*bl, src='z')
    # generate proper PB
    abs = ab.BeamAbsorber(freqs=freq, horizon_angle=theta_cutoff)
    beam = np.abs(abs.response(txyz, smooth=smooth, data_file = abs_file, flat = False, use_abs=True))**2
    # attenuate sky signal and visibility by primary beam
    obs_sky = beam[0] * sky.map
    print obs_sky.max()
    if make_plot == True:
        pl_map = a.healpix.HealpixMap(nside=nside, interp=True)
        x,y,z = xyz = pl_map.px2crd(np.arange(pl_map.npix()))
        top2eq = np.linalg.inv(aa.eq2top_m)
        ex,ey,ez = np.dot(top2eq, xyz)
        # no beam in plot because I can't figure out how to make it work lmao
        pl_map.map = sky[ex,ey,ez] 
        # disappearing pixel at zenith even with interpolation on
        uvtools.plot.plot_hmap_ortho(pl_map, res=1, mx=2.5, drng=2.5)
        pl.colorbar()
        pl.show()
    phs = np.exp(np.complex128(-2j*np.pi*freq*np.dot(bxyz, txyz)))
    vis = (4*np.pi/sky.npix())*np.sum(np.where(tz>0, obs_sky*phs, 0))
    return vis

def monopole_coeffs(aa, nside, bl, freq):
    # does this need to have absorber abilities?
    # I bet it does
    I_sky = makeFlatMap(nside=nside, freq=freq, Tsky=1.)
    txyz = I_sky.px2crd(np.arange(I_sky.npix()))
    coefficient = calcVis(aa=aa, sky=I_sky, txyz = txyz, nside=nside, bl=bl, freq=freq, smooth=0.01, theta_cutoff=np.pi/4, abs_file = False, flat=0., make_plot=False)
    return coefficient

def noise(shape, scale):
    n_real = np.random.normal(size=shape, scale=scale/np.sqrt(2))
    n_imag = np.random.normal(size=shape, scale=scale/np.sqrt(2))
    return n_real + 1j*n_imag

if __name__ == '__main__':

    import optparse, sys

    o = optparse.OptionParser()
    o.add_option('--calfile', default='hyperion_deployment_aug2017')
    o.add_option('--absfile', default='DIP_S11_LIDON_FER_DB.csv')
    o.add_option('--sim_dir', default='/home/kara/capo/kmk/scripts/')
    o.add_option('--smooth', default=0.1)
    o.add_option('--gsm_dir', default='/home/kara/capo/kmk/gsm/gsm_raw/')

    #o.add_option('--fileout', default='sim_results.uv')

    opts,args = o.parse_args(sys.argv[1:])
    calfile = opts.calfile
    absfile = opts.absfile
    sim_dir = opts.sim_dir
    smooth = float(opts.smooth)
    gsm_dir = opts.gsm_dir

    #fileout = opts.fileout

    # set characteristic frequencies of simulation
    freqs = np.linspace(0.050, 0.150, 50)
    wvlens = a.const.c / (freqs*1e9) # cm
        
    aa = a.cal.get_aa(calfile, freqs)
    jd_init = 2458080
    jds = np.linspace(jd_init, jd_init + 0.15, 4)

    N = 32

    # number of baselines in sim, in array form
    ij = np.arange(0,7,1)

    B = 1e9*(np.max(freqs)-np.min(freqs))/len(freqs) # channel bandwidth in Hz
    int_time = 10 # integration time in seconds
    T_sys = 300 # K, room temperature system
    #T_rms = T_sys / np.sqrt(B*int_time)
    T_rms = 0.

    sim_data = []
    vis_data = np.zeros(len(freqs), dtype=complex)
    temp_sum = []
    temp_wgt = []
    temp = []
    pl.figure(0)
    
    I_sky_eq = a.healpix.HealpixMap(nside=N)
    for i in xrange(len(ij)):
        bl = (0, ij[i])
        for j in xrange(len(freqs)):
            # make sky with point source at zenith in equatorial coordinates
            aa.set_jultime(jds[0])
            flatMap = makeFlatMap(nside=N, freq=freqs[j], Tsky=1.)
            point_source_loc = np.array([np.cos(aa.lat)*np.cos(aa.sidereal_time()), np.cos(aa.lat)*np.sin(aa.sidereal_time()), np.sin(aa.lat)])
            pointMap = makePointSourceMap(nside=N, xyz=point_source_loc, freq = freqs[j], Tsource=100.)
            I_sky_eq.map = flatMap.map + pointMap.map
            exyz = ex,ey,ez = I_sky_eq.px2crd(np.arange(I_sky_eq.npix())) # topocentric
            for k in xrange(len(jds)):
                # convert maps from equatorial to topocentric coordinates
                aa.set_jultime(jds[k])
                eq2top = aa.eq2top_m
                txyz = tx,ty,tz = np.dot(eq2top, exyz)
                cal_coeffs = monopole_coeffs(aa=aa, nside=N, bl=bl, freq=freqs[j])
                obs_vis = calcVis(aa=aa, sky=I_sky_eq, txyz = txyz, nside=N, bl=bl, freq=freqs[j], smooth=smooth, theta_cutoff=np.pi/4, abs_file = False, flat=0., make_plot=True)
                noise_amp = 1e23 * 2*a.const.k*T_rms / (wvlens[j]**2) # convert thermal noise level to Janskys
                #obs_vis += noise(shape=obs_vis.size, scale=noise_amp)
                #obs_vis = calcVis(aa=aa, sky=I_sky, nside=N, bl=bl, freq=freqs[j], smooth=smooth, theta_cutoff=np.pi/4, abs_file = absfile, flat=0, make_plot=False)
                weight = (np.abs(cal_coeffs) / noise_amp)**2
                temp_sum.append(weight * obs_vis / cal_coeffs)
                temp_wgt.append(weight)
            vis_data[j] = obs_vis
        bx,by,bz = bxyz = aa.get_baseline(*bl, src='z')
        uv = np.sqrt(bx**2 + by**2 + bz**2) / wvlens
        pl.plot(freqs, np.abs(vis_data), label="Baseline %d" % i) 

    #pl.title("Absorber Smoothing Factor = %f" % smooth)
    pl.title("Flat Sky with Absorber, Flat Brightness") 
    #pl.title("Synchrotron Sky with No Absorber") 
    pl.xlabel("uv-plane Baseline Separation")
    pl.xlabel("Frequency (GHz)")
    pl.ylabel("Visibility Amplitude")
    pl.legend()

    temp_sum = np.array(temp_sum)
    temp_wgt = np.array(temp_wgt)
    #temp = np.array(temp)

    temp_sum.shape = (-1,50)
    temp_wgt.shape = (-1,50)
    #temp.shape = (-1, 50)

    # solve for the observed temperature using inverse variance weighting
    temp = np.sum(temp_sum,axis=0)/np.sum(temp_wgt,axis=0)
    temp_rms = 1./np.sqrt(np.sum(temp_wgt,axis=0)) #rms noise of inverse variance weighted values 
    #temp_rms_gen = np.sqrt(np.sum(temp_wgt**2,axis=0))/np.sqrt(np.sum(temp_wgt,axis=0)) #rms noise of inverse variance weighted values 
    print temp
    #print np.mean(temp), np.mean(temp_rms)


    pl.figure(1)
    #pl.plot(freqs, temp, label="Measured Temperature")
    #pl.plot(freqs, temp_rms)
    pl.errorbar(freqs, temp, yerr=temp_rms)
    pl.title("RMS Noise")
    pl.xlabel("Frequency (GHz)")
    pl.ylabel("RMS Noise Amplitude")
    #pl.legend()
    pl.show()

