#! /usr/bin/env python
import aipy as a
import uvtools
import numpy as np
from scipy import interpolate
import os
import hypersim.absorber as ab
import matplotlib.pylab as pl
import sys
import ares
import maps
from init_uv import init_uv

pl.rc('text', usetex=True)
pl.rc('font', family='serif')
pl.rc('xtick', labelsize=14)
pl.rc('ytick', labelsize=14)

def makeGlobalSignal(freqs):
    """ return the theoretical differential brightness temperatures across chosen frequency range"""
    sim = ares.simulations.Global21cm()
    sim.run()
    print "simulation successfully run"
    nu = sim.history['nu']/1.e3 # convert frequencies from MHz to GHz
    dTb = sim.history['dTb']/1.e3 # convert temperatures from mK to K
    gs_func = interpolate.interp1d(nu, dTb)
    gs = gs_func(freqs)
    return gs 

def calcVis(aa, sky, txyz, nside, bl, freq, smooth, theta_cutoff, absorber, make_plot = True):
    """ simulate sky visibilities for a given baseline and primary beam, 
        provided with a sky map at a known frequency and its coordinate system """
    tx,ty,tz = txyz
    bxyz = aa.get_baseline(*bl, src='z')
    # generate proper PB
    abs = ab.BeamAbsorber(freqs=freq, horizon_angle=theta_cutoff)
    beam = np.abs(abs.response(txyz, smooth=smooth, absorber = absorber, use_abs=True))**2
    # attenuate sky signal and visibility by primary beam
    obs_sky = beam[0] * sky.map
    if make_plot == True:
        pl_map = a.healpix.HealpixMap(nside=nside, interp=True)
        x,y,z = xyz = pl_map.px2crd(np.arange(pl_map.npix()))
        top2eq = np.linalg.inv(aa.eq2top_m)
        ex,ey,ez = np.dot(top2eq, xyz)
        old_map = sky.map
        sky.map = obs_sky
        sky.set_interpol(True)
        pl_map.map = sky[ex,ey,ez] 
        sky.map = old_map
        sky.set_interpol(False)
        uvtools.plot.plot_hmap_ortho(pl_map, res=0.25, mx=2.5, drng=2.5)
        pl.colorbar()
        pl.show()
    phs = np.exp(np.complex128(-2j*np.pi*freq*np.dot(bxyz, txyz)))
    vis = (4*np.pi/sky.npix())*np.sum(np.where(tz>0, obs_sky*phs, 0))
    return vis

def monopole_coeffs(aa, map, nside, bl, freq, txyz, theta_cutoff, absorber, smooth):
    coefficient = calcVis(aa=aa, sky=map, txyz = txyz, nside=nside, bl=bl, freq=freq, smooth=smooth, theta_cutoff=theta_cutoff, absorber = absorber, make_plot=False)
    return coefficient

def noise(shape, scale):
    n_real = np.random.normal(size=shape, scale=scale/np.sqrt(2))
    n_imag = np.random.normal(size=shape, scale=scale/np.sqrt(2))
    return n_real + 1j*n_imag

if __name__ == '__main__':

    import optparse, sys

    o = optparse.OptionParser()
    o.add_option('--calfile', default='hyperion_deployment_aug2017')
    o.add_option('--abs', default='DIP_S11_LIDON_FER_DB.csv')
    o.add_option('--sim_dir', default='/home/kara/capo/kmk/scripts/')
    o.add_option('--smooth', default=0.01)
    o.add_option('--gsm_dir', default='/home/kara/capo/kmk/gsm/gsm_raw/')
    o.add_option('--plot_uv', default=False)
    o.add_option('--fileout', default='sim_results.uv')

    opts,args = o.parse_args(sys.argv[1:])
    calfile = opts.calfile
    abs = opts.abs
    sim_dir = opts.sim_dir
    smooth = float(opts.smooth)
    gsm_dir = opts.gsm_dir
    plot_uv = opts.plot_uv
    fileout = opts.fileout

    sys_opts = ' '.join(sys.argv)

    file_path = "/home/kara/src/python/berkeley/simulations/maps/"

    # set characteristic frequencies of simulation
    freqs = np.linspace(0.050, 0.150, 50)
    wvlens = a.const.c / (freqs*1e9) # cm

    # calculate global signal temperature values over our frequencies
    T_gs = makeGlobalSignal(freqs)
        
    # set time of observation
    aa = a.cal.get_aa(calfile, freqs)
    jd_init = 2458080
    #jds = np.linspace(jd_init, jd_init + 1, 10)
    jds = np.array([jd_init])

    N = 64

    # number of baselines in sim, in array form
    ij = np.arange(0,7,1)
    uv_sep = np.array(['Auto-Correlation',r'$\lambda/2$',r'$\lambda$',r'$3\lambda/2$',r'$2\lambda$',r'$3\lambda$',r'$4\lambda$'])

    B = 1e9*(np.max(freqs)-np.min(freqs))/len(freqs) # channel bandwidth in Hz
    int_time = 10 # integration time in seconds
    T_sys = 300 # K, room temperature system
    T_rms = T_sys / np.sqrt(B*int_time)
    
    uv = init_uv(filename=fileout, aa=aa, freqs=freqs, inttime=int_time, sys_opts=sys_opts)
    cal_uv = init_uv(filename='coefficients.uv', aa=aa, freqs=freqs, inttime=int_time, sys_opts=sys_opts)

    print uv.vars()
    print uv['nchan']

    sim_data = []
    spec = np.zeros(len(freqs), dtype=complex)
    cal_spec = np.zeros(len(freqs), dtype=complex)
    temp_sum = np.zeros((len(ij),len(freqs)), dtype=complex)
    temp_wgt = np.zeros((len(ij),len(freqs)), dtype=complex)
    pl.figure(0)
    
    # make sky maps for all frequencies, stick them into one object
    coeffSky = maps.FlatMap(nside=N, freqs=freqs, filepath=file_path, Tsky=1.)
    globalSignalSky = maps.GlobalSignalMap(freqs=freqs, nside=N, filepath=file_path, Tgs=T_gs)
    #point_source_loc = np.array([np.cos(aa.lat)*np.cos(aa.sidereal_time()), np.cos(aa.lat)*np.sin(aa.sidereal_time()), np.sin(aa.lat)])
    #pointMap = maps.PointMap(nside=N, xyz=point_source_loc, file_path=file_path, freq = freqs[j], Tsource=100.)
    #sky = globalSignalSky
    sky = coeffSky

    for i in xrange(len(ij)):
        bl = (0, ij[i])
        for k in xrange(len(jds)):
            # set starting time of simulated observation
            aa.set_jultime(jds[k])
            # set necessary parameters for each input spectrum for data file
            uv['lst'] = cal_uv['lst'] = aa.sidereal_time()
            uv['ra'] = cal_uv['ra'] = aa.sidereal_time()
            uv['obsra'] = cal_uv['obsra'] = aa.sidereal_time()
            uv['pol'] = cal_uv['pol'] = a.miriad.str2pol['xx'] # not dealing with polarization in this sim, so set them all to same value
            crd = aa.get_baseline(0, ij[i])
            preamble = (crd, jd_init, bl)
            f = 0
            # generate visibility spectra for observed sky and calibration
            for freq,map in sky:
                # convert maps from equatorial to topocentric coordinates
                eq2top = aa.eq2top_m
                exyz = ex,ey,ez = map.px2crd(np.arange(map.npix())) # equatorial
                txyz = tx,ty,tz = np.dot(eq2top, exyz)
                # calculate monopole calibration coefficients
                cal_coeffs = monopole_coeffs(aa=aa, map=coeffSky.maps[f], txyz = txyz, nside=N, bl=bl, freq=freq, smooth=smooth, theta_cutoff=np.pi/4, absorber=abs)
                # calculate monopole visibilities
                obs_vis = calcVis(aa=aa, sky=map, txyz = txyz, nside=N, bl=bl, freq=freq, smooth=smooth, theta_cutoff=np.pi/4, absorber=abs, make_plot=False)
                noise_amp = T_rms # convert thermal noise level to Janskys
                # add noise to simulated visibility
                #obs_vis += noise(shape=obs_vis.size, scale=noise_amp)
                weight = (np.abs(cal_coeffs) / noise_amp)**2
                spec[f] = obs_vis
                cal_spec[f] = cal_coeffs
                temp_sum[i][f] = weight * obs_vis / cal_coeffs
                temp_wgt[i][f] = weight
                f += 1
            uv.write(preamble, data=spec, flags = np.zeros(len(freqs)))
            cal_uv.write(preamble, data=cal_spec, flags = np.zeros(len(freqs)))
        bx,by,bz = bxyz = aa.get_baseline(*bl, src='z')
        uv_crds = np.sqrt(bx**2 + by**2 + bz**2) / wvlens
        if plot_uv == True:
            pl.plot(uv_crds, np.abs(spec), label="Baseline %d" % i) 
        else:
            pl.plot(freqs, np.abs(spec), label=uv_sep[i]) 
    
    # delete pointers to uv files or they won't close properly and things get 
    # fucked up real quick
    del(uv)
    del(cal_uv)

    #pl.suptitle(r"\huge{Monopole Sky without Absorber, Flat Brightness}") 
    pl.suptitle(r"\huge{Monopole Sky with Ferrite Absorber, Flat Brightness}") 
    #pl.title(r"\LARGE{Absorber Smoothing Factor %f}" % smooth) 
    if plot_uv == True:
        pl.xlabel(r"\Large{uv-plane Baseline Separation}")
    else:
        pl.xlabel(r"\Large{Frequency (GHz)}")
    pl.ylabel(r"\Large{Visibility Amplitude}")
    pl.legend()

    # solve for the observed temperature using inverse variance weighting
    temp_sum *= 1. / len(jds)
    temp_wgt *= 1. / len(jds)
    temp = np.sum(temp_sum,axis=0)/np.sum(temp_wgt,axis=0)
    temp_rms = 1./np.sqrt(np.sum(temp_wgt,axis=0)) #rms noise of inverse variance weighted values 
    print np.mean(temp), np.mean(temp_rms)

    # show that data output into UV file properly, is readable
    pl.figure(1)
    uv_data = a.miriad.UV(fileout)
    for i in xrange(len(ij)):
        preamble,uv_spec,flags = uv_data.read(raw=True)
        bl = preamble[2]
        pl.plot(freqs, np.abs(uv_spec),label="Baseline %d" % bl[1])
    pl.show()

