#! /usr/bin/env python
import aipy as a
import uvtools
import numpy as np
from scipy import interpolate
import os
import hypersim.absorber as ab
import matplotlib.pylab as pl
import sys
sys.path.append('/home/kara/capo/kmk/scripts')
import ares
from GlobalSkyModel import GlobalSkyModel

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
    dTb = sim.history['dTb']/1.e3 #convert temperatures from mK to K
    gs_func = interpolate.interp1d(nu, dTb)
    gs = gs_func(freqs)
    return gs 

def makeFlatMap(nside, freq, file_path, Tsky=1.0):
    """ fill sky map of given size and frequency with flat temperature across whole sky,
        returns map in Janskys. """
    filename = file_path+"flat_map_%dK_%dMHz.fits" % (int(Tsky), int(1e3*freq))
    hpm = a.healpix.HealpixMap(nside=nside)
    if os.path.isfile(filename) == True:
        hpm.from_fits(filename)
    else:
        hpm.map = Tsky*np.ones(shape = hpm.map.shape)
        hpm.to_fits(filename)
    return hpm

def makeSynchMap(nside, file_path, freq=0.150):
    """ fill sky map of given size and frequency with flat temperature across whole sky based on synchrotron spectrum,
        returns map in Janskys. """
    filename = file_path+"synch_map_%dK_%dMHz.fits" % (int(Tsky), int(1e3*freq))
    hpm = a.healpix.HealpixMap(nside=nside)
    Tsky = 237 * (freq/0.150)**-2.5
    if os.path.isfile(filename) == True:
        hpm.from_fits(filename)
    else:
        hpm.map = Tsky*np.ones(shape = hpm.map.shape)
        hpm.to_fits(filename)
    return hpm

def makePointSourceMap(nside, xyz, file_path, freq=0.150, Tsource=1.):
    """ make sky map with point source with brightness temperature Tsource at zenith, for testing purposes"""
    filename = file_path+"point_map_%dK_%dMHz.fits" % (int(Tsky), int(1e3*freq))
    hpm = a.healpix.HealpixMap(nside=nside)
    if os.path.isfile(filename) == True:
        hpm.from_fits(filename)
    else:
        x,y,z = xyz
        hpm[x,y,z] = Tsource
        hpm.to_fits(filename)
    return hpm

def makeGSM(nside, freq, file_path, path='/home/kara/capo/kmk/gsm/gsm_raw/'):
    """ make sky map with Global Sky Model, in Kelvin"""
    gsmMap = a.healpix.HealpixMap(nside=nside)
    filename = file_path+"gsm_map_%dMHz.fits" % int(1000*freq)
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
    beam = np.abs(abs.response(txyz, smooth=smooth, data_file = abs_file, flat = flat, use_abs=True))**2
    # attenuate sky signal and visibility by primary beam
    obs_sky = beam[0] * sky.map
    #obs_sky = sky.map
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

def monopole_coeffs(aa, file_path, nside, bl, freq, txyz, flat, theta_cutoff, abs_file, smooth):
    I_sky = makeFlatMap(nside=nside, freq=freq, file_path=file_path, Tsky=1.)
    coefficient = calcVis(aa=aa, sky=I_sky, txyz = txyz, nside=nside, bl=bl, freq=freq, smooth=smooth, theta_cutoff=theta_cutoff, abs_file = abs_file, flat=flat, make_plot=False)
    return coefficient

def noise(shape, scale):
    n_real = np.random.normal(size=shape, scale=scale/np.sqrt(2))
    n_imag = np.random.normal(size=shape, scale=scale/np.sqrt(2))
    return n_real + 1j*n_imag

def init_uv(filename, aa, freqs, inttime, sys_opts):
    uv = a.miriad.UV(filename, status='new')
    uv._wrhd('obstype','mixed-auto-cross')
    uv._wrhd('history','globalSignalVis: created file.\nglobalSignalVis: ' + sys_opts + '\n')
    uv.add_var('telescop','a'); uv['telescop'] = 'HYPERION'
    uv.add_var('operator','a'); uv['operator'] = 'AIPY'
    uv.add_var('version' ,'a'); uv['version'] = '0.0.1'
    uv.add_var('epoch'   ,'r'); uv['epoch'] = 2000.
    uv.add_var('source'  ,'a'); uv['source'] = 'zenith'
    uv.add_var('latitud' ,'d'); uv['latitud'] = aa.lat
    uv.add_var('dec'     ,'d'); uv['dec'] = aa.lat
    uv.add_var('obsdec'  ,'d'); uv['obsdec'] = aa.lat
    uv.add_var('longitu' ,'d'); uv['longitu'] = aa.long
    uv.add_var('npol'    ,'i'); uv['npol'] = 1
    uv.add_var('nspect'  ,'i'); uv['nspect'] = 1
    uv.add_var('nants'   ,'i'); uv['nants'] = len(aa)
    uv.add_var('antpos'  ,'d')
    antpos = np.array([ant.pos for ant in aa], dtype=np.double)
    uv['antpos'] = antpos.transpose().flatten()
    uv.add_var('sfreq'   ,'d'); uv['sfreq'] = freqs[0]
    uv.add_var('freq'    ,'d'); uv['freq'] = freqs[0]
    uv.add_var('restfreq','d'); uv['restfreq'] = freqs[0]
    uv.add_var('sdf'     ,'d'); uv['sdf'] = freqs[1]-freqs[0]
    uv.add_var('nchan'   ,'i'); uv['nchan'] = len(freqs)
    uv.add_var('nschan'  ,'i'); uv['nschan'] = len(freqs)
    uv.add_var('inttime' ,'r'); uv['inttime'] = float(inttime)
    # These variables just set to dummy values
    uv.add_var('vsource' ,'r'); uv['vsource'] = 0.
    uv.add_var('ischan'  ,'i'); uv['ischan'] = 1
    uv.add_var('tscale'  ,'r'); uv['tscale'] = 0.
    uv.add_var('veldop'  ,'r'); uv['veldop'] = 0.
    # These variables will get updated every spectrum
    uv.add_var('coord'   ,'d')
    uv.add_var('time'    ,'d')
    uv.add_var('lst'     ,'d')
    uv.add_var('ra'      ,'d')
    uv.add_var('obsra'   ,'d')
    uv.add_var('baseline','r')
    uv.add_var('pol'     ,'i')
    return uv

if __name__ == '__main__':

    import optparse, sys

    o = optparse.OptionParser()
    o.add_option('--calfile', default='hyperion_deployment_aug2017')
    o.add_option('--absfile', default='DIP_S11_LIDON_FER_DB.csv')
    o.add_option('--sim_dir', default='/home/kara/capo/kmk/scripts/')
    o.add_option('--smooth', default=0.01)
    o.add_option('--gsm_dir', default='/home/kara/capo/kmk/gsm/gsm_raw/')
    o.add_option('--plot_uv', default=False)
    o.add_option('--fileout', default='sim_results.uv')

    opts,args = o.parse_args(sys.argv[1:])
    calfile = opts.calfile
    absfile = opts.absfile
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
    T_rms = 1.
    
    uv = init_uv(filename=fileout, aa=aa, freqs=freqs, inttime=int_time, sys_opts=sys_opts)

    sim_data = []
    spec = np.zeros(len(freqs), dtype=complex)
    temp_sum = []
    temp_wgt = []
    temp = []
    pl.figure(0)
    
    flat = 0.
    I_sky_eq = a.healpix.HealpixMap(nside=N)
    for i in xrange(len(ij)):
        bl = (0, ij[i])
        # set starting time of simulated observation
        aa.set_jultime(jds[0])
        # set necessary parameters for each input spectrum for data file
        # NOTE: THIS WILL NEED MASSAGING WHEN WE REINCORPORATE TIME
        uv['lst'] = aa.sidereal_time()
        uv['ra'] = aa.sidereal_time()
        uv['obsra'] = aa.sidereal_time()
        uv['pol'] = a.miriad.str2pol['xx'] # not dealing with polarization in this sim, so set them all to same value
        crd = aa.get_baseline(0, ij[i])
        preamble = (crd, jd_init, bl)
        for j in xrange(len(freqs)):
            vis_sum_time = 0
            flatMap = makeFlatMap(nside=N, freq=freqs[j], file_path=file_path, Tsky=1.)
            #point_source_loc = np.array([np.cos(aa.lat)*np.cos(aa.sidereal_time()), np.cos(aa.lat)*np.sin(aa.sidereal_time()), np.sin(aa.lat)])
            #pointMap = makePointSourceMap(nside=N, xyz=point_source_loc, file_path=file_path, freq = freqs[j], Tsource=100.)
            #I_sky_eq.map = flatMap.map + pointMap.map
            I_sky_eq.map = flatMap.map
            #exyz = ex,ey,ez = I_sky_eq.px2crd(np.arange(I_sky_eq.npix())) # topocentric
            txyz = tx,ty,tz = I_sky_eq.px2crd(np.arange(I_sky_eq.npix())) # topocentric
            # NOTE: NEED TO PROBABLY PULL TIME TO BE OUTER LOOP TO MAKE UV 
            # FILES
            #for k in xrange(len(jds)):
            #    # convert maps from equatorial to topocentric coordinates
            #    aa.set_jultime(jds[k])
            #    eq2top = aa.eq2top_m
            #    txyz = tx,ty,tz = np.dot(eq2top, exyz)
            cal_coeffs = monopole_coeffs(aa=aa, file_path=file_path, txyz = txyz, nside=N, bl=bl, freq=freqs[j], smooth=smooth, theta_cutoff=np.pi/4, abs_file = False, flat=flat)
            obs_vis = calcVis(aa=aa, sky=I_sky_eq, txyz = txyz, nside=N, bl=bl, freq=freqs[j], smooth=smooth, theta_cutoff=np.pi/4, abs_file = False, flat=flat, make_plot=False)
            noise_amp = T_rms # convert thermal noise level to Janskys
                #obs_vis += noise(shape=obs_vis.size, scale=noise_amp)
            #obs_vis = calcVis(aa=aa, sky=I_sky_eq, nside=N, bl=bl, freq=freqs[j], txyz = txyz, smooth=smooth, theta_cutoff=np.pi/4, abs_file = absfile, flat=0, make_plot=False)
            avg_vis = obs_vis
            weight = (np.abs(cal_coeffs) / noise_amp)**2
            # no noise for now
            temp_sum.append(weight * avg_vis / cal_coeffs)
            temp_wgt.append(weight)
            spec[j] = avg_vis
        uv.write(preamble, data=spec, flags = np.zeros(len(freqs)))
        bx,by,bz = bxyz = aa.get_baseline(*bl, src='z')
        uv_crds = np.sqrt(bx**2 + by**2 + bz**2) / wvlens
        if plot_uv == True:
            pl.plot(uv_crds, np.abs(spec), label="Baseline %d" % i) 
        else:
            pl.plot(freqs, np.abs(spec), label=uv_sep[i]) 

    #pl.suptitle(r"\huge{Monopole Sky without Absorber, Flat Brightness}") 
    pl.suptitle(r"\huge{Monopole Sky with Ferrite Absorber, Flat Brightness}") 
    #pl.title(r"\LARGE{Absorber Smoothing Factor %f}" % smooth) 
    if uv == True:
        pl.xlabel(r"\Large{uv-plane Baseline Separation}")
    else:
        pl.xlabel(r"\Large{Frequency (GHz)}")
    pl.ylabel(r"\Large{Visibility Amplitude}")
    pl.legend()

    temp_sum = np.array(temp_sum)
    temp_wgt = np.array(temp_wgt)

    temp_sum.shape = (-1,50)
    temp_wgt.shape = (-1,50)

    # solve for the observed temperature using inverse variance weighting
    temp = np.sum(temp_sum,axis=0)/np.sum(temp_wgt,axis=0)
    temp_rms = 1./np.sqrt(np.sum(temp_wgt,axis=0)) #rms noise of inverse variance weighted values 
    temp_rms = np.zeros(temp_rms.shape)
    print temp
    #print np.mean(temp), np.mean(temp_rms)


    #pl.figure(1)
    #pl.plot(freqs, temp, label="Measured Temperature")
    #pl.plot(freqs, temp_rms)
    #pl.errorbar(freqs, temp, yerr=temp_rms)
    #pl.title(r"RMS Noise")
    #pl.xlabel(r"Frequency (GHz)")
    #pl.ylabel(r"RMS Noise Amplitude")
    #pl.legend()
    pl.show()

