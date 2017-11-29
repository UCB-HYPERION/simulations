#! /usr/bin/env python
import aipy as a
#import capo
import uvtools
import numpy as np
import os
import hypersim.absorber as ab
import matplotlib.pylab as pl
#from GlobalSkyModel import GlobalSkyModel

def makeFlatMap(nside, freq, Tsky=1.0, jansky = True):
    """ fill sky map of given size and frequency with flat temperature across whole sky,
        returns map in Janskys. """
    hpm = a.healpix.HealpixMap(nside=nside)
    hpm.map = Tsky*np.ones(shape = hpm.map.shape)
    print "flat map size = " + str(hpm.map.shape)
    if jansky == True:
        hpm = hpm_TtoJy(hpm, freq)
    return hpm

def makeSynchMap(nside, freq=0.150, jansky = True):
    """ fill sky map of given size and frequency with flat temperature across whole sky based on synchrotron spectrum,
        returns map in Janskys. """
    hpm = a.healpix.HealpixMap(nside=nside)
    Tsky = 237 * (freq/0.150)**-2.5
    hpm.map = Tsky*np.ones(shape = hpm.map.shape)
    print "synchrotron map size = " + str(hpm.map.shape)
    if jansky == True:
        hpm = hpm_TtoJy(hpm, freq)
    return hpm

def makePointSourceMap(nside, freq=0.150, Tsource=1., jansky = True):
    """ make sky map with point source with brightness temperature Tsource at zenith, for testing purposes"""
    hpm = a.healpix.HealpixMap(nside=nside)
    xyz = hpm.px2crd(np.arange(hpm.npix()))
    theta = np.arcsin(xyz[2])
    hpm.map = np.where(np.sin(theta) == 0.0, Tsource, 0)
    print "point source map size = " + str(hpm.map.shape)
    if jansky == True:
        hpm = hpm_TtoJy(hpm, freq)
    return hpm

def hpm_TtoJy(hpm, freq):
    """ converts given Healpix map from brightness temperature to Janskys, provided map
        frequency """
    wvlen = a.const.c / (freq*1e9)
    hpm.map *= 1e23*(4*np.pi/hpm.npix())*2*a.const.k/wvlen**2 # convert to Jy to make summable
    return hpm

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

def calcVis(aa, sky, nside, bl, freq, smooth, theta_cutoff, abs_file, flat = False, make_plot = True):
    # TODO: import array + GSMMap, calculate topocentric coordinates on the 
    # fly, generate PB on the fly, include time
    """ simulate sky visibilities for a given baseline and primary beam, 
        provided with a sky map at a known frequency and its coordinate system """
    tx,ty,tz= txyz = sky.px2crd(np.arange(sky.npix()))
    bxyz = aa.get_baseline(*bl, src='z')
    # generate proper PB
    abs = ab.BeamAbsorber(freqs=freq, horizon_angle=theta_cutoff)
    beam = np.abs(abs.response(txyz, smooth=smooth, data_file = abs_file, flat = False, use_abs=True))**2
    # attenuate sky signal and visibility by primary beam
    obs_sky = a.healpix.HealpixMap(nside=nside)
    obs_sky.map = beam[0] * sky.map
    if make_plot == True:
        uvtools.plot.plot_hmap_ortho(obs_sky, res=1, mx=2.5, drng=2.5)
        pl.colorbar()
        pl.show()
    phs = np.exp(np.complex128(-2j*np.pi*freq*np.dot(bxyz, txyz)))
    vis = np.sum(np.where(tz>0, obs_sky.map*phs, 0))
    return vis

def monopole_coeffs(aa, nside, bl, freq):
    # does this need to have absorber abilities?
    # I bet it does
    I_sky = makeFlatMap(nside=nside, freq=freq, Tsky=1.)
    coefficient = calcVis(aa=aa, sky=I_sky, nside=nside, bl=bl, freq=freq, smooth=0.01, theta_cutoff=np.pi/4, abs_file = False, flat=0., make_plot=False)
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
    #o.add_option('--gsm_dir', default='/home/kara/capo/kmk/gsm/gsm_raw/')

    #o.add_option('--fileout', default='sim_results.uv')

    opts,args = o.parse_args(sys.argv[1:])
    calfile = opts.calfile
    absfile = opts.absfile
    sim_dir = opts.sim_dir
    smooth = float(opts.smooth)
    #gsm_dir = opts.gsm_dir

    #fileout = opts.fileout

    # set characteristic frequencies of simulation
    freqs = np.linspace(0.050, 0.150, 50)
    wvlens = a.const.c / (freqs*1e9) # cm
        
    aa = a.cal.get_aa(calfile, freqs)

    N = 32

    # number of baselines in sim, in array form
    ij = np.arange(0,7,1)

    B = 1e9*(np.max(freqs)-np.min(freqs))/len(freqs) # channel bandwidth in Hz
    int_time = 10 # integration time in seconds
    T_sys = 300 # K, room temperature system
    T_rms = T_sys / np.sqrt(B*int_time)

    sim_data = []
    vis_data = np.zeros(len(freqs), dtype=complex)
    temp_sum = []
    temp_wgt = []
    temp = []
    pl.figure(0)
    for i in xrange(len(ij)):
        bl = (0, ij[i])
        for j in xrange(len(freqs)):
            I_sky = makeFlatMap(nside=N, freq=freqs[j], Tsky=1., jansky = True)
            #I_sky = makeSynchMap(nside=N, freq=freqs[j], jansky = True)
            cal_coeffs = monopole_coeffs(aa=aa, nside=N, bl=bl, freq=freqs[j])
            obs_vis = calcVis(aa=aa, sky=I_sky, nside=N, bl=bl, freq=freqs[j], smooth=smooth, theta_cutoff=np.pi/4, abs_file = False, flat=0., make_plot=False)
            noise_amp = 1e23 * 2*a.const.k*T_rms / (wvlens[j]**2) # convert thermal noise level to Janskys
            obs_vis += noise(shape=obs_vis.size, scale=noise_amp)
            #obs_vis = calcVis(aa=aa, sky=I_sky, nside=N, bl=bl, freq=freqs[j], smooth=smooth, theta_cutoff=np.pi/4, abs_file = absfile, flat=0, make_plot=False)
            weight = (np.abs(cal_coeffs) / noise_amp)**2
            temp_sum.append(weight * obs_vis / cal_coeffs)
            temp_wgt.append(weight)
            #temp.append(obs_vis/cal_coeffs)
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

