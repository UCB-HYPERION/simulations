'''Add support for absorbing structures to AIPY AntennaArray objects.'''
import aipy, numpy as np
from scipy import interpolate
import csv, sys

class BeamDipole(aipy.phs.Beam):
    def __init__(self, freqs):
        aipy.phs.Beam.__init__(self, np.array([freqs]))
    def bm_response(self, xyz):
        x,y,z = np.array(xyz)
        ones = np.ones((self.afreqs.size, x.size))
        theta = np.arcsin(xyz[2])
        resp = np.where(np.sin(theta) == 0.0, 0, np.cos(np.pi / 2. * np.cos(theta)) / np.sin(theta))
        return resp*ones

class Absorber:
    def __init__(self, horizon_angle):
        self.horizon_angle = horizon_angle
    def freq_atten(self, data_file):
        with open(data_file, 'rb') as f:
            reader = csv.reader(f, delimiter=',')
            freqs = []
            attens = []
            for row in reader:
                freqs.append(float(row[0])/1.e9) # put frequencies into GHz to match everything else
                attens.append(float(row[1])) # values are negative
            func = interpolate.interp1d(freqs, attens)
            return func
    def abs_response(self, xyz, freq, smooth, data_file, flat):
        theta = np.arccos(xyz[2])
        s = 0.5 + 0.5*np.tanh((theta-(np.pi/2 - self.horizon_angle))/smooth)
        if bool(data_file) == True:
            atten = self.freq_atten(data_file)
            dB = atten(freq)
        else:
            dB = flat
        print "dB = %f" % dB
        return s * 10**(dB/20.) + (1-s)*1.

class BeamAbsorber(BeamDipole,Absorber):
    def __init__(self, freqs, horizon_angle):
        BeamDipole.__init__(self, freqs)
        Absorber.__init__(self, horizon_angle)
    def response(self, xyz, smooth, data_file, flat, use_abs=True):
        resp = self.bm_response(xyz)
        if use_abs: resp *= self.abs_response(xyz, self.freqs, smooth, data_file, flat)
        return resp
    def auto_noise(self, Tabs, nside=64):
        h = aipy.healpix.HealpixBase(nside=nside)
        top = h.px2crd(np.arange(h.npix()))
        resp = self.response(top, use_abs=False) * (1 - self.abs_response(top))
        px_angle = 4 * np.pi / h.npix()
        lam = aipy.const.c / (self.freqs * 1e9)
        return np.sum(resp) * px_angle * 2 * aipy.const.k * Tabs / lam**2 * 1e23 # Jy 
