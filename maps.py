import numpy as np
import aipy as a
import os, sys
sys.path.append('/home/kara/capo/kmk/scripts')
from GlobalSkyModel import GlobalSkyModel

class HealpixMapMF(object):
    """Collection of utilities for mapping data on a sphere. Adds multi-frequency functionality to the data map infrastructure of aipy.healpix.HealpixMap"""
    def __init__(self, freqs, nside, filepath):
        self.freqs = list(freqs)
        self.nside = nside
        self.filepath = filepath
        maps = []
        for f in self.freqs:
            hpm = a.healpix.HealpixMap(nside=self.nside)
            maps.append(hpm)
        self.maps = maps
    def get_map(self, freq):
        index = self.freqs.index(freq)
        return self.maps[index]
    def __iter__(self):
        for f,m in zip(self.freqs, self.maps):
            yield f,m
    def read_maps(self, dirname):
        for f,m in self:
            filename = dirname+"%dMHz.fits" % int(1e3*f)
            m.from_fits(filename)
    def save_maps(self, dirname):
        """ make directory for all maps FITS files, save all maps to 
        disk """
        try:
            os.makedirs(dirname)
            for f,m in self:
                filename = dirname+"%dMHz.fits" % int(1e3*f)
                m.to_fits(filename)
        except (IOError, OSError) as e:
            if e.errno != e.EEXIST:
                raise SystemExit('{}: {}'.format(e.filename, e.strerror))

class FlatMap(HealpixMapMF):
    """ make multi-frequency sky map filled with spatially flat emission 
        with brightness temperature Tsky, in Kelvin """
    def __init__(self, freqs, nside, filepath, Tsky=1.):
        super(FlatMap, self).__init__(freqs, nside, filepath)
        dirname = "%sflat_map_%dK_%d-%dMHz_%d/" % (self.filepath, Tsky, 1e3*min(self.freqs), 1e3*max(self.freqs), len(self.freqs))
        if os.path.isdir(dirname):
            self.read_maps(dirname)
        else:
            for f,m in self:
                m.map = Tsky*np.ones(shape=m.map.shape)
            self.save_maps(dirname)

class GlobalSignalMap(HealpixMapMF):
    """ make multi-frequency sky map filled with spatially flat emission 
        with brightness temperature of predicted hydrogen global signal, in Kelvin """
    def __init__(self, freqs, nside, filepath, Tgs):
        super(GlobalSignalMap, self).__init__(freqs, nside, filepath)
        dirname = "%sglobal_signal_map_%d-%dMHz_%d/" % (self.filepath, 1e3*min(self.freqs), 1e3*max(self.freqs), len(self.freqs))
        if os.path.isdir(dirname):
            self.read_maps(dirname)
        else:
            for f,m in self:
                index = self.freqs.index(f)
                m.map = Tgs[index]*np.ones(shape=m.map.shape)
            self.save_maps(dirname)

class SynchrotronMap(HealpixMapMF):
    """ make multi-frequency sky map filled with spatially flat synchrotron 
        emission, in Kelvin """
    def __init__(self, freqs, nside, filepath):
        super(SynchrotronMap, self).__init__(freqs, nside, filepath)
        dirname = "%ssynch_map_%d-%dMHz_%d/" % (self.filepath, 1e3*min(self.freqs), 1e3*max(self.freqs), len(self.freqs))
        if os.path.isdir(dirname):
            self.read_maps(dirname)
        else:
            for f,m in self:
                Tsky = 237 * (f/0.150)**-2.5
                m.map = Tsky*np.ones(shape=m.map.shape)
            self.save_maps(dirname)

class PointMap(HealpixMapMF):
    """ make multi-frequency sky map with point source with brightness 
        temperature Tsource at given location, for testing purposes """
    def __init__(self, freqs, nside, filepath, xyz, Tsource=1.):
        super(PointMap, self).__init__(freqs, nside, filepath)
        # include coordinates in dirname?
        dirname = "%spoint_map_%dK_%d-%dMHz_%d/" % (self.filepath, Tsky, 1e3*min(self.freqs), 1e3*max(self.freqs), len(self.freqs))
        if os.path.isdir(dirname):
            self.read_maps(dirname)
        else:
            for f,m in self:
                x,y,z = xyz
                m[x,y,z] = Tsource
            self.save_maps(dirname)

class GSMMap(HealpixMapMF):
    """ make multi-frequency sky map with Global Sky Model, in Kelvin """
    def __init__(self, freqs, nside, filepath, path='/home/kara/capo/kmk/gsm/gsm_raw/', Tsky=1.):
        super(GSMMap, self).__init__(freqs, nside, filepath)
        dirname = "%sgsm_map_%d-%dMHz_%d/" % (self.filepath, 1e3*self.freqs.min(), 1e3*self.freqs.max(), len(self.freqs))
        if os.path.isdir(dirname):
            self.read_maps(dirname)
        else:
            g = a.healpix.HealPix(nside=512)
            for f,m in self:
                g.map = GlobalSkyModel(freq=1e3*f, GSMlocation=path, GSMNSIDE=512).map
                m.fromhpm(g)
            self.save_maps(dirname)

