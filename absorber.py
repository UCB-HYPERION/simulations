'''Add support for absorbing structures to AIPY AntennaArray objects.'''
import aipy, numpy as np

class AbsorberBaffle:
    def __init__(self, dB, horizon_angle):
        self.dB = dB # XXX add freq-dependence to dB
        self.horizon_angle = horizon_angle
    def abs_response(self, xyz):
        theta = np.arccos(xyz[2])
        return np.where(theta < np.pi/2 - self.horizon_angle, 1., 10**(-self.dB/20.))

class BeamAbsorber(aipy.amp.Beam,AbsorberBaffle):
    def __init__(self, freqs, dB, horizon_angle):
        aipy.amp.Beam.__init__(self, freqs)
        AbsorberBaffle.__init__(self, dB, horizon_angle)
    def response(self, xyz, use_abs=True):
        resp = aipy.amp.Beam.response(self, xyz) 
        if use_abs: resp *= self.abs_response(xyz)
        return resp
    def auto_noise(self, Tabs, nside=64):
        h = aipy.healpix.HealpixBase(nside=nside)
        top = h.px2crd(np.arange(h.npix()))
        resp = self.response(top, use_abs=False) * (1 - self.abs_response(top))
        px_angle = 4 * np.pi / h.npix()
        lam = aipy.const.c / (self.freqs * 1e9)
        return np.sum(resp) * px_angle * 2 * aipy.const.k * Tabs / lam**2 * 1e23 # Jy 
