import numpy as np, aipy as a
import unittest
import globalSignalVis
import hypersim as h
import matplotlib.pylab as pl

class calcVisTest(unittest.TestCase):
    def setUp(self):
        self.N = 512
        self.freqs = np.array([0.100])
        self.smooth = 0.00001
        self.aa = a.cal.get_aa('hyperion_deployment_aug2017', self.freqs)
    def testCalcVis(self):
        bl = (0,0)
        dB = 0.
        theta_cutoff = np.pi / 6
        wvlens = a.const.c / (self.freqs * 1e9)
        sky = globalSignalVis.makeFlatMap(nside=self.N, freq=self.freqs[0])
        # pre-divide by beam to negate it in calcVis
        bm = h.absorber.BeamDipole(freqs=self.freqs[0])
        beam = np.abs(bm.beam_response(xyz = sky.px2crd(np.arange(sky.npix()))))**2
        obs = a.healpix.HealpixMap(nside=self.N)
        obs.map = sky.map / beam[0]
        v = globalSignalVis.calcVis(aa=self.aa, sky=obs, nside=self.N, bl=bl, freq=self.freqs[0], smooth = self.smooth, abs_file = False, flat = dB, theta_cutoff=theta_cutoff, make_plot=False)
        vis_predict = 1e23 * (2 * a.const.k / wvlens[0]**2) * 2 * np.pi
        # verify that visibility matches predicted value
        self.assertAlmostEqual(vis_predict, np.abs(v), -1)
        # verify that visibility has no phase
        self.assertEqual(np.imag(v), 0)

        # verify point source characteristics
        sky = globalSignalVis.makePointSourceMap(nside=self.N, freq=self.freqs[0], Tsource=1.)
        ij = np.arange(0,7,1)
        vis = 0
        for i in xrange(len(ij)):
            bl = (0, ij[i])
            v = globalSignalVis.calcVis(aa=self.aa, sky=sky, nside=self.N, bl=bl, freq=self.freqs[0], smooth = self.smooth, abs_file = False, flat = dB, theta_cutoff=theta_cutoff, make_plot=False)
            if i == 0:
                vis = v
            self.assertEqual(np.imag(v), 0)
            self.assertEqual(v, vis)


if __name__ == '__main__':
    unittest.main()
        


