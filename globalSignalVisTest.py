import numpy as np, aipy as a
import unittest
import globalSignalVis

class calcVisTest(unittest.TestCase):
    def setUp(self):
        self.N = 512
        self.freqs = np.array([0.100])
        self.aa = a.cal.get_aa('hyperion_deployment_aug2017', self.freqs)
    def testCalcVis(self):
        bl = (0,0)
        dB = 0.
        theta_cutoff = np.pi / 6
        wvlens = a.const.c / (self.freqs * 1e9)
        sky = globalSignalVis.makeFlatMap(nside=self.N, freq=self.freqs[0])
        v = globalSignalVis.calcVis(aa=self.aa, sky=sky, nside=self.N, bl=bl, freq=self.freqs[0], dB=dB, theta_cutoff=theta_cutoff)
        vis_predict = 1e23 * (2 * a.const.k / wvlens[0]**2) * 2 * np.pi
        self.assertAlmostEqual(vis_predict, v, -1)

if __name__ == '__main__':
    unittest.main()
        


