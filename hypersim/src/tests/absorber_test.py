import unittest
import aipy as a, numpy as np, hypersim as h

class TestBeamDipole(unittest.TestCase):
    def setUp(self):
        #self.fq = np.arange(0,1,0.1) #GHz
        self.fq = 0.100 #GHz
        self.bm = h.absorber.BeamDipole(self.fq)
    def test_beam_esponse(self):
        xyz = np.array([0,0,1])
        self.assertEqual(self.bm.beam_response(xyz),1)
        xyz = np.array([0,1,0])
        self.assertEqual(self.bm.beam_response(xyz),0)
class TestAbsorber(unittest.TestCase):
    def setUp(self):
        self.fq = 0.100 #GHz
        self.smooth = 1e-6
        self.dB = -15.
        self.ang = np.pi / 6
        self.abs = h.absorber.Absorber(self.ang)
    def test_attributes(self):
        """Test that Absorber has correct attributes."""
        self.assertTrue(np.all(self.abs.horizon_angle == self.ang))
    def test_abs_response(self):
        xyz = np.array([0,0,1])
        self.assertEqual(self.abs.abs_response(xyz, self.fq, self.smooth, data_file=False, flat=self.dB),1)
        xyz = np.array([0,1,0])
        self.assertEqual(self.abs.abs_response(xyz, self.fq, self.smooth, data_file=False, flat=self.dB),10**(self.dB/20.))
        ang1 = self.ang + 0.01
        ang2 = self.ang - 0.01
        xyz = np.array([0,np.cos(ang1),np.sin(ang1)])
        self.assertEqual(self.abs.abs_response(xyz, self.fq, self.smooth, data_file=False, flat=self.dB),1)
        xyz = np.array([0,np.cos(ang2),np.sin(ang2)])
        self.assertEqual(self.abs.abs_response(xyz, self.fq, self.smooth, data_file=False, flat=self.dB),10**(self.dB/20.))

class TestBeamAbsorber(unittest.TestCase):
    def setUp(self):
        self.fq = np.arange(0,1,0.1)
        self.smooth = 1e-6
        self.ang = np.pi / 4
        self.bm = h.absorber.BeamAbsorber(self.fq, self.ang)
    def test_attributes(self):
        """Test that BeamAbsorber has correct attributes."""
        self.assertTrue(np.all(self.bm.freqs == self.fq))
        self.assertTrue(np.all(self.bm.horizon_angle == self.ang))
    def test_response(self):
        """Test that BeamAbsorber gives same results for no absorber and 0 dB 
        attenuation absorber."""
        dB = 0.
        h = a.healpix.HealpixBase(nside=64)
        top = h.px2crd(np.arange(h.npix()))
        no_abs = self.bm.response(top, self.smooth, data_file=False, flat=dB, use_abs=False)
        abs_0dB = self.bm.response(top, self.smooth, data_file=False, flat=dB, use_abs=True)
        np.testing.assert_equal(no_abs, abs_0dB)
        dB = -15.
        abs_15dB = self.bm.response(top, self.smooth, data_file=False, flat=dB, use_abs=True)
        self.assertFalse(np.all(no_abs == abs_15dB))
        ang1 = self.ang + 0.01
        ang2 = self.ang - 0.01
        xyz = np.array([0,np.cos(ang1),np.sin(ang1)])
        np.testing.assert_equal(self.bm.response(xyz, self.smooth, data_file=False, flat=dB, use_abs=True),self.bm.beam_response(xyz))
        xyz = np.array([0,np.cos(ang2),np.sin(ang2)])
        np.testing.assert_equal(self.bm.response(xyz, self.smooth, data_file=False, flat=dB, use_abs=True),self.bm.beam_response(xyz)*10**(-15/20.))


if __name__ == '__main__':
    unittest.main()

