import unittest, absorber
import aipy as a, numpy as np

class TestAbsorberBaffle(unittest.TestCase):
    def setUp(self):
        self.atten = 15
        self.ang = np.pi / 6
        self.abs = absorber.AbsorberBaffle(self.atten, self.ang)
    def test_attributes(self):
        """Test that AbsorberBaffle has correct attributes."""
        self.assertTrue(np.all(self.abs.dB == self.atten))
        self.assertTrue(np.all(self.abs.horizon_angle == self.ang))
    def test_abs_response(self):
        xyz = np.array([0,0,1])
        self.assertEqual(self.abs.abs_response(xyz),1)
        xyz = np.array([0,1,0])
        self.assertEqual(self.abs.abs_response(xyz),10**(-self.abs.dB/20.))
        ang1 = self.ang + 0.01
        ang2 = self.ang - 0.01
        xyz = np.array([0,np.cos(ang1),np.sin(ang1)])
        self.assertEqual(self.abs.abs_response(xyz),1)
        xyz = np.array([0,np.cos(ang2),np.sin(ang2)])
        self.assertEqual(self.abs.abs_response(xyz),10**(-self.abs.dB/20.))

class TestBeamAbsorber(unittest.TestCase):
    def setUp(self):
        self.fq = np.arange(0,1,0.1)
        self.atten = 0 # dB
        self.ang = np.pi / 4
        self.bm = absorber.BeamAbsorber(self.fq, self.atten, self.ang)
    def test_attributes(self):
        """Test that BeamAbsorber has correct attributes."""
        self.assertTrue(np.all(self.bm.freqs == self.fq))
        self.assertTrue(np.all(self.bm.dB == self.atten))
        self.assertTrue(np.all(self.bm.horizon_angle == self.ang))
    def test_response(self):
        """Test that BeamAbsorber gives same results for no absorber and 0 dB 
        attenuation absorber."""
        h = a.healpix.HealpixBase(nside=64)
        top = h.px2crd(np.arange(h.npix()))
        no_abs = self.bm.response(top, use_abs=False)
        abs_0dB = self.bm.response(top, use_abs=True)
        np.testing.assert_equal(no_abs, abs_0dB)
        abs_15dB = absorber.BeamAbsorber(self.fq, 15, self.ang)
        self.assertFalse(np.all(no_abs == abs_15dB))
        ang1 = self.ang + 0.01
        ang2 = self.ang - 0.01
        xyz = np.array([0,np.cos(ang1),np.sin(ang1)])
        np.testing.assert_equal(abs_15dB.response(xyz, use_abs=True),1)
        xyz = np.array([0,np.cos(ang2),np.sin(ang2)])
        np.testing.assert_equal(abs_15dB.response(xyz, use_abs=True),10**(-15/20.))


if __name__ == '__main__':
    unittest.main()

