import aipy as a
import hypersim

def get_aa(freqs):
# Define the location of your instrument
     lat, lon = '37.2314', '-118.2941'
     #lat, lon = '90.', '-118.2941'
     # Create a model of the primary beam.  Beam is a minimal model that has 
     # unity gain in all directions.  
     beam = hypersim.absorber.BeamDipole(freqs)   
     # Make a list of antennas with requisite nanosecond locations, primary 
     # beams, and any other calibration parameters you wish to provide.
     ants = [
         a.fit.Antenna(  0,   0,  0, beam),
         a.fit.Antenna(  0,   5,  0, beam), # 100 MHz lambda/2
         a.fit.Antenna(  0,   10,  0, beam), # 100 MHz 1 lambda
         a.fit.Antenna(  0,   15,  0, beam), # 100 MHz 3 lambda/2
         a.fit.Antenna(  0,   20,  0, beam), # 100 MHz 2 lambda
         a.fit.Antenna(  0,   30,  0, beam), # 100 MHz 3 lambda
         a.fit.Antenna(  0,   40,  0, beam), # 100 MHz 4 lambda
     ]
     
     # Create an AntennaArray at the specified location with the listed 
     # antennas
     aa = a.fit.AntennaArray((lat,lon), ants)
     return aa

src_prms = {
'cen':{ 'jys':10**3.282102, 'index':  0.235166 , },
'cyg':{ 'jys':10**3.566410, 'index':  -0.266315 , },
'hyd':{ 'jys':10**2.448816, 'index':  -0.866462 , },
'pic':{ 'jys':10**2.714456, 'index':  -0.436361 , },
'vir':{ 'jys':10**2.200725, 'index':  0.202425 , },
'Sun': {'a1': 0.00644, 'index': 1.471, 'a2': 0.00586, 'jys': 55701.96, 'th': -0.000512},
'for': {'a1': 0.00851, 'a2': 0.00413, 'jys': 907.09, 'th': 0.230},
}
                         
def get_catalog(srcs=None, cutoff=None, catalogs=['helms','misc']):
# Pass off the request for sources to the AIPY source catalog. If desired, you 
# can substitute your own sources or source calibrations.
    if srcs is None:
        cat = a.src.get_catalog(srcs=srcs, cutoff=cutoff, catalogs=catalogs)
    else:
        cat = a.src.get_catalog(srcs=[s for s in srcs], cutoff=cutoff, catalogs=catalogs)
    cat.set_params(src_prms)
    return cat
