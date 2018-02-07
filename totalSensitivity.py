#! /bin/env python

import matplotlib.pylab as pl
import numpy as np
import aipy as a

# compute weighted average of all observed monopole measurements from multiple 
# baselines

def read_spec(filename):
    uv = a.miriad.UV(filename)
    bl = []
    spec = []
    nchan = uv['nchan']
    start_freq = uv['sfreq']
    df = uv['sdf']
    stop_freq = start_freq - (df * (nchan-1))
    freqs = np.linspace(start_freq, stop_freq, nchan)
    for i in xrange(uv['nants']):
        preamble,uv_spec,flags = uv.read(raw=True)
        bl.append(preamble[2])
        spec.append(uv_spec)
    bl = np.array(bl)
    spec = np.array(spec)
    spec.shape = (-1,nchan)
    return bl,freqs,spec

def weighted_avg(vals, weights):
    num = np.sum(weights*vals,axis=0)
    denom = np.sum(weights**2,axis=0)
    return num / denom

if __name__ == '__main__':

    bl,freqs,data = read_spec('sim_results.uv')
    c_bl,c_freqs,coeffs = read_spec('coefficients.uv')

    # sanity check
    print np.array_equal(bl,c_bl)
    print np.array_equal(freqs,c_freqs)

    total_vis = weighted_avg(data,coeffs)
    print total_vis

    pl.plot(freqs,total_vis)
    pl.show()


    




