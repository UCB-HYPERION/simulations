import aipy as a
import re
import numpy as np

def init_uv(filename, aa, freqs, inttime, sys_opts):
    uv = a.miriad.UV(filename, status='new')
    uv._wrhd('obstype', 'mixed-auto-cross')
    uv._wrhd('history', 'globalSignalVis: created file.\nglobalSignalVis: ' + sys_opts + '\n')
    antpos = np.array([ant.pos for ant in aa], dtype=np.double)
    parameters = '''
        telescop (a) = 'HYPERION'
        operator (a) = 'AIPY'
        version (a) = '0.0.1'
        epoch (r) = 2000.
        source (a) = 'zenith'
        latitud (d) = aa.lat
        dec (d) = aa.lat
        obsdec (d) = aa.lat
        longitu (d) = aa.long
        npol (i) = 1
        nspect (i) = 1
        nants (i) = len(aa)
        antpos (d) = antpos.transpose().flatten()
        sfreq (d) = freqs[0]
        freq (d) = freqs[0]
        restfreq (d) = freqs[0]
        sdf (d) = freqs[1]-freqs[0]
        nchan (i) = len(freqs)
        nschan (i) = len(freqs)
        inttime (r) = float(inttime)

        # These variables just set to dummy values
        vsource (r) = 0.
        ischan (i) = 1
        tscale (r) = 0.
        veldop (r) = 0.

        # These variables will get updated every spectrum
        coord (d)
        time (d)
        lst (d)
        ra (d)
        obsra (d)
        baseline (r)
        pol (i)
    '''
    regex = r'(?P<name>\w+)\s+\((?P<type>\w+)\)(\s*=\s*(?P<value>.*))?'
    for i, line in enumerate(parameters.split('\n')):
        declaration = line.split('#')[0].strip()  # discard comments
        if not declaration:
            continue  # ignore blank lines
        match = re.match(regex, declaration)
        if match:
            name = match.group('name')
            type = match.group('type')
            value = match.group('value')
            uv.add_var(name, type)
            if value:
                try:
                    evaluated = eval(value, None, locals())
                except Exception as e:
                    print 'init_uv %d: %s: %s' % (i, value, e)
                uv[name] = evaluated
        else:
            print 'init_uv %d: garbled parameter definition: %s' % (i, declaration)

    return uv
