
import os
import optparse
import numpy as np
from scipy.interpolate import interp1d
from astropy import cosmology
from astropy.io.misc import fnpickle
import sncosmo
from sncosmo.photdata import standardize_data
from sncosmo.utils import Interp1D
from collections import OrderedDict as odict
from os.path import join, isfile

models = odict()
amplitude0 = {}

# Define the dust model used everywhere.
dust = sncosmo.F99Dust(3.1)

# Create a fast function for distance modulus (used in tying absolute
# magnitude to amplitude parameters)
cosmo = cosmology.FlatLambdaCDM(H0=70., Om0=0.3)
zgrid = np.linspace(0.001, 2., 400)
dmgrid = cosmo.distmod(zgrid)
dm = interp1d(zgrid, dmgrid)
tomabs = Interp1D(0., 1., np.array([-16., -20.]))

def add_salt2():
  model = sncosmo.Model(source='salt2-extended', effects=[dust],
                        effect_names=['mw'],
                        effect_frames=['obs'])
  param_names = ['z', 't0', 'x0', 'x1', 'c']
  amplitude0['salt2-extended'] = model.get('x0')
  ppfs = {'x0': lambda amp, n='salt2-extended': amplitude0[n] * 10.**(-0.4 * (tomabssalt(amp) + dm(models[n]['model'].get('z'))))}
  bounds = {'x1': (-3., 3.),
            'c': (-0.3, 0.3),
            'z': (0.001, 1.2),
            'mabs':(-17.5, -20.)}
  tomabssalt = Interp1D(0., 1., np.array([-16., -20.]))
  model.source.set_peakmag(0., 'bessellb', 'ab')
  models['salt2-extended'] = {'type': 'SN Ia',
                              'mprior': 0.5,
                              'model': model,
                              'param_names': param_names,
                              'bounds': bounds,
                              'ppfs': ppfs}

def add_builtin():
  for name, sntype in [('s11-2004hx', 'SN IIL/P'),
                       ('s11-2005lc', 'SN IIP'),
                       ('s11-2005hl', 'SN Ib'),
                       ('s11-2005hm', 'SN Ib'),
                       ('s11-2005gi', 'SN IIP'),
                       ('s11-2006fo', 'SN Ic'),
                       ('s11-2006jo', 'SN Ib'),
                       ('s11-2006jl', 'SN IIP')]:
    model = sncosmo.Model(source=name,
                          effects=[dust, dust],
                          effect_names=['host', 'mw'],
                          effect_frames=['rest', 'obs'])

    amplitude0[name] = model.get('amplitude')
    ppfs = {'amplitude': lambda amp, n=name: amplitude0[n] * 10**(-0.4 * (tomabs(amp) + dm(models[n]['model'].get('z'))))}
    bounds = {'hostebv': (-0.2, 0.5)}

    model.source.set_peakmag(0., 'bessellb', 'ab')

    models[name] = {'type': sntype,
                    'mprior': 0.5/8.,
                    'model': model,
                    'param_names': ['t0', 'amplitude', 'hostebv'], # print out these params from pickle
                    'bounds': bounds,
                    'ppfs': ppfs}


def add_others():
  templatedir = '/fusion/gpfs/home/kuhlmann/snana/root_v201204/snsed/non1a/'

  if isfile('config.py'):
    from config import config
    templatedir = os.path.expanduser(config['templateDirectory']) if 'templateDirectory' in config else templatedir

  customTemplates = []

  for name, type, file in [
      ('CSP-2006ep', 'SN Ib', join(templatedir, 'CSP-2006ep.SED')),
      ('SDSS-017548', 'SN Ic', join(templatedir, 'SDSS-017548.SED')),
      ('SDSS-000018', 'SN IIP', join(templatedir, 'SDSS-000018.SED')),
      ('SDSS-004012', 'SN Ic', join(templatedir, 'SDSS-004012.SED')),
      ('CSP-2004fe', 'SN Ic', join(templatedir, 'CSP-2004fe.SED')),
      ('SDSS-018457', 'SN IIP', join(templatedir, 'SDSS-018457.SED')),
      ('SDSS-014492', 'SN Ib', join(templatedir, 'SDSS-014492.SED'))
      ]:
    phase, wave, flux = sncosmo.io.read_griddata_ascii(file)
    min = np.min(wave)
    max = np.max(wave)
    if min > 1000:
      padding = np.arange(1000.0, min, 10.0)
      fluxpad = np.zeros(len(padding))
      wave = np.append(padding, wave)
      newflux = []
      for f in flux:
        newflux.append(np.concatenate((fluxpad, f)))
      flux = np.array(newflux)
    if max < 14900:
      padding = np.arange(max+10.0, 14910.0, 10.0)
      fluxpad = np.zeros(len(padding))
      wave = np.append(wave, padding)
      newflux = []
      for f in flux:
        newflux.append(np.concatenate((f, fluxpad)))
      flux = np.array(newflux)

    timeSeries = sncosmo.TimeSeriesSource(phase, wave, flux)
    model = sncosmo.Model(source=timeSeries,
                          effects=[dust, dust],
                          effect_names=['host', 'mw'],
                          effect_frames=['rest', 'obs'])

    amplitude0[name] = model.get('amplitude')
    ppfs = {'amplitude': lambda amp, n=name: amplitude0[n] * 10**(-0.4 * (tomabs(amp) + dm(models[n]['model'].get('z'))))}
    bounds = {'hostebv': (-0.2, 0.5)}

    model.source.set_peakmag(0., 'bessellb', 'ab')

    models[name] = {'type': type,
                    'mprior': 0.5/8.,
                    'model': model,
                    'param_names': ['t0', 'amplitude', 'hostebv'], # print out these params from pickle
                    'bounds': bounds,
                    'ppfs': ppfs}
    customTemplates.append(name)
 
add_salt2()
add_builtin()
add_others()

