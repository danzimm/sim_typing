
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

# Define the dust model used everywhere.
dust = sncosmo.F99Dust(3.1)

def add_salt2():
  model = sncosmo.Model(source='salt2-extended', effects=[dust],
                        effect_names=['mw'],
                        effect_frames=['obs'])
  param_names = ['t0', 'x1', 'c', 'x0', 'z']
  bounds = {'x1': (-3., 3.),
            'c': (-0.3, 0.3)}
  model.source.set_peakmag(0., 'bessellb', 'ab')
  models['salt2-extended'] = {'type': 'SN Ia',
                              'mprior': 0.5,
                              'model': model,
                              'param_names': param_names,
                              'bounds': bounds}

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

    #print "time series for {} min phase: {}".format(name, model.source.minphase())
    bounds = {'hostebv': (-0.2, 0.5)}

    model.source.set_peakmag(0., 'bessellb', 'ab')

    models[name] = {'type': sntype,
                    'mprior': 0.5/8.,
                    'model': model,
                    'param_names': ['t0', 'hostebv', 'amplitude', 'z'], # print out these params from pickle
                    'bounds': bounds}


customTemplates = []
def add_others():
  global customTemplates
  templatedir = '/fusion/gpfs/home/kuhlmann/snana/root_v201204/snsed/non1a/'

  if isfile('config.py'):
    from config import config
    templatedir = os.path.expanduser(config['templateDirectory']) if 'templateDirectory' in config else templatedir

  for name, type in [
      ('CSP-2006ep', 'SN Ib'),
      ('CSP-2004gq', 'SN Ic'),
      ('SDSS-017548', 'SN Ic'),
      ('SDSS-004012', 'SN Ic'),
      ('CSP-2004fe', 'SN Ic'),
      ('SDSS-018457', 'SN IIP'),
      ('Nugent+Scolnic_IIL', 'SN IIL'),
      ('CSP-2004gv', 'SN Ib'),
      ('CSP-2007Y', 'SN Ib'),
      ('SDSS-000018', 'SN IIP'),
      ('SDSS-000020', 'SN Ib'),
      ('SDSS-002744', 'SN Ib'),
      ('SDSS-003818', 'SN IIP'),
      ('SDSS-012842', 'SN IIn'),
      ('SDSS-013195', 'SN Ic'),
      ('SDSS-013376', 'SN IIP'),
      ('SDSS-013449', 'SN IIn'),
      ('SDSS-014450', 'SN IIP'),
      ('SDSS-014475', 'SN Ic'),
      ('SDSS-014492', 'SN Ib'),
      ('SDSS-014599', 'SN IIP'),
      ('SDSS-015031', 'SN IIP'),
      ('SDSS-015320', 'SN IIP'),
      ('SDSS-015339', 'SN IIP'),
      ('SDSS-015475', 'SN Ic'),
      ('SDSS-017564', 'SN IIP'),
      ('SDSS-017862', 'SN IIP'),
      ('SDSS-018109', 'SN IIP'),
      ('SDSS-018297', 'SN IIP'),
      ('SDSS-018408', 'SN IIP'),
      ('SDSS-018441', 'SN IIP'),
      ('SDSS-018590', 'SN IIP'),
      ('SDSS-018596', 'SN IIP'),
      ('SDSS-018700', 'SN IIP'),
      ('SDSS-018713', 'SN IIP'),
      ('SDSS-018734', 'SN IIP'),
      ('SDSS-018793', 'SN IIP'),
      ('SDSS-018834', 'SN IIP'),
      ('SDSS-018892', 'SN IIP'),
      ('SDSS-019323', 'SN Ib'),
      ('SDSS-020038', 'SN IIP'),
      ('SNLS-04D1la', 'SN Ibc'),
      ('SNLS-04D4jv', 'SN Ic'),
      ]:
    if name in models:
      print "Warning: Adding " + str(name) + " to models again" # my bad I made a syntax error
    file = join(templatedir, name + '.SED')
    phase, wave, flux = sncosmo.io.read_griddata_ascii(file)
    # investigate why minphase is negative for some of these???
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
    #print "time series for {} min phase: {}".format(name, timeSeries.minphase())
    model = sncosmo.Model(source=timeSeries,
                          effects=[dust, dust],
                          effect_names=['host', 'mw'],
                          effect_frames=['rest', 'obs'])

    bounds = {'hostebv': (-0.2, 0.5)}

    model.source.set_peakmag(0., 'bessellb', 'ab')

    models[name] = {'type': type,
                    'mprior': 0.5/8.,
                    'model': model,
                    'param_names': ['t0', 'hostebv', 'amplitude', 'z'], # print out these params from pickle
                    'bounds': bounds}
    customTemplates.append(name)
 
add_salt2()
add_builtin()
add_others()

