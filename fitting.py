#!/bin/env python
import sys
import argparse
import numpy as np
from astropy.table import Table
from astropy.io.misc import fnpickle, fnunpickle
import sncosmo
from modeldefs import models
from sncosmo.fitting import _nest_lc
from astropy.io import fits
import os
from os import listdir
from os.path import isfile, join, splitext

_cacheDirectory = './.cache'

def load_summary(file):
  return [int(sn) for sn in np.genfromtxt(file, usecols=(2), dtype=None)]

def filter_meta(meta, fro, to, summary):
  snids = load_summary(summary)
  metas = [met for met in meta if int(met['SNID']) in snids]
  if fro == -1:
    fro = 0
  if to == -1:
    to = len(metas)+1
  metas = metas[fro:to]
  return metas

def fit_and_save(metas, datas):
  model = sncosmo.Model(source='salt2-extended',
                        effects=[sncosmo.F99Dust()],
                        effect_names=['mw'],
                        effect_frames=['obs'])
  for meta in metas:
    results = {}
    
    dat = datas[(meta['PTROBS_MIN'] - 1):meta['PTROBS_MAX']]
    dat['FLT'][:] = np.char.strip(dat['FLT'])
    data = Table(dat)
    
    pikname = join(_cacheDirectory, meta['SNID']) + '.pik'
    if os.path.exists(pikname):
      results = fnunpickle(pikname)
      continue

    print "Fitting " + meta['SNID']

    data.rename_column('FLUXCAL', 'flux')
    data.rename_column('FLUXCALERR', 'fluxerr')
    data.rename_column('FLT', 'band')
    data.rename_column('MJD', 'time')
    data['zp'] = 27.5
    data['zpsys'] = 'ab'

    model.set(mwebv=meta['MWEBV'])
    model.set(z=meta['REDSHIFT_FINAL'])
    res, fitchisqmodel = sncosmo.fit_lc(data, model, ['t0', 'x0', 'x1', 'c'],
                                        bounds={'x1':(-5., 5.),
                                                'c':(-0.5, 0.5),
                                                't0':(np.min(data['time'])-15, np.max(data['time'])+15)},
                                        flatten=False)

    mask1 = (data['time'] < fitchisqmodel.maxtime()) & (data['time'] > fitchisqmodel.mintime()) # total epochs
    mask2 = (data['time'] < fitchisqmodel.get('t0')) & (data['time'] > fitchisqmodel.mintime())  #  before peak
    mask3 = (data['time'] > fitchisqmodel.get('t0')+10) & (data['time'] < fitchisqmodel.maxtime()) # after peak+10

    if mask1.sum()<5:
      continue
    if mask2.sum()<1:
      continue
    if mask3.sum()<1:
      continue

    data = data[mask1]
    # Get min and max data times
    dtmin = np.min(data['time'])
    dtmax = np.max(data['time'])

    results['meta'] = meta
    
    for name, m in models.iteritems():
      if name in results:
        continue
      print "\tTo model " + name
      m['model'].set(z=meta['REDSHIFT_FINAL'])      # set model redshift to host specz
      m['model'].set(mwebv=meta['MWEBV'])                 # set mwebv of model.
      t0off = m['model'].get('t0') - m['model'].mintime() # t0 offset from mintime
      t0min = dtmin - 30. + t0off
      t0max = dtmax - 30. + t0off
      m['bounds']['t0'] = (t0min, t0max)        # set t0 bounds
      m['bounds']['z'] = (meta['REDSHIFT_FINAL'], meta['REDSHIFT_FINAL'])
      res = _nest_lc(data, m['model'], m['param_names'], False,
                     bounds=m['bounds'], tied=m['tied'], nobj=50, verbose=False)
      res.chisq = -2. * res.loglmax
      res.chisqdof = res.chisq / res.ndof
      results[name] = res
    
    fnpickle(results, pikname)

def open_fits(dir):
  fitsfiles = [f for f in listdir(dir) if isfile(join(dir, f)) and f.lower().endswith('.fits')]
  headfile = photfile = None
  heads = [f for f in fitsfiles if splitext(f)[0].lower().endswith("head")]
  phots = [f for f in fitsfiles if splitext(f)[0].lower().endswith("phot")]
  if len(heads) == 0:
    raise Exception('Could not find metadata \'head\' file in {0}'.format(dir))
  if len(phots) == 0:
    raise Exception('Could not find data \'head\' file in {0}'.format(dir))
  hdrfits = fits.open(join(dir, heads[0]))
  datfits = fits.open(join(dir, phots[0]))
  filterdir = [f for f in listdir(dir) if not isfile(join(dir, f)) and f.lower().endswith('.filters')][0]
  banddir = join(dir, filterdir)
  for name in ['g', 'r', 'i', 'z']:
    filename = os.path.join(banddir, name + '.dat')
    band = sncosmo.read_bandpass(filename, name=name)
    sncosmo.registry.register(band)
  return hdrfits[1].data, datfits[1].data

def main(args):
  
  summary = "PassAllCuts500.summary"
  fro = -1
  to = -1
  directory = "/home/kuhlmann/snana/root_v201204/SIM/DES_5years_CC_v1033f/"

  parser = argparse.ArgumentParser(description='Analyze SN Data Simulated from SNANA')
  parser.add_argument('-f', '--fro', type=int, nargs=1, help="from where")
  parser.add_argument('-t', '--to', type=int, nargs=1, help="to where")
  parser.add_argument('-s', '--summary', nargs=1, help="location of .summary file to consult")
  parser.add_argument('-d', '--dir', nargs=1, help="directory of the FITS files to analyze")
  opts = parser.parse_args(args)
  if opts.fro is not None:
    fro = opts.fro[0]
  if opts.to is not None:
    to = opts.to[0]
  if opts.dir is not None:
    directory = opts.dir[0]
  if opts.summary is not None:
    summary = opts.summary[0]
  meta, data = open_fits(directory)
  meta = filter_meta(meta, fro, to, summary)
  if not os.path.exists(_cacheDirectory):
    os.makedirs(_cacheDirectory)
  fit_and_save(meta, data)

if __name__ == "__main__":
  main(sys.argv[1:])

