
import numpy as np
from os.path import isfile, join, splitext
from os import listdir
from astropy.io import fits
import sncosmo

def load_summary(file):
  return [int(sn) for sn in np.genfromtxt(file, usecols=(2), dtype=None)]

def open_sim_fits(dir):
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
    filename = join(banddir, name + '.dat')
    band = sncosmo.read_bandpass(filename, name=name)
    sncosmo.registry.register(band)
  meta = hdrfits[1].data
  data = datfits[1].data
  hdrfits.close()
  datfits.close()
  return meta, data

