#!/bin/env python
import sys
from scipy import stats
from os import listdir
from os.path import join
from astropy.io.misc import fnunpickle
from modeldefs import models
import argparse
from modelmap import *

_cacheDirectory = "./.cache"

def chisqdof_to_prob(chisq, dof):
  return 1 - stats.chi2.cdf(chisq, dof)

def analyze_data(data, include_specials):
  chisqdofs = {}
  probs = {}
  lowestchisqdofs = {}
  for results in data:
    if not include_specials:
      for special in specials():
        results.pop(special)
    meta = results['meta']
    chisqbundle = {}
    probbundle = {}
    for name, result in results.iteritems():
      if name == 'meta':
        continue
      chisq = result['chisq']
      dof = result['ndof']
      chisqdof = results['chisqdof']
      prob = chisqdof_to_prob(chisq, dof)
      chisqbundle[name] = chisqdof
      probbundle[name] = prob

    lowchi = sorted([[name, chisqdof] for name, chisqdof in chisqbundle.iteritems()],
                    cmp = lambda a,b: a[1] - b[1])[0]
    chisqbundle['lowest'] = name
    probbundle['lowest'] = name
    chisqdofs[meta['SNID']] = chisqbundle
    probs[meta['SNID']] = probbundle
    if lowchi[0] in lowestchisqdofs:
      lowestchisqdofs[lowchi[0]] += 1
    else:
      lowestchisqdofs[lowchi[0]] = 1
  return chisqdofs, probs, lowestchisqdofs

def load_data():
  files = [f for f in listdir(dir) if isfile(join(dir, f)) and f.lower().endswith('.pik')]
  return [fnunpickle(join(_cacheDirectory, file)) for file in files]

def plot_types(lowestchisqdofs, show, outname):
  types = []
  for name in model_names():
    type = type_for_name(name)
    if type not in types:
      types.append(type)
  data = [0 for type in types]
  for name, val in lowestchisqdofs:
    type = type_for_name(name)
    if type == 'SN IIP':
      data[types.indexOf('SN IIP')] += val
    data[types.indexOf(type)] += val

  fig, ax = plt.subplots()
  width = 0.35
  ind = np.arange(len(types))
  ax.bar(ind, data, width, color='#8ef3af')
  ax.set_ylabel('# SN')
  ax.set_xticks(ind+width)
  ax.set_xticklabels(types)
  if show:
    plt.show()
  else:
    savefig(outname)

def main(args):
  include_special = True
  show = False
  outname = 'out.png'
  parser = argparse.ArgumentParser(description='Analyze SN Data Simulated from SNANA and fit with SNCosmo')
  parser.add_argument('-x', '--extra', action='store_const', const=False, default=True, help="exclude extra templates in analysis")
  parser.add_argument('-s', '--show', action='store_const', const=True, default=False, help="show plot instead of saving")
  parser.add_argument('-o', '--out', nargs=1, default='out.png', help='name of png to create')
  opts = parser.parse_args(args)
  include_special = opts.extra
  show = opts.show
  outname = opts.out
  data = load_data()
  chisqdofs, probs, lowestchisqdofs = analyze_data(data, include_special)
  plot_types(lowestchisqdofs, show, outname)

if __name__ == "__main__":
  main(sys.argv[1:])

