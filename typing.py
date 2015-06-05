#!/usr/bin/env python
import sys
from scipy import stats
from os import listdir
from os.path import join, isfile
from astropy.io.misc import fnunpickle
from modeldefs import models
import argparse
from modelmap import *
import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.io import fits
import sncosmo
from astropy.table import Table
from snutils import open_sim_fits
from triangle import corner

_cacheDirectory = "./.cache"
_isMCMC = False

def chisqdof_to_prob(chisq, dof):
  return 1 - stats.chi2.cdf(chisq, dof)

def sntypes_equiv(typea, typeb):
  typea = typea.replace('SN ', '')
  typeb = typeb.replace('SN ', '')
  if (typea == 'IIL/P' and (typeb == 'IIP' or typeb == 'IIL')) or (typeb == 'IIL/P' and (typea == 'IIP' or typea == 'IIL')):
    return True
  return typea == typeb

def analyze_data(data, include_specials):
  chisqdofs = {}
  probs = {}
  lowestchisqdofs = {}
  lowestcorrect = {}
  lowestincorrect = {}
  for results in data:
    if not include_specials:
      for special in specials():
        results.pop(special, None)
    meta = results['meta']
    chisqbundle = {}
    probbundle = {}
    if not _isMCMC:
      for name, result in results.iteritems():
        if name == 'meta':
          continue
        chisq = result['chisq']
        dof = result['ndof']
        chisqdof = result['chisqdof']
        prob = chisqdof_to_prob(chisq, dof)
        chisqbundle[name] = chisqdof
        probbundle[name] = prob

      orderedchis = sorted([[name, chisqdof] for name, chisqdof in chisqbundle.iteritems()],
                           cmp = lambda a,b: -1 if a[1] - b[1] < 0 else (0 if a[1] == b[1] else 1))
      orderedprobs = sorted([[name, prob] for name, prob in probbundle.iteritems()],
                            cmp = lambda a,b: 1 if a[1] - b[1] < 0 else (0 if a[1] == b[1] else -1))

      lowchi = orderedchis[0]
      realtype = meta['SIM_TYPE_NAME']
      correctlytyped = sntypes_equiv(realtype, type_for_name(lowchi[0]))

      chisqbundle['correct'] = probbundle['correct'] = correctlytyped
      chisqbundle['ordered'] = orderedchis
      probbundle['ordered'] = orderedprobs

      if lowchi[0] in lowestchisqdofs:
        lowestchisqdofs[lowchi[0]] += 1
      else:
        lowestchisqdofs[lowchi[0]] = 1

      adder = lowestcorrect if correctlytyped else lowestincorrect
      if lowchi[0] in adder:
        adder[lowchi[0]] += 1
      else:
        adder[lowchi[0]] = 1

    chisqbundle['meta'] = probbundle['meta'] = meta
    chisqbundle['results'] = probbundle['results'] = results
    chisqdofs[meta['SNID']] = chisqbundle
    probs[meta['SNID']] = probbundle

  return chisqdofs, probs, lowestchisqdofs, lowestcorrect, lowestincorrect

def load_data():
  files = [f for f in listdir(_cacheDirectory) if isfile(join(_cacheDirectory, f)) and f.lower().endswith('.pik')]
  return [fnunpickle(join(_cacheDirectory, file)) for file in files]

def plot_types(lowestchisqdofs, show, outname):
  types = []
  for name in model_names():
    type = type_for_name(name)
    if type not in types:
      types.append(type)
  data = [0 for type in types]
  for name, val in lowestchisqdofs.iteritems():
    type = type_for_name(name)
    #if type == 'SN IIP':
    #  data[types.index('SN IIL/P')] += val
    data[types.index(type)] += val

  fig, ax = plt.subplots()
  width = 0.35
  ind = np.arange(len(types))
  ax.bar(ind, data, width, color='#418654')
  ax.set_ylabel('# SN')
  ax.set_xticks(ind+(width/2))
  ax.set_xticklabels(types)
  if show:
    plt.show()
  else:
    plt.savefig(outname, dpi=300)

def histo_probs(probs, show, outname):
  lowestprobs = [prob['ordered'][0][1] for snid, prob in probs.iteritems()]
  histodata(lowestprobs, 'Probability', '# SN', show, outname)

def histo_chisqdiff(chisqs, show, outname):
  toptwodiff = [chisq['ordered'][1][1] - chisq['ordered'][0][1] for snid, chisq in chisqs.iteritems()]
  histodata(toptwodiff, 'Difference in top two chisq/dof', '# SN', show, outname)

def histo_probdiff(probs, show, outname):
  toptwodiff = [prob['ordered'][0][1] - prob['ordered'][1][1] for snid, prob in probs.iteritems()]
  histodata(toptwodiff, 'Difference in top two probabilities', '# SN', show, outname)

def histodata(data, xl, yl, show, outname):
  fig, ax = plt.subplots()
  ax.hist(data, 30, color='#418654')
  ax.set_ylabel(yl)
  ax.set_xlabel(xl)
  if show:
    plt.show()
  else:
    plt.savefig(outname, dpi=300)

def filter_probabilities(probs, corrects, incorrects, greater, cutoff, condition=None):
  cutter = (lambda x: x >= cutoff) if greater else (lambda y: y <= cutoff)
  subcutter = condition if condition is not None else (lambda bundle: True)
  retval = {}
  for snid, prob in probs.iteritems():
    if not cutter(prob['ordered'][0][1] - prob['ordered'][1][1]):
      continue
    if not subcutter(prob):
      continue
    if prob['correct'] and corrects:
      retval[snid] = prob
    elif (not prob['correct']) and incorrects:
      retval[snid] = prob
  return retval

def print_prob_info(probs):
  for snid, prob in probs.iteritems():
    print "\tLooking at {}:".format(snid)
    print "\t\tMeta information:"
    meta = prob['meta']
    for key in meta.array.names:
      print "\t\t\t{}: {}".format(key, meta[key])
    print "\t\tCorrect: {}".format(prob['correct'])
    print "\t\tOrdered: {}".format(prob['ordered'])
    """
    print "\t\tResults Bundle:"
    for key, val in prob['results'].iteritems():
      print "\t\t\t{}: {}".format(key, val)
    """
    print "\t\tOther information:"
    for key, val in prob.iteritems():
      if key == 'meta' or key == 'correct' or key == 'ordered' or key == 'results':
        continue
      print "\t\t\t{}: {}".format(key, val)

def print_snids(datas):
  for snid, val, in datas.iteritems():
    print snid

def count_zero_prob(probs):
  zeros = 0
  for snid, prob in probs.iteritems():
    for key, val in prob.iteritems():
      if key == 'meta' or key == 'correct' or key == 'ordered' or key == 'results':
        continue
      if val == 0.0:
        zeros += 1
        break
  print "Got {} / {} zeros".format(zeros, len(probs))

def load_snid_data(snid, directory):
  metas, datas = open_sim_fits(directory)
  meta = [m for m in metas if int(m['SNID']) == int(snid)][0]
  return datas[(meta['PTROBS_MIN'] - 1):meta['PTROBS_MAX']]

def plot_lc(snid, directory, show, outname, result, mms):
  bestname = mms[0]
  ms = [m for m in mms if m in models] # used to eliminate errors
  specialm = ms[1] # this is the model name of the model that was used to simulate the SN
  print "Plotting {} on {}{}".format(ms, snid, (" instead of " + str(mms) + " because the required models aren't loaded" if len(mms) != len(ms) else ""))
  mods = [(m, models[m]['model']) for m in ms]
  dat = load_snid_data(snid, directory)
  dat['FLT'][:] = np.char.strip(dat['FLT'])
  data = Table(dat)
  data.rename_column('FLUXCAL', 'flux')
  data.rename_column('FLUXCALERR', 'fluxerr')
  data.rename_column('FLT', 'band')
  data.rename_column('MJD', 'time')
  data['zp'] = 27.5
  data['zpsys'] = 'ab'
  for m, mod in mods:
    if m == specialm:
      res = result[m]
      param_dict = {'hostebv': res['param_dict']['hostebv'], 'z': result['meta']['SIM_REDSHIFT_CMB'], 'mwebv': result['meta']['SIM_MWEBV'], 't0': result['meta']['SIM_PEAKMJD'], 'amplitude': res['param_dict']['amplitude']}
      mod.set(**param_dict)
    else:
      res = result[m]
      mod.set(**res['param_dict'])
  fig = sncosmo.plot_lc(data, [mod for m, mod in mods])
  if show:
    plt.show()
  else:
    fig.savefig(outname, dpi=300)

def plot_lcs(probs, data, figuresDirectory, directory):
  if not os.path.exists(figuresDirectory):
    os.mkdir(figuresDirectory)
  snids = [snid for snid, val in probs.iteritems()]
  for snid in snids:
    result = [dat for dat in data if int(dat['meta']['SNID']) == int(snid)][0]
    plot_lc(snid, directory, False, join(figuresDirectory, str(snid) + '.png'), result, [probs[snid]['ordered'][0][0], SNANAidx_to_model(probs[snid]['meta']['SIM_NON1a'])])

def plot_corner(prob, figuresDirectory, model, tag=None):
  if not os.path.exists(figuresDirectory):
    os.mkdir(figuresDirectory)
  #model = prob['ordered'][0][0]
  result = prob['results'][model]
  snid = prob['meta']['SNID']
  weights = None
  if not _isMCMC:
    weights = result['weights']
    if np.max(weights) >= 0.9999999999:
      print "Failed to create corner plots for {} - max(weights) >= 0.9999999999".format(snid)
      return
  samples = result['samples']
  param_names = result['vparam_names'] # this will break when ran with data created by sncosmo.__version__ < 1
  extents = len(param_names) * [0.9999999999]
  nbins = 15
  if not _isMCMC:
    fig = corner(samples, labels=param_names, weights=weights, extents=extents, bins=nbins)
  else:
    fig = corner(samples, labels=param_names, extents=extents, bins=nbins)

  fig.savefig(join(figuresDirectory, str(snid) + (tag + '.png' if tag != None else '.png')))

def plot_corners(prob, figuresDirectory):
  if not os.path.exists(figuresDirectory):
    os.mkdir(figuresDirectory)
  if not _isMCMC:
    for pair in prob['ordered']:
      model = pair[0]
      plot_corner(prob, figuresDirectory, model, '_' + model + '_' + str(pair[1]))
  else:
    for model in prob['results'].keys():
      if model != 'meta':
        plot_corner(prob, figuresDirectory, model, '_' + model)

def main(args):
  global _cacheDirectory, _isMCMC

  show = False
  outname = 'out.png'
  
  conf = {}
  if isfile('config.py'):
    from config import config
    conf = config
  directory = os.path.expanduser(conf['fitsDirectory']) if 'fitsDirectory' in conf else '/home/kuhlmann/snana/root_v201204/SIM/DES_5years_CC_v1033f/'
  _cacheDirectory = conf['cacheDirectory'] if 'cacheDirectory' in conf else _cacheDirectory
  _isMCMC = conf['mcmc'] if 'mcmc' in conf else _isMCMC
  
  parser = argparse.ArgumentParser(description='Analyze SN Data Simulated from SNANA and fit with SNCosmo')
  parser.add_argument('-x', '--extra', action='store_const', const=False, default=True, help="exclude extra templates in analysis")
  parser.add_argument('-s', '--show', action='store_const', const=True, default=False, help="show plot instead of saving")
  parser.add_argument('-o', '--out', nargs=1, default='out.png', help='name of png to create')
  opts = parser.parse_args(args)
  include_special = opts.extra
  show = opts.show
  outname = opts.out[0]
  dirout = outname if not outname.endswith('.png') else 'figures'
  data = load_data()
  chisqdofs, probs, lowestchisqdofs, lowestcorrect, lowestincorrect = analyze_data(data, include_special)
  #lowprobs = filter_probabilities(probs, False, True, True, 0.1, condition=lambda bundle: bundle['meta']['SIM_NON1a'] == 104)

  #aprob = {snid: p for snid, p in probs.iteritems() if int(snid) == 338990}
  #print_prob_info(aprob)
  #plot_corners([p for snid, p in aprob.iteritems()][0], dirout)
  #plot_lcs(aprob, data, dirout, directory)
  plot_types(lowestchisqdofs, show, outname)

  
if __name__ == "__main__":
  main(sys.argv[1:])

