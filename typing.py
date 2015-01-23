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

_cacheDirectory = "./.cache"

def chisqdof_to_prob(chisq, dof):
  return 1 - stats.chi2.cdf(chisq, dof)

def sntypes_equiv(typea, typeb):
  typea = typea.replace('SN ', '')
  typeb = typeb.replace('SN ', '')
  if typea == 'IIL/P' and (typeb == 'IIP' or typeb == 'IIL'):
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
    correctlytyped = sntypes_equiv(realtype, lowchi[0])
    chisqbundle['correct'] = probbundle['correct'] = correctlytyped
    chisqbundle['ordered'] = orderedchis
    probbundle['ordered'] = orderedprobs
    chisqbundle['meta'] = probbundle['meta'] = meta
    chisqdofs[meta['SNID']] = chisqbundle
    probs[meta['SNID']] = probbundle
    if lowchi[0] in lowestchisqdofs:
      lowestchisqdofs[lowchi[0]] += 1
    else:
      lowestchisqdofs[lowchi[0]] = 1
    adder = lowestcorrect if correctlytyped else lowestincorrect
    if lowchi[0] in adder:
      adder[lowchi[0]] += 1
    else:
      adder[lowchi[0]] = 1

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

def filter_probabilities(probs, corrects, incorrects, greater, cutoff):
  cutter = (lambda x: x >= cutoff) if greater else (lambda y: y <= cutoff)
  retval = {}
  for snid, prob in probs.iteritems():
    if not cutter(prob['ordered'][0][1] - prob['ordered'][1][1]):
      continue
    if prob['correct'] and corrects:
      retval[snid] = prob
    elif not prob['correct'] and incorrects:
      retval[snid] = prob
  return retval

def print_prob_info(probs):
  for snid, prob in probs.iteritems():
    print "\tLooking at {}:".format(snid)
    #print "\t\tMeta information:"
    #meta = prob['meta']
    #for key in meta.array.names:
    #  print "\t\t\t{}: {}".format(key, meta[key])
    print "\t\tOther information:"
    for key, val in prob.iteritems():
      if key == 'meta':
        continue
      print "\t\t\t{}: {}".format(key, val)


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
  outname = opts.out[0]
  data = load_data()
  chisqdofs, probs, lowestchisqdofs, lowestcorrect, lowestincorrect = analyze_data(data, include_special)
  #plot_types(lowestchisqdofs, show, outname)
  #plot_types(lowestcorrect, show, outname)
  #plot_types(lowestincorrect, show, outname)
  #histo_probs(probs, show, outname)
  #histo_chisqdiff(chisqdofs, show, outname)
  #histo_probdiff(probs, show, outname)
  lowprobs = filter_probabilities(probs, False, True, True, 0.5)
  #histo_probdiff(lowprobs, show, outname)
  #print "Probability info for SN with false typing with diff > 50%:"
  print_prob_info(probs)

if __name__ == "__main__":
  main(sys.argv[1:])

