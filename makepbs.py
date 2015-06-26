#!/usr/bin/env python

import sys
from os import mkdir
from os.path import join, exists
from snutils import load_summary, open_sim_fits
import argparse

def make(nnodes, ncores, pbsfile, cmdlist, workingDirectory, logDirectory, cache, summaryFile, scriptName, jobName, projName, walltime):
  ncomm = nnodes*ncores
  processors = 'nodes={}:ppn={}'.format(nnodes, ncores)

  if not exists(join(workingDirectory, cache)):
    mkdir(join(workingDirectory, cache))
  if not exists(join(workingDirectory, logDirectory)):
    mkdir(join(workingDirectory, logDirectory))

  nsnids = len(load_summary(summaryFile))
  nper = nsnids / ncomm
  nleft = nsnids % ncomm
  current = 0
  f = open(join(workingDirectory, cmdlist), 'w')

  for i in range(ncomm):
    ni = nper
    if nleft > 0:
      nleft -= 1
      ni += 1
    f.write('{} -f {} -t {} 2>&1 | tee {}/output_{}-{}.log\n'.format(join(workingDirectory, scriptName), current, current + ni, join(workingDirectory, logDirectory), current, current + ni))
    current += ni

  f.close()

  if not exists(join(workingDirectory, logDirectory, 'output')):
    mkdir(join(workingDirectory, logDirectory, 'output'))
  if not exists(join(workingDirectory, logDirectory, 'error')):
    mkdir(join(workingDirectory, logDirectory, 'error'))
  
  f = open(pbsfile, 'w')
  command = 'mpirun -n {} minions < {}'.format(ncomm, join(workingDirectory, cmdlist))
  f.write("""#!/bin/sh
#PBS -N {}
#PBS -A {}
#PBS -q shared
#PBS -m bea
#PBS -l walltime={}
#PBS -l {}
#PBS -o {}
#PBS -e {}
cd {}
{}""".format(jobName, projName, walltime, processors, join(workingDirectory, logDirectory, 'output', jobName + '.out'), join(workingDirectory, logDirectory, 'error', jobName + '.err'), workingDirectory, command))
  f.close()

def parseArgs(args):
  nnodes = 3
  ncores = 8

  pbsfile = '.fitting.pbs'
  cmdlist = '.fitting.list'
  workingDirectory = '/lcrc/project/Supernova_DES/sim_typing/'
  cache = './.cache/'
  logDirectory = './.cache/.logs'
  summaryFile = './PassAllCuts500.summary'
  scriptName = './fitting.py'
  jobName = 'fitting'
  projName = 'Supernova_DES'
  walltime = '00:55:00'
  
  return [
      nnodes, ncores, pbsfile, cmdlist, workingDirectory, logDirectory,
      cache, summaryFile, scriptName, jobName, projName, walltime
  ]

def main(args):
  ar = parseArgs(args)
  make(*ar)

if __name__ == "__main__":
  main(sys.argv[1:])

