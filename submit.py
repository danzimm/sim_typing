#!/usr/bin/env python

import makepbs
import sys
from subprocess import call

def main(args):
  ar = makepbs.parseArgs(args)
  makepbs.make(*ar)
  call(["qsub", str(ar[2])])

if __name__ == "__main__":
  main(sys.argv[1:])

