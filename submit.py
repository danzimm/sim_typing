#!/usr/bin/env python

import makepbs
import sys
from subprocess import call
import time

def main(args):
  ar = makepbs.parseArgs(args)
  makepbs.make(*ar)
  #time.sleep(1)
  call(["qsub", str(ar[2])])

if __name__ == "__main__":
  main(sys.argv[1:])

