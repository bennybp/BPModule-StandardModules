#!/usr/bin/env python3

#######################################
# Generates grid for the
# boys function to arbitrary precision
#######################################

import argparse
import sys
from mpmath import mp # arbitrary-precision math


# A slow but accurate way of calculating the boys function
# via the incomplete gamma function
def BoysValue(n, x):
    F = []

    if x == mp.mpf("0"):
        for i in range(0, n+1):
            F.append(mp.mpf(1.0)/(mp.mpf(2.0*i+1)))
    else:
        for i in range(0, n+1):
            N = i+mp.mpf("0.5")
            F.append(mp.gammainc(N, 0, x) * 1.0/(2.0 * mp.power(x, N)))

    return F

                                                                         

# Note that types are being stored as a string. This is so mpmath
# can parse it without turning it into a (possibly) lesser-precision float 
parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True,               help="Output file name base (no extension)")
parser.add_argument("--max-n",    type=int, required=True,               help="Maximum n value to go to")
parser.add_argument("--max-x",    type=str, required=True,               help="Cutoff for the x value")
parser.add_argument("--spacing",  type=str, required=True,               help="Space between pre-computed points")
parser.add_argument("--dps",      type=int, required=False, default=256, help="Decimal precision/sig figs to use/calculate")
parser.add_argument("--eps",      type=int, required=False, default=32,  help="Overall precision / accuracy")
args = parser.parse_args()

# Set the dps option
mp.dps = args.dps

# Convert stuff to mpmath
inc = mp.mpf(args.spacing)
maxx = mp.mpf(args.max_x)
eps = mp.power(10, -args.eps)
maxn = args.max_n
npoints = int(maxx / inc) + 1

# Start at x=0 and increment up
x = mp.mpf(0)
F = []
pts = []

while x < maxx or (x-maxx < inc):
  F2 = [None]*(maxn+1)  # Creates a list of maxn+1 elements

  x2 = 2*x
  ex = mp.exp(-x)
 
  F2 = BoysValue(maxn, x) 

  F.append(F2)
  pts.append(x)
  x += inc

# Output to file
with open(args.filename + ".cpp", 'w') as f:
  f.write("/*\n")
  f.write(" Generated with:\n")
  f.write("   " + " ".join(sys.argv[:]))
  f.write("\n")
  f.write("------------------------------------\n")
  f.write("Options for Boys function Fn(x):\n")
  f.write("    Max n: {}\n".format(maxn))
  f.write("    Max x: {}\n".format(maxx))
  f.write("  Spacing: {}\n".format(inc))
  f.write("  npoints: {}\n".format(npoints))
  f.write("      DPS: {}\n".format(args.dps))
  f.write("      EPS: 10^({})\n".format(-args.eps))
  f.write("         = {}\n".format(eps))
  f.write("------------------------------------\n")
  f.write("*/\n\n")

  f.write("namespace psr_modules {\n")
  f.write("namespace integrals {\n")
  f.write("namespace lut {\n")

  f.write("\n")
  f.write("extern const double boys_shortgrid[{}][{}] = \n".format(npoints, maxn+1))
  f.write("{\n")

  for p,x in zip(F,pts):
    f.write("/* x = {:12}*/  {{".format(mp.nstr(x, 4)))
    for n in p:
      f.write("{:32}, ".format(mp.nstr(n, 18)))
    f.write("},\n")
  f.write("};\n")

  f.write("\n")
  f.write("} // closing namespace lut\n")
  f.write("} // closing namespace integrals\n")
  f.write("} // closing namespace psr_modules\n")
  f.write("\n")

with open(args.filename + ".hpp", 'w') as f: 
  f.write("#pragma once\n")
  f.write("\n")
  f.write("#define PSR_MODULES_BOYS_SHORTGRID_MAXN {}\n".format(maxn))
  f.write("#define PSR_MODULES_BOYS_SHORTGRID_MAXX {}\n".format(maxx))
  f.write("#define PSR_MODULES_BOYS_SHORTGRID_SPACE {}\n".format(inc))
  f.write("#define PSR_MODULES_BOYS_SHORTGRID_NPOINT {}\n".format(npoints))
  f.write("#define PSR_MODULES_BOYS_SHORTGRID_LOOKUPFAC {}\n".format(1.0/inc))
  f.write("#define PSR_MODULES_BOYS_SHORTGRID_LOOKUPFAC2 {}\n".format(0.5*inc))
  f.write("\n")

  f.write("namespace psr_modules {\n")
  f.write("namespace integrals {\n")
  f.write("namespace lut {\n")
  f.write("\n")

  f.write("/* A grid of precomputed values of the Boys function */\n")
  f.write("extern const double boys_shortgrid[PSR_MODULES_BOYS_SHORTGRID_NPOINT][PSR_MODULES_BOYS_SHORTGRID_MAXN+1];\n")

  f.write("\n")
  f.write("} // closing namespace lut\n")
  f.write("} // closing namespace integrals\n")
  f.write("} // closing namespace psr_modules\n")
  f.write("\n")
