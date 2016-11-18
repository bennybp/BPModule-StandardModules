#!/usr/bin/env python3

#######################################
# Generates an array for the constant
# used in the long-range formula
#
# Eqn. 9.8.9 from Helgaker et al
# 
# Fn(x) = (2n-1)!!/(2**(n+1)) sqrt(pi/(x**(2n+1)))
#
# So we can pre-compute (2n-1)!!/(2**(n+1)) sqrt(pi) for
# all n that we might need
# After that, Fn(x) = longfac[n] * 1/sqrt(x**(2n+1))
#######################################

import argparse
import sys
from mpmath import mp # arbitrary-precision math

# Note that types are being stored as a string. This is so mpmath
# can parse it without turning it into a (possibly) lesser-precision float 
parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True,               help="Output file name base (no extension)")
parser.add_argument("--max-n",    type=int, required=True,               help="Maximum n value to go to")
parser.add_argument("--dps",      type=int, required=False, default=256, help="Decimal precision/sig figs to use/calculate")
args = parser.parse_args()

# Set the dps option
mp.dps = args.dps

# Convert stuff to mpmath
maxn = args.max_n

# Generate fac = sqrt(pi)*(2n-1)!!/(2**(n+1))
longfac = [None] * (maxn+1)
longfac[0] = mp.sqrt(mp.pi)/mp.mpf(2)

# each iteration just adds a 1/2 and a (2n-1)
for n in range(1, maxn+1):
 longfac[n] = longfac[n-1] * 0.5 * (2.0*n - 1)

# Output to file
with open(args.filename + ".cpp", 'w') as f:
  f.write("/*\n")
  f.write(" Generated with:\n")
  f.write("   " + " ".join(sys.argv[:]))
  f.write("\n")
  f.write("------------------------------------\n")
  f.write("Options for asymptotic Boys function factor sqrt(pi)*(2n-1)!!/(2**(n+1)):\n")
  f.write("    Max n: {}\n".format(maxn))
  f.write("      DPS: {}\n".format(args.dps))
  f.write("------------------------------------\n")
  f.write("*/\n\n")

  f.write("namespace psr_modules {\n")
  f.write("namespace integrals {\n")
  f.write("namespace lut {\n")
  f.write("\n")

  f.write("extern const double boys_longfac[{}] = \n".format(maxn+1))
  f.write("{\n")

  for i,n in enumerate(longfac):
      f.write("/* n = {:4} */  {:32},\n".format(i, mp.nstr(n, 18)))
  f.write("};\n")

  f.write("\n")
  f.write("} // closing namespace lut\n")
  f.write("} // closing namespace integrals\n")
  f.write("} // closing namespace psr_modules\n")
  f.write("\n")

with open(args.filename + ".hpp", 'w') as f: 
  f.write("#pragma once\n")
  f.write("\n")
  f.write("#define PSR_MODULES_BOYS_LONGFAC_MAXN {}\n".format(maxn))
  f.write("\n")
  f.write("namespace psr_modules {\n")
  f.write("namespace integrals {\n")
  f.write("namespace lut {\n")

  f.write("\n")
  f.write("/*! Array of factors in the asymptotic Boys function approximation\n")
  f.write(" *\n")
  f.write(" * boys_longfac[n] = sqrt(pi)*(2n-1)!!/(2**(n+1))\n")
  f.write(" *\n")
  f.write(" * Then, evaluating the asymptotic Boys function is easy:\n")
  f.write(" *\n")
  f.write(" * F(n,x) = boys_longfac[n] * sqrt(1.0/(x^(2n+1)))\n")
  f.write(" */\n")
  f.write("extern const double boys_longfac[PSR_MODULES_BOYS_LONGFAC_MAXN+1];\n")
  f.write("\n")

  f.write("} // closing namespace lut\n")
  f.write("} // closing namespace integrals\n")
  f.write("} // closing namespace psr_modules\n")
  f.write("\n")
