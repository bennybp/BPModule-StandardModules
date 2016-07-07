#pragma once

#define PSR_MODULES_BOYS_SHORTGRID_MAXN 40
#define PSR_MODULES_BOYS_SHORTGRID_MAXX 43.0
#define PSR_MODULES_BOYS_SHORTGRID_SPACE 0.1
#define PSR_MODULES_BOYS_SHORTGRID_NPOINT 431
#define PSR_MODULES_BOYS_SHORTGRID_LOOKUPFAC 10.0
#define PSR_MODULES_BOYS_SHORTGRID_LOOKUPFAC2 0.05

namespace psr_modules {
namespace integrals {
namespace lut {

/* A grid of precomputed values of the Boys function */
extern const double boys_shortgrid[PSR_MODULES_BOYS_SHORTGRID_NPOINT][PSR_MODULES_BOYS_SHORTGRID_MAXN+1];

} // closing namespace lut
} // closing namespace integrals
} // closing namespace psr_modules

