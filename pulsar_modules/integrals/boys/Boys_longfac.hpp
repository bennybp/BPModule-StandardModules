#pragma once

#define PSR_MODULES_BOYS_LONGFAC_MAXN 40

namespace psr_modules {
namespace integrals {
namespace lut {

/*! Array of factors in the asymptotic Boys function approximation
 *
 * boys_longfac[n] = sqrt(pi)*(2n-1)!!/(2**(n+1))
 *
 * Then, evaluating the asymptotic Boys function is easy:
 *
 * F(n,x) = boys_longfac[n] * sqrt(1.0/(x^(2n+1)))
 */
extern const double boys_longfac[PSR_MODULES_BOYS_LONGFAC_MAXN+1];

} // closing namespace lut
} // closing namespace integrals
} // closing namespace psr_modules

