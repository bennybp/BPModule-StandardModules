#pragma once

#include "Integrals/boys/Boys_longfac.hpp"
#include "Integrals/boys/Boys_shortgrid.hpp"

#include <cmath>

namespace psr_modules {
namespace integrals {
namespace detail {

/*! \brief Calculate the value of the Boys function for a small value of x
 *         via a Taylor series
 *
 *  The Boys function for all [0, n] are evaluated via the Taylor series
 *
 * \warning The output buffer \p F must be large enough to hold \p n + 1 elements,
 *          since the output includes F(0, x) and F(n, x)

 * \param [out] F Output buffer. Size must be n+1 elements
 * \param [in] n Maximum order of the Boys function
 * \param [in] x The value of x for the Boys function
 */
inline void boys_f_taylor(double * const RESTRICT F, int n, double x)
{
    // calculate the index of the nearest grid point
    const int lookup_idx = (int)((1.0/PSR_MODULES_BOYS_SHORTGRID_SPACE)*(x+(PSR_MODULES_BOYS_SHORTGRID_SPACE/2)));

    // the value of x at the nearest grid point
    const double xi = ((double)lookup_idx * PSR_MODULES_BOYS_SHORTGRID_SPACE);

    // distance from the nearest grid point (with negative sign)
    const double dx = xi-x;   // -delta x

    // Get all the values of of F(n,x) for all n at the nearest grid point
    double const * const RESTRICT gridpts = &(lut::boys_shortgrid[lookup_idx][0]);

    for(int i = 0; i <= n; ++i)
    {
        const double f0xi = gridpts[i];
        const double f1xi = gridpts[i+1];
        const double f2xi = gridpts[i+2];
        const double f3xi = gridpts[i+3];
        const double f4xi = gridpts[i+4];
        const double f5xi = gridpts[i+5];
        const double f6xi = gridpts[i+6];
        const double f7xi = gridpts[i+7];

        F[i] = f0xi
               + dx * (                  f1xi
               + dx * ( (1.0/2.0   )   * f2xi
               + dx * ( (1.0/6.0   )   * f3xi
               + dx * ( (1.0/24.0  )   * f4xi
               + dx * ( (1.0/120.0 )   * f5xi
               + dx * ( (1.0/720.0 )   * f6xi
               + dx * ( (1.0/5040.0)   * f7xi
               )))))));
    }
}


/*! \brief Calculate the value of the Boys function for a large value of x
 *         via an asymptotic equation
 *
 *  The Boys function for all [0, n] are evaluated via
 *
 *  \f[
 *   F(n, x) \approx \frac{(2n-1)!!}{2^{n+1}} \sqrt{\frac{\pi}{x^{2n+1}}}
 *  \f]
 *
 * \warning The output buffer \p F must be large enough to hold \p n + 1 elements,
 *          since the output includes F(0, x) and F(n, x)

 * \param [out] F Output buffer. Size must be n+1 elements
 * \param [in] n Maximum order of the Boys function
 * \param [in] x The value of x for the Boys function
 */
inline void boys_f_long(double * const RESTRICT F, int n, double x)
{
    
    // boys_longfac[n] = sqrt(pi)*(2n-1)!!/(2**(n+1))
    const double x1 = 1.0/x;
    double x2 = sqrt(x1);

    for(int i = 0; i <= n; i++)
    {
        F[i] = lut::boys_longfac[i] * x2;
        x2 *= x1;
    }
}



/*! \brief Calculate the value of the Boys function for a given value
 *
 * Depending on the magnitude of \p x, the value may be calculated via a 
 * Taylor series or by an asymptotic equation

 * \warning The output buffer \p F must be large enough to hold \p n + 1 elements,
 *          since the output includes F(0, x) and F(n, x)

 * \param [out] F Output buffer. Size must be n+1 elements
 * \param [in] n Maximum order of the Boys function
 * \param [in] x The value of x for the Boys function
 */
inline void calculate_f(double * const RESTRICT F, int n, double x)
{
    if(x < PSR_MODULES_BOYS_SHORTGRID_MAXX)
        detail::boys_f_taylor(F, n, x);
    else
        detail::boys_f_long(F, n, x);
}


} // close namespace detail
} // close namespace integrals
} // close namespace psr_modules

