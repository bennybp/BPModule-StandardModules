#ifndef STANDARD_COMMON_BOYS_HPP
#define STANDARD_COMMON_BOYS_HPP

#include <cmath>

////////////////////////////////////////////////////////////////////////
// Info for boys function of large x (asymptotic equation)
// stored in Boys.cpp
#define BOYS_LONGFAC_MAXN 31  // store [0, BOYS_LONGFAC_MAXN] inclusive
extern double const boys_longfac[BOYS_LONGFAC_MAXN+1];
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Info for boys function of short x (taylor series)
#define BOYS_SHORTGRID_MAXN 31
#define BOYS_SHORTGRID_MAXX 43.0
#define BOYS_SHORTGRID_SPACE 0.1
#define BOYS_SHORTGRID_NPOINT 431
#define BOYS_SHORTGRID_LOOKUPFAC 10.0
#define BOYS_SHORTGRID_LOOKUPFAC2 0.05
extern double boys_shortgrid[BOYS_SHORTGRID_NPOINT][BOYS_SHORTGRID_MAXN+1];
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////
// Calculate the value of the boys function for a small value of x
// via a taylor series
////////////////////////////////////////////////////////////////////
inline void Boys_F_taylor(double * const RESTRICT F, int n, double x)
{
    const int lookup_idx = (int)(BOYS_SHORTGRID_LOOKUPFAC*(x+BOYS_SHORTGRID_LOOKUPFAC2));
    const double xi = ((double)lookup_idx * BOYS_SHORTGRID_SPACE);
    const double dx = xi-x;   // -delta x

    double const * const RESTRICT gridpts = &(boys_shortgrid[lookup_idx][0]);

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


////////////////////////////////////////////////////////////////////
// Calculate the value of the boys function for a small value of x
// via the asymptotic equation
////////////////////////////////////////////////////////////////////
inline void Boys_F_long(double * const RESTRICT F, int n, double x)
{
    const double x1 = 1.0/x;
    double x2 = sqrt(x1);

    for(int i = 0; i <= n; i++)
    {
        F[i] = boys_longfac[i] * x2;
        x2 *= x1;
    }
}


///////////////////////////////////////////////////////////////////
// Main interface for calculation of the boys function
///////////////////////////////////////////////////////////////////
inline void CalculateF(double * const RESTRICT F, int n, double x)
{
    if(x < BOYS_SHORTGRID_MAXX)
        Boys_F_taylor(F, n, x);
    else
        Boys_F_long(F, n, x);
}

#endif
