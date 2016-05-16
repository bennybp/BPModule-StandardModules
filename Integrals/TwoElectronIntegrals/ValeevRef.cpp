///////////////////////////////////////////
// Reference calculation of ERI
// Adapted from libint 2.0.5, by E. Valeev
// https://github.com/evaleev/libint
// http://www.valeevgroup.chem.vt.edu
//
// Originally released under GPL v2
///////////////////////////////////////////
//
// LIBINT (version 2) - a library for the evaluation of molecular integrals of many-body
// operators over Gaussian functions
// Copyright (C) 2004-2013 Edward F. Valeev
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program (see file LICENSE); if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
///////////////////////////////////////////


#include <cmath>
#include <pulsar/math/Factorial.hpp>
#include <pulsar/math/Binomial.hpp>
#include <pulsar/constants.h>


#define EPS 1.0e-17
#define MAXFAC 100

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


using namespace pulsar::math;


static double* init_array(unsigned long int size)
{
    double* result = (double *)malloc(size * sizeof(double));
    for (unsigned long int i = 0; i < size; i++)
        result[i] = 0.0;
    return result;
}

static void free_array(double* array)
{
    free(array);
}


static void Valeev_F(double *F, int n, double x)
{
    int i, m;
    int m2;
    double t2;
    double num;
    double sum;
    double term1;
    const double K = 0.8862269254527580136490837416705725913987747280611935641069038949264556422955160906874753283692723327l;
    double et;


    if (x > 20.0)   /* For big t's do upward recursion */
    {
        t2 = 2 * x;
        et = exp(-x);
        x = sqrt(x);
        F[0] = K * erf(x) / x;
        for (m = 0; m <= n - 1; m++)
        {
            F[m + 1] = ((2 * m + 1) * F[m] - et) / (t2);
        }
    }
    else
    {
        /* For smaller t's compute F with highest n using
         asymptotic series (see I. Shavitt in
         Methods in Computational Physics, ed. B. Alder eta l,
         vol 2, 1963, page 8) */
        et = exp(-x);
        t2 = 2 * x;
        m2 = 2 * n;
        num = DoubleFactorialD(m2-1);
        i = 0;
        sum = 1.0 / (m2 + 1);
        do
        {
            i++;
            num = num * t2;
            term1 = num / DoubleFactorialD(m2 + 2 * i + 1);
            sum += term1;
        }
        while (fabsl(term1) > EPS && i < MAXFAC);
        F[n] = sum * et;
        for (m = n - 1; m >= 0; m--)   /* And then do downward recursion */
        {
            F[m] = (t2 * F[m + 1] + et) / (2 * m + 1);
        }
    }
}

double ValeevRef_eri(int l1, int m1, int n1, double alpha1, const double* A,
                     int l2, int m2, int n2, double alpha2, const double* B,
                     int l3, int m3, int n3, double alpha3, const double* C,
                     int l4, int m4, int n4, double alpha4, const double* D)
{

    const double gammap = alpha1 + alpha2;
    const double Px = (alpha1 * A[0] + alpha2 * B[0]) / gammap;
    const double Py = (alpha1 * A[1] + alpha2 * B[1]) / gammap;
    const double Pz = (alpha1 * A[2] + alpha2 * B[2]) / gammap;
    const double PAx = Px - A[0];
    const double PAy = Py - A[1];
    const double PAz = Pz - A[2];
    const double PBx = Px - B[0];
    const double PBy = Py - B[1];
    const double PBz = Pz - B[2];
    const double AB2 = (A[0] - B[0]) * (A[0] - B[0]) + (A[1] - B[1]) * (A[1]
                       - B[1]) + (A[2] - B[2]) * (A[2] - B[2]);

    const double gammaq = alpha3 + alpha4;
    const double gammapq = gammap * gammaq / (gammap + gammaq);
    const double Qx = (alpha3 * C[0] + alpha4 * D[0]) / gammaq;
    const double Qy = (alpha3 * C[1] + alpha4 * D[1]) / gammaq;
    const double Qz = (alpha3 * C[2] + alpha4 * D[2]) / gammaq;
    const double QCx = Qx - C[0];
    const double QCy = Qy - C[1];
    const double QCz = Qz - C[2];
    const double QDx = Qx - D[0];
    const double QDy = Qy - D[1];
    const double QDz = Qz - D[2];
    const double CD2 = (C[0] - D[0]) * (C[0] - D[0]) + (C[1] - D[1]) * (C[1]
                       - D[1]) + (C[2] - D[2]) * (C[2] - D[2]);

    const double PQx = Px - Qx;
    const double PQy = Py - Qy;
    const double PQz = Pz - Qz;
    const double PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;

    int u1, u2, v1, v2, w1, w2, tx, ty, tz, txmax, tymax, tzmax;
    int i, j, k;
    int lp, lq, mp, mq, np, nq;
    int zeta;
    double *flp, *flq, *fmp, *fmq, *fnp, *fnq;
    double *F;
    double K1, K2;
    double Gx, Gy, Gz;
    double pfac;
    double result = 0.0;
    double tmp;
    int u1max, u2max, v1max, v2max, w1max, w2max;

    K1 = exp(-alpha1 * alpha2 * AB2 / gammap);
    K2 = exp(-alpha3 * alpha4 * CD2 / gammaq);
    pfac = 2 * pow(PI, 2.5) * K1 * K2 / (gammap * gammaq
            * sqrt(gammap + gammaq));

    F = init_array(l1 + l2 + l3 + l4 + m1 + m2 + m3 + m4 + n1 + n2 + n3 + n4 + 1);
    Valeev_F(F, l1 + l2 + l3 + l4 + m1 + m2 + m3 + m4 + n1 + n2 + n3 + n4,
           PQ2 * gammapq);

    flp = init_array(l1 + l2 + 1);
    for (k = 0; k <= l1 + l2; k++)
        for (i = 0; i <= MIN(k,l1); i++)
        {
            j = k - i;
            if (j > l2)
                continue;
            tmp = BinomialCoefficient(l1,i) * BinomialCoefficient(l2,j);
            if (l1 - i > 0)
                tmp *= pow(PAx, l1 - i);
            if (l2 - j > 0)
                tmp *= pow(PBx, l2 - j);
            flp[k] += tmp;
        }
    fmp = init_array(m1 + m2 + 1);
    for (k = 0; k <= m1 + m2; k++)
        for (i = 0; i <= MIN(k,m1); i++)
        {
            j = k - i;
            if (j > m2)
                continue;
            tmp = BinomialCoefficient(m1,i) * BinomialCoefficient(m2,j);
            if (m1 - i > 0)
                tmp *= pow(PAy, m1 - i);
            if (m2 - j > 0)
                tmp *= pow(PBy, m2 - j);
            fmp[k] += tmp;
        }
    fnp = init_array(n1 + n2 + 1);
    for (k = 0; k <= n1 + n2; k++)
        for (i = 0; i <= MIN(k,n1); i++)
        {
            j = k - i;
            if (j > n2)
                continue;
            tmp = BinomialCoefficient(n1,i) * BinomialCoefficient(n2,j);
            if (n1 - i > 0)
                tmp *= pow(PAz, n1 - i);
            if (n2 - j > 0)
                tmp *= pow(PBz, n2 - j);
            fnp[k] += tmp;
        }
    flq = init_array(l3 + l4 + 1);
    for (k = 0; k <= l3 + l4; k++)
        for (i = 0; i <= MIN(k,l3); i++)
        {
            j = k - i;
            if (j > l4)
                continue;
            tmp = BinomialCoefficient(l3,i) * BinomialCoefficient(l4,j);
            if (l3 - i > 0)
                tmp *= pow(QCx, l3 - i);
            if (l4 - j > 0)
                tmp *= pow(QDx, l4 - j);
            flq[k] += tmp;
        }
    fmq = init_array(m3 + m4 + 1);
    for (k = 0; k <= m3 + m4; k++)
        for (i = 0; i <= MIN(k,m3); i++)
        {
            j = k - i;
            if (j > m4)
                continue;
            tmp = BinomialCoefficient(m3,i) * BinomialCoefficient(m4,j);
            if (m3 - i > 0)
                tmp *= pow(QCy, m3 - i);
            if (m4 - j > 0)
                tmp *= pow(QDy, m4 - j);
            fmq[k] += tmp;
        }
    fnq = init_array(n3 + n4 + 1);
    for (k = 0; k <= n3 + n4; k++)
        for (i = 0; i <= MIN(k,n3); i++)
        {
            j = k - i;
            if (j > n4)
                continue;
            tmp = BinomialCoefficient(n3,i) * BinomialCoefficient(n4,j);
            if (n3 - i > 0)
                tmp *= pow(QCz, n3 - i);
            if (n4 - j > 0)
                tmp *= pow(QDz, n4 - j);
            fnq[k] += tmp;
        }

    for (lp = 0; lp <= l1 + l2; lp++)
        for (lq = 0; lq <= l3 + l4; lq++)
        {
            u1max = lp / 2;
            u2max = lq / 2;
            for (u1 = 0; u1 <= u1max; u1++)
                for (u2 = 0; u2 <= u2max; u2++)
                {
                    Gx = pow(-1, lp) * flp[lp] * flq[lq] * FactorialD(lp) * FactorialD(lq)
                         * pow(gammap, u1 - lp) * pow(gammaq, u2 - lq) * FactorialD(lp + lq - 2 * u1 - 2 * u2) * pow(gammapq, lp + lq - 2 * u1 - 2 * u2)
                         / (FactorialD(u1) * FactorialD(u2) * FactorialD(lp - 2 * u1) * FactorialD(lq - 2 * u2));
                    for (mp = 0; mp <= m1 + m2; mp++)
                        for (mq = 0; mq <= m3 + m4; mq++)
                        {
                            v1max = mp / 2;
                            v2max = mq / 2;
                            for (v1 = 0; v1 <= v1max; v1++)
                                for (v2 = 0; v2 <= v2max; v2++)
                                {
                                    Gy = pow(-1, mp) * fmp[mp] * fmq[mq] * FactorialD(mp) * FactorialD(mq)
                                         * pow(gammap, v1 - mp) * pow(gammaq, v2 - mq) * FactorialD(mp + mq - 2 * v1 - 2 * v2) * pow(gammapq,
                                                         mp + mq - 2 * v1 - 2 * v2)
                                         / (FactorialD(v1) * FactorialD(v2) * FactorialD(mp - 2 * v1)
                                            * FactorialD(mq - 2 * v2));
                                    for (np = 0; np <= n1 + n2; np++)
                                        for (nq = 0; nq <= n3 + n4; nq++)
                                        {
                                            w1max = np / 2;
                                            w2max = nq / 2;
                                            for (w1 = 0; w1 <= w1max; w1++)
                                                for (w2 = 0; w2 <= w2max; w2++)
                                                {
                                                    Gz = pow(-1, np) * fnp[np] * fnq[nq] * FactorialD(np)
                                                         * FactorialD(nq) * pow(gammap, w1 - np) * pow(gammaq,
                                                                 w2 - nq)
                                                         * FactorialD(np + nq - 2 * w1 - 2 * w2)
                                                         * pow(gammapq, np + nq - 2 * w1 - 2 * w2)
                                                         / (FactorialD(w1) * FactorialD(w2) * FactorialD(np - 2 * w1) * FactorialD(nq - 2 * w2));
                                                    txmax = (lp + lq - 2 * u1 - 2 * u2) / 2;
                                                    tymax = (mp + mq - 2 * v1 - 2 * v2) / 2;
                                                    tzmax = (np + nq - 2 * w1 - 2 * w2) / 2;
                                                    for (tx = 0; tx <= txmax; tx++)
                                                        for (ty = 0; ty <= tymax; ty++)
                                                            for (tz = 0; tz <= tzmax; tz++)
                                                            {
                                                                zeta = lp + lq + mp + mq + np + nq - 2 * u1 - 2
                                                                       * u2 - 2 * v1 - 2 * v2 - 2 * w1 - 2 * w2
                                                                       - tx - ty - tz;
                                                                result += Gx * Gy * Gz * F[zeta]
                                                                          * pow(-1, tx + ty + tz) * pow(
                                                                              PQx,
                                                                              lp + lq - 2
                                                                              * u1 - 2
                                                                              * u2 - 2
                                                                              * tx)
                                                                          * pow(PQy,
                                                                                mp + mq - 2 * v1 - 2 * v2 - 2 * ty)
                                                                          * pow(PQz,
                                                                                np + nq - 2 * w1 - 2 * w2 - 2 * tz)
                                                                          / (pow(
                                                                                 4,
                                                                                 u1 + u2 + tx + v1 + v2 + ty + w1
                                                                                 + w2 + tz) * pow(gammapq, tx)
                                                                             * pow(gammapq, ty) * pow(gammapq, tz)
                                                                             * FactorialD(lp + lq - 2 * u1 - 2 * u2 - 2 * tx) * FactorialD(tx) * FactorialD(mp + mq - 2 * v1 - 2 * v2 - 2 * ty) * FactorialD(ty)
                                                                             * FactorialD(np + nq - 2 * w1 - 2 * w2 - 2 * tz) * FactorialD(tz));
                                                            }
                                                }
                                        }
                                }
                        }
                }
        }

    free_array(F);
    free_array(flp);
    free_array(fmp);
    free_array(fnp);
    free_array(flq);
    free_array(fmq);
    free_array(fnq);

    return result * pfac;
}


