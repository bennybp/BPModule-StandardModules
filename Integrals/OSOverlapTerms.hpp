#pragma once

namespace psr_modules {
namespace integrals {
namespace detail {


/*! \brief Calculation of Obara-Saika overlap terms
 *
 * For each cartesian direction \f$ d \in x,y,z\f$ and two angular momentum
 * components \f$ i,j \f$,
 * this funtion calculates \f$ \left< G_1 | G_2 \right> \f$, where
 *
 * \f[
 *    G_1 = d^i e^{-\alpha_1 r^2}
 * \f]
 *
 * and similarly for \f$G_2\f$.
 *
 * All pairs of i, j are calculated, with i in the range [0, \p nam1)
 * and j in the range [0, \p nam2).
 *
 * The results are stored as \p outbuffer[d][i*nam2+j]
 *
 * \param [in] alpha1 The exponent parameter of the first gaussian
 * \param [in] xyz1   The coordinates of the first gaussian
 * \param [in] alpha2 The exponent parameter of the second gaussian
 * \param [in] xyz2   The coordinates of the second gaussian
 * \param [in] nam1   Number of angular momentum terms to calculate
 *                    (for each direction) for the first center
 * \param [in] nam2   Number of angular momentum terms to calculate
 *                    (for each direction) for the second center
 * \param [out] outbuffer Results buffer. Dimentions should be [3][nam1*nam2] 
 */
void os_overlap_terms(const double alpha1, const double xyz1[3],
                      const double alpha2, const double xyz2[3],
                      int nam1, int nam2,
                      double ** outbuffer);


} // close namespace detail
} // close namespace integrals
} // close namespace psr_modules
