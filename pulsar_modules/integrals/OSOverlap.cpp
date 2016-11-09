#include <pulsar/system/AOOrdering.hpp>
#include <pulsar/system/SphericalTransformIntegral.hpp>

#include "Common/BasisSetCommon.hpp"
#include "Integrals/OSOverlapTerms.hpp"
#include "Integrals/OSOverlap.hpp"


using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::datastore;

namespace psr_modules {
namespace integrals {


uint64_t OSOverlap::calculate_(uint64_t shell1, uint64_t shell2,
                               double * outbuffer, size_t bufsize)
{
    const BasisSetShell & sh1 = bs1_->shell(shell1);
    const BasisSetShell & sh2 = bs2_->shell(shell2);

    const size_t nfunc = sh1.n_functions() * sh2.n_functions();

    if(bufsize < nfunc)
        throw GeneralException("Buffer is too small", "size", bufsize, "required", nfunc);


    // degree of general contraction
    const size_t ngen1 = sh1.n_general_contractions();
    const size_t ngen2 = sh2.n_general_contractions();

    ////////////////////////////////////////
    // These vectors store the ordering of
    // cartesian basis functions for each
    // general contraction. Remember that
    // general contractions include combined AM
    // contractions.
    //////////////////////////////////////////
    std::vector<const std::vector<IJK> *> sh1_ordering, sh2_ordering;

    for(size_t g1 = 0; g1 < ngen1; g1++)
        sh1_ordering.push_back(&cartesian_ordering(sh1.general_am(g1)));
    for(size_t g2 = 0; g2 < ngen2; g2++)
        sh2_ordering.push_back(&cartesian_ordering(sh2.general_am(g2)));


    // The total AM of the shell. May be negative
    const int am1 = sh1.am();
    const int am2 = sh2.am();

    // Used for dimensioning and loops. Storage goes from
    // [0, am], so we need to add one.
    const int nam1 = std::abs(am1) + 1;
    const int nam2 = std::abs(am2) + 1;

    // coordinates from each shell
    const double * xyz1 = sh1.coords_ptr();
    const double * xyz2 = sh2.coords_ptr();

    // We need to zero the workspace. Actually, not all of it,
    // but this is easier
    std::fill(work_.begin(), work_.end(), 0.0);

    // loop over primitives
    const size_t nprim1 = sh1.n_primitives();
    const size_t nprim2 = sh2.n_primitives();

    for(size_t a = 0; a < nprim1; a++)
    for(size_t b = 0; b < nprim2; b++)
    {
        // Calculate all the S_ij terms
        detail::os_overlap_terms(sh1.alpha(a), xyz1,
                                 sh2.alpha(b), xyz2,
                                 nam1, nam2, xyzwork_);

        // general contraction and combined am
        size_t outidx = 0;
        for(size_t g1 = 0; g1 < ngen1; g1++)
        for(size_t g2 = 0; g2 < ngen2; g2++)
        {
            const double prefac = sh1.coef(g1, a) * sh2.coef(g2, b);

            // go over the orderings for this AM
            for(const IJK & ijk1 : *(sh1_ordering[g1]))
            for(const IJK & ijk2 : *(sh2_ordering[g2]))
            {
                const int xidx = ijk1[0]*nam2 + ijk2[0];
                const int yidx = ijk1[1]*nam2 + ijk2[1];
                const int zidx = ijk1[2]*nam2 + ijk2[2];

                const double val = xyzwork_[0][xidx] *
                                   xyzwork_[1][yidx] *
                                   xyzwork_[2][zidx];

                // remember: a and b are indices of primitives
                sourcework_[outidx++] += prefac * val;
            }
        }
    }

    // performs the spherical transform, if necessary
    CartesianToSpherical_2Center(sh1, sh2, sourcework_, outbuffer, transformwork_, 1);

    return nfunc;
}



void OSOverlap::initialize_(unsigned int deriv,
                           const Wavefunction & wfn,
                           const BasisSet & bs1,
                           const BasisSet & bs2)
{
    if(deriv != 0)
        throw NotYetImplementedException("Not Yet Implemented: Overlap integral with deriv != 0");

    // from common components
    bs1_ = NormalizeBasis(cache(), out, bs1);
    bs2_ = NormalizeBasis(cache(), out, bs2);

    ///////////////////////////////////////
    // Determine the size of the workspace
    ///////////////////////////////////////

    // storage size for each x,y,z component
    int max1 = bs1_->max_am();
    int max2 = bs2_->max_am();
    size_t worksize = (max1+1)*(max2+1);  // for each component, we store [0, am]

    // find the maximum number of cartesian functions, not including general contraction
    size_t maxsize1 = bs1_->max_property(n_cartesian_gaussian_for_shell_am);
    size_t maxsize2 = bs2_->max_property(n_cartesian_gaussian_for_shell_am);
    size_t transformwork_size = maxsize1 * maxsize2; 

    // find the maximum number of cartesian functions, including general contraction
    maxsize1 = bs1_->max_property(n_cartesian_gaussian_in_shell);
    maxsize2 = bs2_->max_property(n_cartesian_gaussian_in_shell);
    size_t sourcework_size = maxsize1 * maxsize2;

    // allocate all at once, then partition
    work_.resize(3*worksize + transformwork_size + sourcework_size);
    xyzwork_[0] = work_.data();
    xyzwork_[1] = xyzwork_[0] + worksize;
    xyzwork_[2] = xyzwork_[1] + worksize;
    transformwork_ = xyzwork_[2] + worksize;
    sourcework_ = transformwork_ + transformwork_size;
}


} // close namespace integrals
} // close namespace psr_modules
