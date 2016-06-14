#include <pulsar/system/AOOrdering.hpp>
#include <pulsar/system/SphericalTransformIntegral.hpp>

#include "../Common.hpp"
#include "OSOverlap.hpp"
#include "Overlap.hpp"


// Get a value of S_IJ
#define S_IJ(i,j) (s_ij[((i)*(nam2) + j)])

using namespace pulsar::modulemanager;
using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::datastore;



uint64_t Overlap::Calculate_(uint64_t shell1, uint64_t shell2,
                             double * outbuffer, size_t bufsize)
{
    if(work_.size() == 0)
        throw GeneralException("Workspace not allocated. Did you set the bases?");

    const BasisSetShell & sh1 = bs1_->Shell(shell1);
    const BasisSetShell & sh2 = bs2_->Shell(shell2);

    const size_t nfunc = sh1.NFunctions() * sh2.NFunctions();

    if(bufsize < nfunc)
        throw GeneralException("Buffer is too small", "size", bufsize, "required", nfunc);



    // degree of general contraction
    const size_t ngen1 = sh1.NGeneral();
    const size_t ngen2 = sh2.NGeneral();

    ////////////////////////////////////////
    // These vectors store the ordering of
    // cartesian basis functions for each
    // general contraction. Remember that
    // general contractions include combined AM
    // contractions.
    //////////////////////////////////////////
    std::vector<const std::vector<IJK> *> sh1_ordering, sh2_ordering;

    for(size_t g1 = 0; g1 < ngen1; g1++)
        sh1_ordering.push_back(&CartesianOrdering(sh1.GeneralAM(g1)));
    for(size_t g2 = 0; g2 < ngen2; g2++)
        sh2_ordering.push_back(&CartesianOrdering(sh2.GeneralAM(g2)));


    // The total AM of the shell. May be negative
    const int am1 = sh1.AM();
    const int am2 = sh2.AM();

    // Used for dimensioning and loops. Storage goes from
    // [0, am], so we need to add one.
    const int nam1 = std::abs(am1) + 1;
    const int nam2 = std::abs(am2) + 1;

    // coordinates from each shell
    const double * xyz1 = sh1.CoordsPtr();
    const double * xyz2 = sh2.CoordsPtr();

    // We need to zero the workspace. Actually, not all of it,
    // but this is easier
    std::fill(work_.begin(), work_.end(), 0.0);

    // loop over primitives
    const size_t nprim1 = sh1.NPrim();
    const size_t nprim2 = sh2.NPrim();

    for(size_t a = 0; a < nprim1; a++)
    for(size_t b = 0; b < nprim2; b++)
    {
        OSOverlap(sh1.Alpha(a), xyz1,
                  sh2.Alpha(b), xyz2,
                  nam1, nam2, xyzwork_);

        // general contraction and combined am
        size_t outidx = 0;
        for(size_t g1 = 0; g1 < ngen1; g1++)
        for(size_t g2 = 0; g2 < ngen2; g2++)
        {
            const double prefac = sh1.Coef(g1, a) * sh2.Coef(g2, b);

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



void Overlap::Initialize_(unsigned int deriv,
                          const Wavefunction & wfn,
                          const BasisSet & bs1,
                          const BasisSet & bs2)
{
    if(deriv != 0)
        throw NotYetImplementedException("Not Yet Implemented: Overlap integral with deriv != 0");

    // from common components
    bs1_ = NormalizeBasis(Cache(), out, bs1);
    bs2_ = NormalizeBasis(Cache(), out, bs2);

    ///////////////////////////////////////
    // Determine the size of the workspace
    ///////////////////////////////////////

    // storage size for each x,y,z component
    int max1 = bs1_->MaxAM();
    int max2 = bs2_->MaxAM();
    size_t worksize = (max1+1)*(max2+1);  // for each component, we store [0, am]

    // find the maximum number of cartesian functions, not including general contraction
    size_t maxsize1 = bs1_->MaxProperty(NCartesianGaussianForShellAM);
    size_t maxsize2 = bs2_->MaxProperty(NCartesianGaussianForShellAM);
    size_t transformwork_size = maxsize1 * maxsize2; 

    // find the maximum number of cartesian functions, including general contraction
    maxsize1 = bs1_->MaxProperty(NCartesianGaussianInShell);
    maxsize2 = bs2_->MaxProperty(NCartesianGaussianInShell);
    size_t sourcework_size = maxsize1 * maxsize2;

    // allocate all at once, then partition
    work_.resize(3*worksize + transformwork_size + sourcework_size);
    xyzwork_[0] = work_.data();
    xyzwork_[1] = xyzwork_[0] + worksize;
    xyzwork_[2] = xyzwork_[1] + worksize;
    transformwork_ = xyzwork_[2] + worksize;
    sourcework_ = transformwork_ + transformwork_size;
}
