#include <pulsar/output/Output.hpp>
#include <pulsar/modulemanager/ModuleManager.hpp>
#include <pulsar/system/BasisSet.hpp>
#include <pulsar/system/AOOrdering.hpp>
#include <pulsar/system/NShellFunction.hpp>
#include <pulsar/system/SphericalTransformIntegral.hpp>
#include <pulsar/constants.h>

#include "../Common.hpp"
#include "Overlap.hpp"


using namespace pulsar::modulemanager;
using namespace pulsar::exception;
using namespace pulsar::system;


Overlap::Overlap(ID_t id)
    : OneElectronIntegral(id)
{ }

Overlap::~Overlap()
{ }


uint64_t Overlap::Calculate_(uint64_t deriv,
                             uint64_t shell1, uint64_t shell2,
                             double * outbuffer, size_t bufsize)
{
    if(!work_)
        throw GeneralException("Workspace not allocated. Did you set the bases?");
    if(deriv != 0)
        throw NotYetImplementedException("Not Yet Implemented: Overlap integral with deriv != 0");

    const BasisSetShell & sh1 = bs1_->Shell(shell1);
    const BasisSetShell & sh2 = bs2_->Shell(shell2);

    const size_t nfunc = sh1.NFunctions() * sh2.NFunctions();

    if(bufsize < nfunc)
        throw GeneralException("Buffer is to small", "size", bufsize, "required", nfunc);

    const int am1 = sh1.AM();
    const int am2 = sh2.AM();

    // Used for dimensioning and loops. Storage goes from
    // [0, am], so we need to add one.
    const int nam1 = std::abs(am1) + 1;
    const int nam2 = std::abs(am2) + 1;

    // coordinates from each shell
    const double * xyz1 = sh1.CoordsPtr();
    const double * xyz2 = sh2.CoordsPtr();

    // degree of general contraction
    size_t ngen1 = sh1.NGeneral();
    size_t ngen2 = sh2.NGeneral();

    // Number of cartesian gaussians in each shell
    const size_t ncart1 = NCartesianGaussianInShell(sh1);
    const size_t ncart2 = NCartesianGaussianInShell(sh2);
    const size_t ncart12 = ncart1*ncart2;

    std::fill(sourcework_, sourcework_ + ncart12, 0.0);

    const double AB[3] = { xyz1[0] - xyz2[0], xyz1[1] - xyz2[1], xyz1[2] - xyz2[2] };
    const double AB2 = AB[0]*AB[0] + AB[1]*AB[1] + AB[2]*AB[2];


    // loop over primitives
    const size_t nprim1 = sh1.NPrim();
    const size_t nprim2 = sh2.NPrim();

    for(size_t i = 0; i < nprim1; i++)
    {
        const double a1 = sh1.Alpha(i);
        const double a1xyz[3] = { a1*xyz1[0], a1*xyz1[1], a1*xyz1[2] };

        for(size_t j = 0; j < nprim2; j++)
        {
            const double a2 = sh2.Alpha(j);
            const double a2xyz[3] = { a2*xyz2[0], a2*xyz2[1], a2*xyz2[2] };

            const double oop = 1.0/(a1 + a2); // = 1/p = 1/(a1 + a2)
            const double oo2p = 0.5*oop;

            const double mu = a1*a2*oop; // (a1+a2)/(a1*a2)

            const double P[3] = { (a1xyz[0]+a2xyz[0])*oop,
                                  (a1xyz[1]+a2xyz[1])*oop,
                                  (a1xyz[2]+a2xyz[2])*oop };

            const double PA[3] = { P[0] - xyz1[0], P[1] - xyz1[1], P[2] - xyz1[2] };
            const double PB[3] = { P[0] - xyz2[0], P[1] - xyz2[1], P[2] - xyz2[2] };

            const double S00 = sqrt(PI * oop) * exp(-mu * AB2);

            // three cartesian directions
            for(int d = 0; d < 3; d++)
            {
                // the workspace for this direction
                double * const RESTRICT ptr = xyzwork_[d];
                ptr[0] = S00;

                for(int b = 1; b < nam1; b++)
                {
                    size_t dest = b*nam2;

                    // form b0 via bra recurrence
                    ptr[dest] = PA[d]*ptr[(b-1)*nam2];
                    if(b > 2)
                        ptr[dest] += oo2p*(b-1)*ptr[(b-2)*nam2];

                    for(int k = 1; k <= b; k++)
                    {
                        ptr[dest+k] = PB[d]*ptr[dest + k - 1];
                        ptr[dest+k] += oo2p*b*ptr[dest - nam2 + k - 1];  // (b-1)*nam2 + k - 1

                        if(k > 1)
                            ptr[dest+k] += oo2p*(k-1)*ptr[dest + k - 2];
                    }
                }
            }

            // general contraction and combined am
            size_t outidx = 0;
            for(size_t g1 = 0; g1 < ngen1; g1++)
            for(size_t g2 = 0; g2 < ngen2; g2++)
            {
                const int gam1 = sh1.GeneralAM(g1);
                const int gam2 = sh2.GeneralAM(g2);
                const auto & ijk1vec = CartesianOrdering(gam1);
                const auto & ijk2vec = CartesianOrdering(gam2);

                for(const auto & ijk1 : ijk1vec)
                for(const auto & ijk2 : ijk2vec)
                {
                    const int xidx = ijk1[0]*nam2 + ijk2[0];
                    const int yidx = ijk1[1]*nam2 + ijk2[1];
                    const int zidx = ijk1[2]*nam2 + ijk2[2];

                    const double val = xyzwork_[0][xidx] *
                                       xyzwork_[1][yidx] *
                                       xyzwork_[2][zidx];

                    sourcework_[outidx++] += val * sh1.Coef(g1, i) * sh2.Coef(g2, j);
                }
            }
        }
    }

    CartesianToSpherical_OneElectron(sh1, sh2, sourcework_,
                                     outbuffer, transformwork_);

    return nfunc;
}



void Overlap::SetBases_(const std::string & bs1, const std::string & bs2)
{
    out.Debug("Overlap: Initializing with bases %? %?\n", bs1, bs2);

    if(!(InitialWfn().system))
        throw GeneralException("Error - not given a system in the initial wavefunction");

    const BasisSet basisset1 = InitialWfn().system->GetBasisSet(bs1);
    const BasisSet basisset2 = InitialWfn().system->GetBasisSet(bs2);

    // from common components
    bs1_ = NormalizeBasis(Cache(), out, basisset1);
    bs2_ = NormalizeBasis(Cache(), out, basisset2);

    ///////////////////////////////////////
    // Determine the size of the workspace
    ///////////////////////////////////////

    // storage size for each x,y,z component
    int max1 = bs1_->MaxAM();
    int max2 = bs2_->MaxAM();
    size_t worksize = (max1+1)*(max2+1);  // for wach component, we store [0, am]

    // find the maximum number of cartesian functions, not including general contraction
    size_t maxsize1 = bs1_->MaxProperty(NCartesianGaussianForShellAM);
    size_t maxsize2 = bs2_->MaxProperty(NCartesianGaussianForShellAM);
    size_t transformwork_size = maxsize1 * maxsize2; 

    // find the maximum number of cartesian functions, including general contraction
    maxsize1 = bs1_->MaxProperty(NCartesianGaussianInShell);
    maxsize2 = bs2_->MaxProperty(NCartesianGaussianInShell);
    size_t sourcework_size = maxsize1 * maxsize2;

    // allocate all at once, then partition
    work_ = std::unique_ptr<double[]>(new double[3*worksize + transformwork_size + sourcework_size]);
    xyzwork_[0] = work_.get();
    xyzwork_[1] = xyzwork_[0] + worksize;
    xyzwork_[2] = xyzwork_[1] + worksize;
    transformwork_ = xyzwork_[2] + worksize;
    sourcework_ = transformwork_ + transformwork_size;
}
