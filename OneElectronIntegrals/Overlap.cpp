#include <pulsar/output/Output.hpp>
#include <pulsar/modulemanager/ModuleManager.hpp>
#include <pulsar/system/BasisSet.hpp>
#include <pulsar/system/AOOrdering.hpp>
#include <pulsar/system/NShellFunction.hpp>
#include <pulsar/math/Factorial.hpp>
#include <pulsar/math/MinMax.hpp>
#include <pulsar/system/SphericalTransformIntegral.hpp>
#include <pulsar/constants.h>

#include "Overlap.hpp"


using namespace pulsar::modulemanager;
using namespace pulsar::exception;
using namespace pulsar::system;


// coordinates not used
static BasisShellInfo NormalizeShell_(const BasisShellInfo & shell, const CoordType &)
{
    static const double norm_fac[25] =
    {
    /* l =    0 */  5.56832799683170785             ,
    /* l =    1 */  2.78416399841585392             ,
    /* l =    2 */  4.17624599762378088             ,
    /* l =    3 */  10.4406149940594522             ,
    /* l =    4 */  36.5421524792080827             ,
    /* l =    5 */  164.439686156436372             ,
    /* l =    6 */  904.418273860400048             ,
    /* l =    7 */  5878.71878009260031             ,
    /* l =    8 */  44090.3908506945023             ,
    /* l =    9 */  374768.32223090327              ,
    /* l =   10 */  3560299.06119358106             ,
    /* l =   11 */  37383140.1425326012             ,
    /* l =   12 */  429906111.639124913             ,
    /* l =   13 */  5373826395.48906142             ,
    /* l =   14 */  72546656339.1023291             ,
    /* l =   15 */  1051926516916.98377             ,
    /* l =   16 */  16304861012213.2485             ,
    /* l =   17 */  269030206701518.6               ,
    /* l =   18 */  4708028617276575.5              ,
    /* l =   19 */  87098529419616646.7             ,
    /* l =   20 */  1.69842132368252461e+18         ,
    /* l =   21 */  3.48176371354917545e+19         ,
    /* l =   22 */  7.48579198413072722e+20         ,
    /* l =   23 */  1.68430319642941362e+22         ,
    /* l =   24 */  3.95811251160912202e+23         ,
    };



    // common to all general contractions
    const size_t nprim = shell.NPrim();
    const double * const alpha = shell.AlphaPtr();

    BasisShellInfo newshell(shell);

    for(size_t n = 0; n < shell.NGeneral(); n++)
    {
        const int iam = shell.GeneralAM(n);
        const double am = static_cast<double>(iam);
        const double m = am + 1.5;
        const double m2 = 0.5 * m;

        std::vector<double> coefs = shell.GetCoefs(n);

        double sum = 0.0;

        for(size_t i = 0; i < nprim; ++i)
        {
            const double a1 = alpha[i];
            const double c1 = coefs[i];

            for(size_t j = 0; j < nprim; j++)
            {
                const double a2 = alpha[j];
                const double c2 = coefs[j];
                sum += ( c1 * c2 *  pow(a1*a2, m2) ) / ( pow(a1+a2, m) );
            }
        }

        const double norm = 1.0 / sqrt(sum * norm_fac[iam]);

        // apply the rest of the normalization and store
        for (size_t i = 0; i < nprim; ++i)
            newshell.SetCoef(n, i, coefs[i] * norm * pow(alpha[i], m2));
    }

    return newshell;
}


std::shared_ptr<BasisSet> Overlap::Normalize_(const BasisSet & bs)
{
    std::string cachekey = std::string("bs:") + bs.MyHash().String();
    if(Cache().Count(cachekey)) // options are unimportant
    {
        auto ret = Cache().Get<std::shared_ptr<BasisSet>>(cachekey);
        out.Debug("Overlap: Found normalized basis in cache: %? -> %?\n",
                  bs.MyHash().String(), ret->MyHash().String());
        return ret;
    }

    bs.Print(out);
    std::shared_ptr<BasisSet> newbs = std::make_shared<BasisSet>(bs.Transform(NormalizeShell_));

    // add to the cache
    Cache().Set(cachekey, newbs);

    return newbs;
}




Overlap::Overlap(ID_t id)
    : OneElectronIntegral(id)
{ }

Overlap::~Overlap()
{ }



uint64_t Overlap::Calculate_(uint64_t deriv,
                             uint64_t shell1, uint64_t shell2,
                             double * outbuffer, size_t bufsize)
{
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
    const size_t ncart1 = NCartesianGaussian(sh1);
    const size_t ncart2 = NCartesianGaussian(sh2);
    const size_t ncart12 = ncart1*ncart2;

    double * const RESTRICT srcptr = srcwork_.get();
    std::fill(srcptr, srcptr + ncart12, 0.0);

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
                double * const RESTRICT ptr = work3_[d];
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

                    const double val = work3_[0][xidx] *
                                       work3_[1][yidx] *
                                       work3_[2][zidx];

                    srcptr[outidx++] += val * sh1.Coef(g1, i) * sh2.Coef(g2, j);
                }
            }
        }
    }

    CartesianToSpherical_OneElectron(sh1, sh2, srcptr,
                                     outbuffer, transformwork_.get());

    return nfunc;
}



void Overlap::SetBases_(const std::string & bs1, const std::string & bs2)
{
    out.Debug("Overlap: Initializing with bases %? %?\n", bs1, bs2);

    //! \todo - check if bases are set

    if(!(InitialWfn().system))
        throw GeneralException("Error - not given a system in the initial wavefunction");

    const BasisSet basisset1 = InitialWfn().system->GetBasisSet(bs1);
    const BasisSet basisset2 = InitialWfn().system->GetBasisSet(bs2);

    bs1_ = Normalize_(basisset1);
    bs2_ = Normalize_(basisset2);

    bs1_->Print(out);

    ///////////////////////////////////////
    // Determine the size of the workspace
    ///////////////////////////////////////
    int max1 = bs1_->MaxAM();
    int max2 = bs2_->MaxAM();
    size_t worksize = (max1+1)*(max2+1);

    work_ = std::unique_ptr<double[]>(new double[3*worksize]);
    work3_[0] = work_.get();
    work3_[1] = work3_[0] + worksize;
    work3_[2] = work3_[1] + worksize;

    size_t transformwork_size = 10*NCartesianGaussian(max1)*NCartesianGaussian(max2);
    transformwork_ = std::unique_ptr<double[]>(new double[transformwork_size]);
    srcwork_ = std::unique_ptr<double[]>(new double[transformwork_size]);
}
