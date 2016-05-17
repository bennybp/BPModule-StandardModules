#include <pulsar/output/OutputStream.hpp>
#include <pulsar/system/AOOrdering.hpp>
#include <pulsar/system/NShellFunction.hpp>
#include <pulsar/system/SphericalTransformIntegral.hpp>
#include <pulsar/constants.h>

#include "../Common.hpp"
#include "Overlap.hpp"


// Get a value of S_IJ
#define S_IJ(i,j) (s_ij[((i)*(nam2) + j)])

using namespace pulsar::modulemanager;
using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::datastore;



uint64_t Overlap::Calculate_(uint64_t deriv,
                             uint64_t shell1, uint64_t shell2,
                             double * outbuffer, size_t bufsize)
{
    if(work_.size() == 0)
        throw GeneralException("Workspace not allocated. Did you set the bases?");
    if(deriv != 0)
        throw NotYetImplementedException("Not Yet Implemented: Overlap integral with deriv != 0");

    const BasisSetShell & sh1 = bs1_->Shell(shell1);
    const BasisSetShell & sh2 = bs2_->Shell(shell2);

    const size_t nfunc = sh1.NFunctions() * sh2.NFunctions();

    if(bufsize < nfunc)
        throw GeneralException("Buffer is too small", "size", bufsize, "required", nfunc);



    // degree of general contraction
    size_t ngen1 = sh1.NGeneral();
    size_t ngen2 = sh2.NGeneral();

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


    // coordinates from each shell
    const CoordType xyz1 = sh1.GetCoords();
    const CoordType xyz2 = sh2.GetCoords();

    const double AB[3] = { xyz1[0] - xyz2[0], xyz1[1] - xyz2[1], xyz1[2] - xyz2[2] };
    const double AB2[3] = { AB[0]*AB[0], AB[1]*AB[1], AB[2]*AB[2] };

    // Used for dimensioning and loops. Storage goes from
    // [0, am], so we need to add one.
    int nam1 = std::abs(am1) + 1;
    int nam2 = std::abs(am2) + 1;

    // are we calculating a dipole integral?
    bool isdipole = ( inttype_ == IntegralType_::Dipole_x ||
                      inttype_ == IntegralType_::Dipole_y ||
                      inttype_ == IntegralType_::Dipole_z );

    // offset due to the dipole operator (ie, are we incrementing x, y, or z)
    std::array<int, 3> dipoff{0, 0, 0};

    if(isdipole)
    {
        nam2++; // need one more due to the dipole operator

        if(inttype_ == IntegralType_::Dipole_x)
            dipoff[0]++;
        else if(inttype_ == IntegralType_::Dipole_y)
            dipoff[1]++;
        else if(inttype_ == IntegralType_::Dipole_z)
            dipoff[2]++;
    }

    // We need to zero the workspace. Actually, not all of it,
    // but this is easier
    std::fill(work_.begin(), work_.end(), 0.0);

    /////////////////////////////////////////////////////////
    // General notes about the following
    //
    // This is the OS algorithm taken pretty much verbatim
    // from Helgaker, Jorgensen, & Olsen. For each x,y,z
    // direction, we construct an array of Sij of length nam1*nam2.
    // These are placed in the pre-allocated workspace.
    //
    // We then handle multiplication with coefficients within the
    // general contractions.
    //
    // The indexing of these arrays is pretty straightforward.
    // Sij = ptr[i*nam2+j]. This is in the S_IJ macro, hopefully
    // to make the code clearer.
    /////////////////////////////////////////////////////////
    // loop over primitives
    const size_t nprim1 = sh1.NPrim();
    const size_t nprim2 = sh2.NPrim();

    for(size_t a = 0; a < nprim1; a++)
    {
        const double a1 = sh1.Alpha(a);
        const double a1xyz[3] = { a1*xyz1[0], a1*xyz1[1], a1*xyz1[2] };

        for(size_t b = 0; b < nprim2; b++)
        {
            const double a2 = sh2.Alpha(b);
            const double oop = 1.0/(a1 + a2); // = 1/p = 1/(a1 + a2)
            const double mu = a1*a2*oop; // (a1*a2)/(a1+a2)


            const double oo2p = 0.5*oop;
            const double a2xyz[3] = { a2*xyz2[0], a2*xyz2[1], a2*xyz2[2] };

            const double P[3] = { (a1xyz[0]+a2xyz[0])*oop,
                                  (a1xyz[1]+a2xyz[1])*oop,
                                  (a1xyz[2]+a2xyz[2])*oop };

            const double PA[3] = { P[0] - xyz1[0], P[1] - xyz1[1], P[2] - xyz1[2] };
            const double PB[3] = { P[0] - xyz2[0], P[1] - xyz2[1], P[2] - xyz2[2] };


            // three cartesian directions
            for(int d = 0; d < 3; d++)
            {
                // the workspace for this direction
                // THIS IS THEN ACCESSED THROUGH THE S_IJ MACRO
                double * const RESTRICT s_ij = xyzwork_[d];

                S_IJ(0,0) = sqrt(PI * oop) * exp(-mu * AB2[d]);

                // do j = 0 for all remaining i
                for(int i = 1; i < nam1; i++)
                {
                    S_IJ(i,0) = PA[d]*S_IJ(i-1,0);
                    if(i > 1)
                        S_IJ(i,0) += (i-1)*oo2p*S_IJ(i-2,0);
                }

 
                // now do i = 0 for all remaining j
                for(int j = 1; j < nam2; j++)
                {
                    S_IJ(0,j) = PB[d]*S_IJ(0,j-1);
                    if(j > 1)
                        S_IJ(0,j) += (j-1)*oo2p*S_IJ(0,j-2);
                }

                // now all the rest
                for(int i = 1; i < nam1; i++)
                for(int j = 1; j < nam2; j++)
                {
                    S_IJ(i,j) = PB[d]*S_IJ(i,j-1) + oo2p*i*S_IJ(i-1,j-1);
                    if(j > 1)
                        S_IJ(i,j) += oo2p*(j-1)*S_IJ(i,j-2);
                }
            }

            // general contraction and combined am
            size_t outidx = 0;
            for(size_t g1 = 0; g1 < ngen1; g1++)
            for(size_t g2 = 0; g2 < ngen2; g2++)
            {
                // go over the orderings for this AM
                for(const IJK & ijk1 : *(sh1_ordering[g1]))
                for(const IJK & ijk2 : *(sh2_ordering[g2]))
                {
                    const int xidx = ijk1[0]*nam2 + ( ijk2[0] + dipoff[0] );
                    const int yidx = ijk1[1]*nam2 + ( ijk2[1] + dipoff[1] );
                    const int zidx = ijk1[2]*nam2 + ( ijk2[2] + dipoff[2] );

                    const double val = xyzwork_[0][xidx] *
                                       xyzwork_[1][yidx] *
                                       xyzwork_[2][zidx];

                    // remember: a and b are indices of primitives
                    sourcework_[outidx++] += val * sh1.Coef(g1, a) * sh2.Coef(g2, b);
                }
            }
        }
    }

    // performs the spherical transform, if necessary
    CartesianToSpherical_2Center(sh1, sh2, sourcework_, outbuffer, transformwork_);

    return nfunc;
}



void Overlap::SetBases_(const System & sys,
                        const std::string & bs1, const std::string & bs2)
{
    // determine the integral we are calculating
    std::string inttype_str = Options().Get<std::string>("TYPE");
    if(inttype_str == "OVERLAP")
        inttype_ = IntegralType_::Overlap;
    else if(inttype_str == "DIPOLE_X")
        inttype_ = IntegralType_::Dipole_x;
    else if(inttype_str == "DIPOLE_Y")
        inttype_ = IntegralType_::Dipole_y;
    else if(inttype_str == "DIPOLE_Z")
        inttype_ = IntegralType_::Dipole_z;
    else
        throw GeneralException("Unknown integral type", "type", inttype_str);

    // from common components
    bs1_ = NormalizeBasis(Cache(), out, sys.GetBasisSet(bs1));
    bs2_ = NormalizeBasis(Cache(), out, sys.GetBasisSet(bs2));

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
