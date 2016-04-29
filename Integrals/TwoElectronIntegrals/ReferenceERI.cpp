#include <pulsar/system/AOOrdering.hpp>
#include <pulsar/system/SphericalTransformIntegral.hpp>

#include "../Common.hpp"
#include "ReferenceERI.hpp"

using namespace pulsar::output;
using namespace pulsar::exception;
using namespace pulsar::system;

////////////////////////////////////////////////////////////////////////////////////////
// In ValeevRef.cpp
double ValeevRef_eri(int l1, int m1, int n1, double alpha1, const double* A,
                     int l2, int m2, int n2, double alpha2, const double* B,
                     int l3, int m3, int n3, double alpha3, const double* C,
                     int l4, int m4, int n4, double alpha4, const double* D);
////////////////////////////////////////////////////////////////////////////////////////




ReferenceERI::ReferenceERI(ID_t id)
    : TwoElectronIntegral(id)
{ }



ReferenceERI::~ReferenceERI()
{ }



uint64_t ReferenceERI::Calculate_(size_t deriv,
                                  size_t shell1, size_t shell2,
                                  size_t shell3, size_t shell4,
                                  double * outbuffer, size_t bufsize)
{
    if(deriv != 0)
        throw NotYetImplementedException("Not Yet Implemented: Overlap integral with deriv != 0");

    const BasisSetShell & sh1 = bs1_->Shell(shell1);
    const BasisSetShell & sh2 = bs2_->Shell(shell2);
    const BasisSetShell & sh3 = bs3_->Shell(shell3);
    const BasisSetShell & sh4 = bs4_->Shell(shell4);

    size_t nfunc = sh1.NFunctions() * sh2.NFunctions() * sh3.NFunctions() * sh4.NFunctions();
    if(bufsize < nfunc)
        throw GeneralException("Buffer to small for ERI", "bufsize", bufsize, "nfunc", nfunc);

    // lots of loops. This isn't really meant to be fast....
    size_t idx = 0;

    for(size_t ng1 = 0; ng1 < sh1.NGeneral(); ng1++)
    for(size_t ng2 = 0; ng2 < sh2.NGeneral(); ng2++)
    for(size_t ng3 = 0; ng3 < sh3.NGeneral(); ng3++)
    for(size_t ng4 = 0; ng4 < sh4.NGeneral(); ng4++)
    {
        // now loop over all cartesian components
        const auto & cartorder1 = CartesianOrdering(sh1.GeneralAM(ng1));
        const auto & cartorder2 = CartesianOrdering(sh2.GeneralAM(ng2));
        const auto & cartorder3 = CartesianOrdering(sh3.GeneralAM(ng3));
        const auto & cartorder4 = CartesianOrdering(sh4.GeneralAM(ng4));

        for(const auto & c1 : cartorder1)
        for(const auto & c2 : cartorder2)
        for(const auto & c3 : cartorder3)
        for(const auto & c4 : cartorder4)
        {
            double myint = 0.0;

            // now the primitives
            for(size_t i = 0; i < sh1.NPrim(); i++)
            for(size_t j = 0; j < sh2.NPrim(); j++)
            for(size_t k = 0; k < sh3.NPrim(); k++)
            for(size_t l = 0; l < sh4.NPrim(); l++)
            {
                // now we can calculate the beast
                double val = ValeevRef_eri(c1[0], c1[1], c1[2], sh1.GetAlpha(i), sh1.CoordsPtr(),
                                           c2[0], c2[1], c2[2], sh2.GetAlpha(j), sh2.CoordsPtr(),  
                                           c3[0], c3[1], c3[2], sh3.GetAlpha(k), sh3.CoordsPtr(),  
                                           c4[0], c4[1], c4[2], sh4.GetAlpha(l), sh4.CoordsPtr()); 

                myint += val * sh1.GetCoef(ng1, i) * sh2.GetCoef(ng2, j) * sh3.GetCoef(ng3, k) * sh4.GetCoef(ng4, l);
            }

            outbuffer[idx++] = myint;
        }
    }

    return nfunc;
}



void ReferenceERI::SetBases_(const std::string & bs1, const std::string & bs2,
                             const std::string & bs3, const std::string & bs4)
{
    out.Debug("ReferenceERI: Initializing with bases %? %? %? %?\n", bs1, bs2, bs3, bs4);

    if(!(InitialWfn().system))
        throw GeneralException("Error - not given a system in the initial wavefunction");

    const BasisSet basisset1 = InitialWfn().system->GetBasisSet(bs1);
    const BasisSet basisset2 = InitialWfn().system->GetBasisSet(bs2);
    const BasisSet basisset3 = InitialWfn().system->GetBasisSet(bs3);
    const BasisSet basisset4 = InitialWfn().system->GetBasisSet(bs4);

    // from common components
    bs1_ = NormalizeBasis(Cache(), out, basisset1);
    bs2_ = NormalizeBasis(Cache(), out, basisset2);
    bs3_ = NormalizeBasis(Cache(), out, basisset3);
    bs4_ = NormalizeBasis(Cache(), out, basisset4);
    bs1_->Print(out);
}

