#include <pulsar/output/OutputStream.hpp>
#include <pulsar/system/AOOrdering.hpp>
#include <pulsar/system/SphericalTransformIntegral.hpp>

#include "Common/BasisSetCommon.hpp"
#include "Integrals/ReferenceERI.hpp"

using namespace pulsar::output;
using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::datastore;

////////////////////////////////////////////////////////////////////////////////////////
// In ValeevRef.cpp
double ValeevRef_eri(int l1, int m1, int n1, double alpha1, const double* A,
                     int l2, int m2, int n2, double alpha2, const double* B,
                     int l3, int m3, int n3, double alpha3, const double* C,
                     int l4, int m4, int n4, double alpha4, const double* D);
////////////////////////////////////////////////////////////////////////////////////////




uint64_t ReferenceERI::calculate_(size_t shell1, size_t shell2,
                                  size_t shell3, size_t shell4,
                                  double * outbuffer, size_t bufsize)
{
    const BasisSetShell & sh1 = bs1_->shell(shell1);
    const BasisSetShell & sh2 = bs2_->shell(shell2);
    const BasisSetShell & sh3 = bs3_->shell(shell3);
    const BasisSetShell & sh4 = bs4_->shell(shell4);

    size_t nfunc = sh1.n_functions() * sh2.n_functions() * sh3.n_functions() * sh4.n_functions();

    if(bufsize < nfunc)
        throw GeneralException("Buffer to small for ERI", "bufsize", bufsize, "nfunc", nfunc);

    // lots of loops. This isn't really meant to be fast....
    size_t idx = 0;

    for(size_t ng1 = 0; ng1 < sh1.n_general_contractions(); ng1++)
    {
        const auto & cartorder1 = cartesian_ordering(sh1.general_am(ng1));

        for(size_t ng2 = 0; ng2 < sh2.n_general_contractions(); ng2++)
        {
            const auto & cartorder2 = cartesian_ordering(sh2.general_am(ng2));

            for(size_t ng3 = 0; ng3 < sh3.n_general_contractions(); ng3++)
            {
                const auto & cartorder3 = cartesian_ordering(sh3.general_am(ng3));

                for(size_t ng4 = 0; ng4 < sh4.n_general_contractions(); ng4++)
                {
                    const auto & cartorder4 = cartesian_ordering(sh4.general_am(ng4));

                    // now loop over all cartesian components
                    for(const auto & c1 : cartorder1)
                    for(const auto & c2 : cartorder2)
                    for(const auto & c3 : cartorder3)
                    for(const auto & c4 : cartorder4)
                    {
                        double myint = 0.0;

                        // now the primitives
                        for(size_t i = 0; i < sh1.n_primitives(); i++)
                        for(size_t j = 0; j < sh2.n_primitives(); j++)
                        for(size_t k = 0; k < sh3.n_primitives(); k++)
                        for(size_t l = 0; l < sh4.n_primitives(); l++)
                        {
                            // now we can calculate the beast
                            double val = ValeevRef_eri(c1[0], c1[1], c1[2], sh1.get_alpha(i), sh1.coords_ptr(),
                                                       c2[0], c2[1], c2[2], sh2.get_alpha(j), sh2.coords_ptr(),  
                                                       c3[0], c3[1], c3[2], sh3.get_alpha(k), sh3.coords_ptr(),  
                                                       c4[0], c4[1], c4[2], sh4.get_alpha(l), sh4.coords_ptr()); 

                            myint += val * sh1.get_coef(ng1, i) * sh2.get_coef(ng2, j) * sh3.get_coef(ng3, k) * sh4.get_coef(ng4, l);
                        }

                        sourcework_[idx++] = myint;
                    }
                }
            }
        }
    }


    CartesianToSpherical_4Center(sh1, sh2, sh3, sh4, sourcework_, outbuffer, transformwork_, 1);

    return nfunc;
}



void ReferenceERI::initialize_(unsigned int deriv,
                               const Wavefunction & wfn,
                               const BasisSet & bs1,
                               const BasisSet & bs2,
                               const BasisSet & bs3,
                               const BasisSet & bs4)
{
    if(deriv != 0)
        throw NotYetImplementedException("Not Yet Implemented: ReferenceERI integral with deriv != 0");

    // from common components
    bs1_ = NormalizeBasis(cache(), out, bs1);
    bs2_ = NormalizeBasis(cache(), out, bs2);
    bs3_ = NormalizeBasis(cache(), out, bs3);
    bs4_ = NormalizeBasis(cache(), out, bs4);

    size_t maxsize1 = bs1_->max_property(n_cartesian_gaussian_for_shell_am);
    size_t maxsize2 = bs2_->max_property(n_cartesian_gaussian_for_shell_am);
    size_t maxsize3 = bs3_->max_property(n_cartesian_gaussian_for_shell_am);
    size_t maxsize4 = bs4_->max_property(n_cartesian_gaussian_for_shell_am);
    size_t transformwork_size = maxsize1*maxsize2*maxsize3*maxsize4;
    
    maxsize1 =  bs1_->max_property(n_cartesian_gaussian_in_shell);
    maxsize2 =  bs2_->max_property(n_cartesian_gaussian_in_shell);
    maxsize3 =  bs3_->max_property(n_cartesian_gaussian_in_shell);
    maxsize4 =  bs4_->max_property(n_cartesian_gaussian_in_shell);
    size_t sourcework_size = maxsize1*maxsize2*maxsize3*maxsize4;

    work_.resize(sourcework_size+transformwork_size);
    sourcework_ = work_.data();
    transformwork_ = sourcework_ + sourcework_size;   
}

