#include <pulsar/system/AOOrdering.hpp>
#include <pulsar/system/SphericalTransformIntegral.hpp>
#include <pulsar/math/Factorial.hpp>
#include <pulsar/constants.h>

#include "Common/BasisSetCommon.hpp"
#include "Integrals/boys/Boys.hpp"
#include "Integrals/OSOneElectronPotential.hpp"
#include "Integrals/OSOneElectronPotential_LUT.hpp"


using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::datastore;
using namespace pulsar::math;


namespace psr_modules {
namespace integrals {


uint64_t OSOneElectronPotential::calculate_with_grid_(uint64_t shell1, uint64_t shell2,
                                                      const Grid & grid,
                                                      double * outbuffer, size_t bufsize)
{
    if(!sys_)
        throw PulsarException("No system given");

    const BasisSetShell & sh1 = bs1_->shell(shell1);
    const BasisSetShell & sh2 = bs2_->shell(shell2);

    const size_t nfunc = sh1.n_functions() * sh2.n_functions();

    if(bufsize < nfunc)
        throw PulsarException("Buffer is too small", "size", bufsize, "required", nfunc);

    // degree of general contraction
    size_t ngen1 = sh1.n_general_contractions();
    size_t ngen2 = sh2.n_general_contractions();

    // The total AM of the shell. May be negative
    const int am1 = sh1.am();
    const int am2 = sh2.am();

    // used for loops
    const int absam1 = std::abs(am1);
    const int absam2 = std::abs(am2);
    const int absam12 = absam1 + absam2;

    // coordinates
    const CoordType xyz1 = sh1.get_coords();
    const CoordType xyz2 = sh2.get_coords();

    // number of primitives
    const size_t nprim1 = sh1.n_primitives();
    const size_t nprim2 = sh2.n_primitives();

    const double AB[3] = { xyz1[0] - xyz2[0], xyz1[1] - xyz2[1], xyz1[2] - xyz2[2] };
    const double AB2 = AB[0]*AB[0] + AB[1]*AB[1] + AB[2]*AB[2];
    
    std::fill(work_.begin(), work_.end(), 0.0);

    for(const auto & gridpt : grid)
    {
        for(size_t a = 0; a < nprim1; a++)
        {
            const double a1 = sh1.alpha(a);
            const double a1xyz[3] = { a1*xyz1[0], a1*xyz1[1], a1*xyz1[2] };

            for(size_t b = 0; b < nprim2; b++)
            {
                const double a2 = sh2.alpha(b);
                const double p = a1 + a2;
                const double oop = 1.0/(a1 + a2); // = 1/p = 1/(a1 + a2)
                const double mu = a1*a2*oop; // (a1+a2)/(a1*a2)

                const double oo2p = 0.5*oop;
                const double a2xyz[3] = { a2*xyz2[0], a2*xyz2[1], a2*xyz2[2] };

                const double P[3] = { (a1xyz[0]+a2xyz[0])*oop,
                                      (a1xyz[1]+a2xyz[1])*oop,
                                      (a1xyz[2]+a2xyz[2])*oop };

                const double PA[3] = { P[0] - xyz1[0], P[1] - xyz1[1], P[2] - xyz1[2] };
                const double PB[3] = { P[0] - xyz2[0], P[1] - xyz2[1], P[2] - xyz2[2] };
                const double PC[3] = { P[0] - gridpt.coords[0], P[1] - gridpt.coords[1], P[2] - gridpt.coords[2] };
                const double PC2 = PC[0]*PC[0] + PC[1]*PC[1] + PC[2]*PC[2];

                // boys function
                const double T = PC2 * p;
                detail::calculate_f(amwork_[0][0], absam12, T); 
                for(int i = 0; i <= absam12; i++)
                    amwork_[0][0][i] *= 2*PI*oop*exp(-mu*AB2);

                // nested recurrence
                // we skip (0,0) since that is the boys function
                for(int i = 0; i <= absam1; i++)
                {
                    // vector of recurrence info
                    const auto & am1info = lut::am_recur_map[i];

                    // number of cartesians in the previous two shells
                    const size_t incart   = n_cartesian_gaussian(i);
                    const size_t incart_1 = (i > 0) ? n_cartesian_gaussian(i-1) : 0;
                    const size_t incart_2 = (i > 1) ? n_cartesian_gaussian(i-2) : 0;

                    // only form if i != 0 (ie, don't do (0,0)
                    if(i > 0)
                    {
                        // form (i,0)
                        double * iwork = amwork_[i][0];
                        double * iwork14 = amwork_[i-1][0];                      // location of the 1st and 4th terms
                        double * iwork25 = (i > 1) ? amwork_[i-2][0] : nullptr;  // location of the 2nd and 5th terms

     
                        // maximum value of (m) to calculate
                        // we need [0, absam12-i] inclusive
                        const int max_m = absam12 - i; 
                        size_t idx = 0;
                        for(int m = 0; m <= max_m; m++)
                        {
                            // commonly used in dimensioning
                            const size_t offset_1 = m*incart_1;
                            const size_t offset_4 = offset_1 + incart_1;  // (m+1)*incart_1
                            const size_t offset_2 = m*incart_2;
                            const size_t offset_5 = offset_2 + incart_2;  // (m+1)*incart_1

                            for(const auto & inf : am1info)
                            {
                                // get the recurrence information for this cartesian
                                const auto d = inf.dir;
                                const auto i_ijk = inf.ijk[d];
                                const size_t idx1 = offset_1 + inf.idx[d][0]; // index for 1st term
                                const size_t idx4 = offset_4 + inf.idx[d][0]; // index for 4th term
                                const size_t idx2 = offset_2 + inf.idx[d][1]; // index for 2nd term
                                const size_t idx5 = offset_5 + inf.idx[d][1]; // index for 5th term

                                iwork[idx] = PA[d]*iwork14[idx1] - PC[d]*iwork14[idx4];  // 1st and 4th terms

                                if(i_ijk > 1)
                                    iwork[idx] += oo2p*(i_ijk-1)*(iwork25[idx2] - iwork25[idx5]); // 2nd and 5th terms

                                idx++;
                            }
                        }
                    }

                    // now (i,j) via second vertical recurrence
                    for(int j = 1; j <= absam2; j++)
                    {
                        // vector of recurrence info
                        const auto & am2info = lut::am_recur_map[j];

                        // number of cartesians in the previous two shells
                        //const size_t jncart   = n_cartesian_gaussian(j);
                        const size_t jncart_1 = n_cartesian_gaussian(j-1);   // j can't be zero (loop starts at 1)
                        const size_t jncart_2 = (j > 1) ? n_cartesian_gaussian(j-2) : 0;

                        double * jwork = amwork_[i][j];

                        double * jwork14 = amwork_[i][j-1];                        // location of the 1st and 4th terms
                        double * jwork36 = (j > 1) ? amwork_[i][j-2] : nullptr;    // location of the 3rd and 6th terms
                        double * jwork25 = (i > 0) ? amwork_[i-1][j-1] : nullptr;  // location of the 2nd and 5th terms

                        const int max_m2 = absam2 - j + 1; // need [0, absam2-1] inclusive

                        size_t cartidx = 0; // index of the pair of cartesians
                        for(int m = 0; m <= max_m2; m++)
                        {
                            size_t cartidx_1 = 0;  // index of just the first cartesian
                            for(const auto & cart1 : am1info)
                            {
                                // precompute some of the offsets
                                // storage is  m, cart1, cart2
                                // so total index would be (m*ncart1*ncart2 + cart1*ncart2 + cart2)
                                const size_t offset1 = jncart_1*(m*incart + cartidx_1);                     // m*incart*jncart_1 + cartidx_1*jncart_1 
                                const size_t offset4 = jncart_1*((m+1)*incart + cartidx_1);                 // (m+1)*incart*jncart_1 + cartidx_1*jncart_1
                                const size_t offset3 = jncart_2*(m*incart + cartidx_1);                     // m*incart*jncart_2 + cartidx_1*jncart_2
                                const size_t offset6 = jncart_2*((m+1)*incart + cartidx_1);                 // (m+1)*incart*jncart_2 + cartidx_1*jncart_2
                                const size_t offset2[3] = { jncart_1*(m*incart_1 + cart1.idx[0][0]),        // m*incart_1*jncart_1 + cart1.idx[0][0]*jncart_1 
                                                            jncart_1*(m*incart_1 + cart1.idx[1][0]),        // m*incart_1*jncart_1 + cart1.idx[1][0]*jncart_1 
                                                            jncart_1*(m*incart_1 + cart1.idx[2][0]) };      // m*incart_1*jncart_1 + cart1.idx[2][0]*jncart_1 
                                const size_t offset5[3] = { jncart_1*((m+1)*incart_1 + cart1.idx[0][0]),    // (m+1)*incart_1*jncart_1 + cart1.idx[0][0]*jncart_1 
                                                            jncart_1*((m+1)*incart_1 + cart1.idx[1][0]),    // (m+1)*incart_1*jncart_1 + cart1.idx[1][0]*jncart_1 
                                                            jncart_1*((m+1)*incart_1 + cart1.idx[2][0]) };  // (m+1)*incart_1*jncart_1 + cart1.idx[2][0]*jncart_1 

                                for(const auto & cart2 : am2info)
                                {
                                    const auto d = cart2.dir; // direction we should recurse
                                    const auto i_ijk = cart1.ijk[d];  // values of i and j in that direction
                                    const auto j_ijk = cart2.ijk[d];
                                    const size_t idx1 = offset1 + cart2.idx[d][0];     // 1st term
                                    const size_t idx4 = offset4 + cart2.idx[d][0];     // 4th term
                                    const size_t idx3 = offset3 + cart2.idx[d][1];     // 3rd term
                                    const size_t idx6 = offset6 + cart2.idx[d][1];     // 6th term
                                    const size_t idx2 = offset2[d] + cart2.idx[d][0];  // 2nd term
                                    const size_t idx5 = offset5[d] + cart2.idx[d][0];  // 5th term


                                    jwork[cartidx] = PB[d]*jwork14[idx1] - PC[d]*jwork14[idx4]; // terms 1 & 4

                                    if(i_ijk > 0)
                                        jwork[cartidx] += oo2p*(i_ijk)*jwork25[idx2] - oo2p*(i_ijk)*jwork25[idx5]; // terms 2 & 5

                                    if(j_ijk > 1)
                                        jwork[cartidx] += oo2p*(j_ijk-1)*(jwork36[idx3] - jwork36[idx6]); // terms 3 & 6

                                    cartidx++;
                                }

                                cartidx_1++;
                            }
                        } // end loop over m
                    } // end loop over j
                } // end loop over i

                // general contraction and combined am
                size_t outidx = 0;
                for(size_t g1 = 0; g1 < ngen1; g1++)
                for(size_t g2 = 0; g2 < ngen2; g2++)
                {
                    const int gam1 = sh1.general_am(g1);
                    const int gam2 = sh2.general_am(g2);

                    const size_t ncart1 = n_cartesian_gaussian(gam1);
                    const size_t ncart2 = n_cartesian_gaussian(gam2);
                    double const * const amptr = amwork_[gam1][gam2];


                    size_t cartidx = 0;
                    // go over the orderings for this AM
                    for(size_t i = 0; i < ncart1; i++)
                    for(size_t j = 0; j < ncart2; j++)
                    {
                        const double val = amptr[cartidx++];

                        // remember: a and b are indices of primitives
                        // Also, the subtraction takes care of the minus sign
                        sourcework_[outidx++] -= val * sh1.coef(g1, a) * sh2.coef(g2, b) * gridpt.value;
                    }
                }
            } // end loop over primitive a
        } // end loop over primitive b
    } // close loop over atoms

    // performs the spherical transform, if necessary
    CartesianToSpherical_2Center(sh1, sh2, sourcework_, outbuffer, transformwork_, 1);

    return nfunc;
}


uint64_t OSOneElectronPotential::calculate_(uint64_t shell1, uint64_t shell2,
                                            double * outbuffer, size_t bufsize)
{
    // what grid are we using?
    std::string gridopt = options().get<std::string>("grid");

    GridUniverse gu;

    if(gridopt == "ATOMS")
    {
        // create the grid from the system
        for(const auto & atom : *sys_)
            if(atom.Z != 0)
                gu.insert({atom.get_coords(), numeric_cast<double>(atom.Z)});
        
    }
    else
        throw PulsarException("Unknown grid", "gridopt", gridopt);

    Grid grid(std::make_shared<GridUniverse>(std::move(gu)), true);

    // will check sizes of buffer, etc
    return calculate_with_grid_(shell1, shell2, grid, outbuffer, bufsize);
}



void OSOneElectronPotential::initialize_(unsigned int deriv,
                                         const Wavefunction & wfn,
                                         const BasisSet & bs1,
                                         const BasisSet & bs2)
{
    if(deriv != 0)
        throw NotYetImplementedException("Not Yet Implemented: OSOneElectronPotential integral with deriv != 0");

    sys_ = wfn.system;

    // from common components
    bs1_ = NormalizeBasis(cache(), out, bs1);
    bs2_ = NormalizeBasis(cache(), out, bs2);

    ///////////////////////////////////////
    // Determine the size of the workspace
    ///////////////////////////////////////

    // storage size for each x,y,z component
    int max1 = bs1_->max_am();
    int max2 = bs2_->max_am();
    size_t worksize = 0;

    // This overestimates a bit
    for(int i = 0; i <= max1; i++)
    for(int j = 0; j <= max2; j++)
        worksize += n_cartesian_gaussian(i)*n_cartesian_gaussian(j);
    int maxm = max1+max2;
    worksize *= (maxm+1);

    // find the maximum number of cartesian functions, not including general contraction
    size_t maxsize1 = bs1_->max_property(n_cartesian_gaussian_for_shell_am);
    size_t maxsize2 = bs2_->max_property(n_cartesian_gaussian_for_shell_am);
    size_t transformwork_size = maxsize1 * maxsize2; 

    // find the maximum number of cartesian functions, including general contraction
    maxsize1 = bs1_->max_property(n_cartesian_gaussian_in_shell);
    maxsize2 = bs2_->max_property(n_cartesian_gaussian_in_shell);
    size_t sourcework_size = maxsize1 * maxsize2;

    // allocate all at once, then partition
    work_.resize(worksize + transformwork_size + sourcework_size);

    amwork_.resize(max1+1);
    for(auto & it : amwork_)
        it.resize(max2+1);

    double * ptr = work_.data();
    for(int i = 0; i <= max1; i++)
    for(int j = 0; j <= max2; j++)
    {
        amwork_[i][j] = ptr;
        ptr += n_cartesian_gaussian(i)*n_cartesian_gaussian(j)*(maxm+1);
    }

    transformwork_ = ptr;
    sourcework_ = transformwork_ + transformwork_size;
}


} // close namespace integrals
} // close namespace psr_modules
