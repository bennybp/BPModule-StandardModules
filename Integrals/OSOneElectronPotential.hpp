#pragma once

#include <pulsar/modulebase/OneElectronIntegral.hpp>
#include <pulsar/math/Grid.hpp>

namespace psr_modules {
namespace integrals {


/*! \brief Calculation of one-electron potential integrals via Obara-Saika recurrence
 */
class OSOneElectronPotential : public pulsar::modulebase::OneElectronIntegral
{
    public:
        using pulsar::modulebase::OneElectronIntegral::OneElectronIntegral;

        virtual void initialize_(unsigned int deriv,
                                 const pulsar::datastore::Wavefunction & wfn,
                                 const pulsar::system::BasisSet & bs1,
                                 const pulsar::system::BasisSet & bs2);

        virtual uint64_t calculate_(uint64_t shell1, uint64_t shell2,
                                    double * outbuffer, size_t bufsize);

    private:
        std::vector<double> work_;

        // amwork_[i][j] = work for am pair i,j
        std::vector<std::vector<double *>> amwork_;
        std::shared_ptr<const pulsar::system::System> sys_;

        double * transformwork_;
        double * sourcework_;

        std::shared_ptr<pulsar::system::BasisSet> bs1_, bs2_;

        uint64_t calculate_with_grid_(uint64_t shell1, uint64_t shell2,
                                      const pulsar::math::Grid & grid,
                                      double * outbuffer, size_t bufsize);
};


} // close namespace integrals
} // close namespace psr_modules
