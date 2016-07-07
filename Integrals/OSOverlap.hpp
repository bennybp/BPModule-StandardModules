#pragma once

#include <pulsar/modulebase/OneElectronIntegral.hpp>

namespace psr_modules {
namespace integrals {

/*! \brief Calculation of overlap integrals via Obara-Saika recurrence
 */
class OSOverlap : public pulsar::modulebase::OneElectronIntegral
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

        double * transformwork_;
        double * sourcework_;
        double * xyzwork_[3];

        std::shared_ptr<pulsar::system::BasisSet> bs1_, bs2_;
};


} // close namespace integrals
} // close namespace psr_modules
