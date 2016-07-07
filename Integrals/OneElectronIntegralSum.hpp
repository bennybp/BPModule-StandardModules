#pragma once

#include <pulsar/modulebase/OneElectronIntegral.hpp>

namespace psr_modules {
namespace integrals {

/*! \brief Calculates a sum of one-electron integrals
 *
 * The keys for which modules to incorporate into the summation. Each call
 * to calculate() will in turn call these modules and sum the results into
 * the output buffer.
 */
class OneElectronIntegralSum : public pulsar::modulebase::OneElectronIntegral
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
        typedef pulsar::modulemanager::ModulePtr<pulsar::modulebase::OneElectronIntegral> OneInt;

        //! The modules to use in constructing the core hamiltonian
        std::map<std::string, OneInt> modules_;
};

} // close namespace integrals
} // close namespace psr_modules
