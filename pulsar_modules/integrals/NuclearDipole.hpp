#pragma once

#include <pulsar/modulebase/SystemIntegral.hpp>

namespace psr_modules {
namespace integrals {


/*! \brief Calculation of simple nuclear dipole integrals
 */
class NuclearDipole : public pulsar::modulebase::SystemIntegral
{
    public:
        using pulsar::modulebase::SystemIntegral::SystemIntegral;

        virtual void initialize_(unsigned int deriv, const pulsar::system::System & sys);

        virtual uint64_t calculate_(double * outbuffer, size_t bufsize);

    private:
        const pulsar::system::System * sys_;
};


} // close namespace integrals
} // close namespace psr_modules

