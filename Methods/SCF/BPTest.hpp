#ifndef BPTEST_HPP
#define BPTEST_HPP

#include <pulsar/modulebase/EnergyMethod.hpp>

namespace pulsarmethods {

class BPTest : public pulsar::modulebase::EnergyMethod
{
    public:
        using pulsar::modulebase::EnergyMethod::EnergyMethod;
        
        virtual DerivReturnType Deriv_(size_t order, const pulsar::datastore::Wavefunction & wfn);
};

}

#endif

