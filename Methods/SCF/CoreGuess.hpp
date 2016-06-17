#ifndef PULSAR_GUARD_SCF__COREGUESS_HPP_
#define PULSAR_GUARD_SCF__COREGUESS_HPP_

#include <pulsar/modulebase/EnergyMethod.hpp>

namespace pulsarmethods {

class CoreGuess : public pulsar::modulebase::EnergyMethod
{
    public:
        using pulsar::modulebase::EnergyMethod::EnergyMethod; 
        
        virtual DerivReturnType Deriv_(size_t order, const pulsar::datastore::Wavefunction & wfn);
};

}

#endif

