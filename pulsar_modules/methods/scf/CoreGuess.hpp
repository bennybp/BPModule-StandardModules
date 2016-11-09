#ifndef PULSAR_GUARD_SCF__COREGUESS_HPP_
#define PULSAR_GUARD_SCF__COREGUESS_HPP_

#include <pulsar/modulebase/EnergyMethod.hpp>

namespace pulsarmethods {

class CoreGuess : public pulsar::EnergyMethod
{
    public:
        using pulsar::EnergyMethod::EnergyMethod; 
        
        virtual pulsar::DerivReturnType deriv_(size_t order, const pulsar::Wavefunction & wfn);
};

}

#endif

