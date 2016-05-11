#ifndef COREGUESS_HPP
#define COREGUESS_HPP

#include <vector>
#include <pulsar/modulebase/All.hpp>

namespace pulsarmethods {

class CoreGuess : public pulsar::modulebase::EnergyMethod
{
    public:
        using pulsar::modulebase::EnergyMethod::EnergyMethod; 
        
        virtual DerivReturnType Deriv_(size_t order, const pulsar::datastore::Wavefunction & wfn);
};

}

#endif

