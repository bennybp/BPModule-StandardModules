#ifndef HFITERATE_HPP
#define HFITERATE_HPP

#include <vector>
#include <pulsar/modulebase/EnergyMethod.hpp>
#include <pulsar/modulebase/All.hpp>

namespace pulsarmethods {

class HFIterate : public pulsar::modulebase::EnergyMethod
{
    public:
        using pulsar::modulebase::EnergyMethod::EnergyMethod;
        
        virtual DerivReturnType Deriv_(size_t order, const pulsar::datastore::Wavefunction & wfn);
};

}

#endif

