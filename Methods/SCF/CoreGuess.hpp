#ifndef COREGUESS_HPP
#define COREGUESS_HPP

#include <vector>
#include <pulsar/modulebase/All.hpp>

namespace pulsarmethods {

class CoreGuess : public pulsar::modulebase::EnergyMethod
{
    public:
        using pulsar::modulebase::EnergyMethod::EnergyMethod; 
        
        virtual std::vector<double> Deriv_(size_t order);
};

}

#endif

