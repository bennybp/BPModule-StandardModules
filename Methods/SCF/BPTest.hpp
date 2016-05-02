#ifndef BPTEST_HPP
#define BPTEST_HPP

#include <vector>
#include <pulsar/modulebase/EnergyMethod.hpp>
#include <pulsar/modulebase/All.hpp>

namespace pulsarmethods {

class BPTest : public pulsar::modulebase::EnergyMethod
{
    public:
        BPTest(ID_t id);
        
        virtual std::vector<double> Deriv_(size_t order);
};

}

#endif

