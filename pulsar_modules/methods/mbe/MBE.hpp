#pragma once

#include <pulsar/modulebase/EnergyMethod.hpp>

class MBE : public pulsar::EnergyMethod{
public:
    //Uses constructor of base class
    using pulsar::EnergyMethod::EnergyMethod;
    ///Returns the \p Order -th derivative of the MBE of system in \p wfn
    pulsar::DerivReturnType deriv_(size_t Order,const pulsar::Wavefunction& wfn);

};
