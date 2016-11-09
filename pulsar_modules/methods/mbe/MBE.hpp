#pragma once

#include <pulsar/modulebase/EnergyMethod.hpp>

namespace pulsarmethods{

/** \brief A common routine for many-body expansions
 
 */
class MBE : public pulsar::EnergyMethod{
    private:
        using Base_t=pulsar::EnergyMethod;
        using Wfn_t=pulsar::Wavefunction;
    public:
        //Uses constructor of base class
        using Base_t::EnergyMethod;
        ///Returns the \p Order -th derivative of the MBE of system in \p wfn
        pulsar::DerivReturnType deriv_(size_t Order,const Wfn_t& wfn);

};

}//End namespace pulsarmethods

