/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MBE.hpp
 * Author: richard
 *
 * Created on March 7, 2016, 11:06 AM
 */

#ifndef MBE_HPP
#define MBE_HPP

#include <pulsar/modulebase/EnergyMethod.hpp>

namespace pulsarmethods{

/** \brief A common routine for many-body expansions
 
 */
class MBE : public pulsar::modulebase::EnergyMethod{
    private:
        using Base_t=pulsar::modulebase::EnergyMethod;
        using Wfn_t=pulsar::datastore::Wavefunction;
    public:
        //Uses constructor of base class
        using Base_t::EnergyMethod;
        ///Returns the \p Order -th derivative of the MBE of system in \p wfn
        DerivReturnType deriv_(size_t Order,const Wfn_t& wfn);

};

}//End namespace pulsarmethods

#endif /* MBE_HPP */

