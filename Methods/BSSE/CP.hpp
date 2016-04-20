/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CP.hpp
 * Author: richard
 *
 * Created on April 4, 2016, 5:55 PM
 */

#ifndef CP_HPP
#define CP_HPP

#include <vector>
#include <pulsar/modulebase/EnergyMethod.hpp>

namespace pulsarmethods{

/** \brief The Boys and Bernardi CP correction
 *
 *  The traditional Boys and Bernardi CP correction to the energy of a system
 *  that has been decomposed into
 *  \f$N\f$ disjoint subsystems is computed via:
 *  \f[
 *     BSSE=\sum_i^N E_i^{Full}-E_i^{Normal}
 *  \f]
 *  where full is the full system's basis set and normal is whatever basis set
 *  the \f$i\f$-th system would normally have.
 *  This function returns the energy of the supersystem less those the
 *  BSSE (and derivatives there of).  Hence to get BSSE-free interaction
 *  energies subtract out the sum of the monomers in their normal basis set.
 *  We insist on returning BSSE-free energies to maintain a consistent
 *  interface.  This means this module will perform the supersystem computation.
 * 
 *  \todo Make a module for a MBE in the full supersystem's basis set
 */
class CP: public pulsar::modulebase::EnergyMethod{
private:
    typedef pulsar::modulebase::EnergyMethod Base_t;
public:
    using Base_t::EnergyMethod;
    virtual std::vector<double> Deriv_(size_t Order);
};

}//End namespace

#endif /* SSFC_HPP */

