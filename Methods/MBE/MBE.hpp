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

#include <vector>
#include <pulsar/modulebase/EnergyMethod.hpp>

namespace pulsarmethods{
class MBE : public pulsar::modulebase::EnergyMethod{
    private:
        typedef pulsar::modulebase::EnergyMethod Base_t;
    public:
        //Pull in energy method's constructors
        using Base_t::EnergyMethod;
        
        virtual std::vector<double> Deriv_(size_t Order);

};
}

#endif /* MBE_HPP */

