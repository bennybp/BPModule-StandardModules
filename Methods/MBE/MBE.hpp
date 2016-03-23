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
#include <bpmodule/modulebase/EnergyMethod.hpp>

namespace bpmethods{
class MBE : public bpmodule::modulebase::EnergyMethod{
    private:
        typedef bpmodule::modulebase::EnergyMethod Base_t;
        std::vector<double> DerivImpl(size_t Order)const;
    public:
        //Pull in energy method's constructors
        using Base_t::EnergyMethod;
        
        virtual std::vector<double> Deriv_(size_t Order){
            if(Order>0)return Base_t::Deriv_(Order);
            return DerivImpl(Order);
        }

};
}

#endif /* MBE_HPP */
