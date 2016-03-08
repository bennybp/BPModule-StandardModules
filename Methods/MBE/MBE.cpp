/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "MBE.hpp"
#include "bpmodule/system/Molecule.hpp"

namespace LibMBE{
    
    std::vector<double> MBE::DerivImpl(size_t Order)const{
        std::cout<<"Taking the "<<Order
                 <<"-th derivative of a Many-Body Expansion."<<std::endl;
        const bpmodule::system::Molecule& Mol=*Wfn().system;
        
        double sum=0.0;
        for(const auto& atom: Mol)
            for(size_t i=0;i<3;i++)sum+=atom[i]*atom[i];
        return std::vector<double>(1,sum);
    }
    
}