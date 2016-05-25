/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Methods/BSSE/CP.hpp"
#include "Methods/BSSE/BSSEKernel.hpp"
#include "Methods/MBE/MBECommon.hpp"

namespace pulsarmethods{
    
using std::vector;    
using std::string;    

EnergyMethod::DerivReturnType CP::Deriv_(size_t Order,const Wavefunction& Wfn){
        const System& Mol=*Wfn.system;
        RealGhostData NewSys=GhostTheSystem(Mol);
        Fragmenter_t Fragger=CreateChild<SystemFragmenter>(
                                   Options().Get<string>("FRAGMENTIZER"));
        SystemMap NMers=Fragger->Fragmentize(Mol);        
        SNList_t SNs=BinNMers(NMers);
        size_t NFrags=SNs[0].size();
        
        //Make supersytem SN and name

        std::stringstream ss;
        for(size_t i=0;i<NFrags;++i){
            ss<<i;
            if(i<NFrags-1)ss<<"_";
        }
        std::map<string,double> Coeffs=SSFCKernel(NMers,NewSys,1,1,
                                                  NFrags-1,NFrags-1);
        
        NMers.emplace(ss.str(),*NewSys.RealSystem);
        Coeffs.emplace(ss.str(),1.0);
                
        return RunCalcs(NMers,Wfn,Coeffs,NewSys,Order,ID(),MManager(),
                        Options().Get<string>("METHOD"),
                        Options().Get<string>("MIM_KEY"));
        
        
        
        
}    
    
    
}//End namespace
