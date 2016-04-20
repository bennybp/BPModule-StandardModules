/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <unordered_map>
#include <bpmodule/math/PowerSetItr.hpp>
#include "Methods/MBE/MBEUtils.hpp"
#include "Methods/BSSE/VMFC.hpp"
#include "Methods/BSSE/BSSEKernel.hpp"
#include "Methods/MBE/MBECommon.hpp"

using std::vector;
using std::string;
using std::unordered_map;

using bpmodule::system::AtomSetUniverse;
using bpmodule::system::MakeGhost;
using bpmodule::system::IsGhost;
using bpmodule::math::PowerSetItr;

typedef std::shared_ptr<const AtomSetUniverse> SharedUniverse;
typedef vector<double> Return_t;

namespace bpmethods{  
    

    
   inline System SwitchBasis(const System& OldRealSys,const System& NewU){
      return NewU.Partition([&OldRealSys](const Atom& AtomI){
                            return OldRealSys.Contains(AtomI);});
   }
   
    vector<double> VMFC::Deriv_(size_t Order){
        const System& Mol=*InitialWfn().system;
        
        size_t NAtoms=Mol.Size();
        size_t DoF=1;
        for(size_t i=0;i<Order;++i)DoF*=3*NAtoms;
        
        RealGhostData NewSys=GhostTheSystem(Mol);
        System NewU(NewSys.NewSystem,true);      
       
        Fragmenter_t Fragger=CreateChildModule<SystemFragmenter>(
                                   Options().Get<string>("FRAGMENTIZER"));
        SystemMap NMers=Fragger->Fragmentize(Mol);
        
        SNList_t SNs=BinNMers(NMers);
        
        //Now we make all the little subsystems
        SystemMap AllFrags;
        std::map<string,double> Coeffs;
        for(size_t n=0;n<SNs.size();++n){//Loop over n in n-mers
            for(const auto& NMerI: SNs[n]){//Loop over n-mers
                const SN_t& RealSN=NMerI.first;
                string RealName=NMerI.second;
                PowerSetItr<SN_t> SubNMers(RealSN,1,n);
                //Add old molecule
                const System& RealNMer=SwitchBasis(NMers.at(RealName),NewU);
                AllFrags.emplace(RealName,RealNMer);
                Coeffs[RealName]=1.0;//Real monomers are always 1
                while(!SubNMers.Done()){//Ghost the sub calculations
                    const SN_t& SN=*SubNMers;
                    System SubNMerI=
                            SwitchBasis(NMers.at(SNs[SN.size()-1][SN]),NewU);
                    System Comp=RealNMer-SubNMerI;
                    for(const Atom AtomI: Comp)
                        SubNMerI.Insert(NewSys.Real2Ghost.at(AtomI));
                    string SubName=MakeName(RealSN,SN);
                    AllFrags.emplace(SubName,SubNMerI);
                    Coeffs[SubName]=(RealSN.size()%2==SN.size()%2?1.0:-1.0);
                    ++SubNMers;
                }
            }
        }

        return RunCalcs(AllFrags,Coeffs,NewSys,Order,ID(),MManager(),
                        Options().Get<string>("METHOD"),
                        Options().Get<string>("MIM_KEY"));
        
    }
    
    
    
}//End namespace