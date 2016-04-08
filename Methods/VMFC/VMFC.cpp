/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <unordered_map>
#include <bpmodule/system/System.hpp>
#include <bpmodule/math/PowerSetItr.hpp>
#include <bpmodule/modulebase/SystemFragmenter.hpp>
#include "Methods/MBE/MBEUtils.hpp"
#include "Methods/VMFC/VMFC.hpp"

using std::vector;
using std::string;

using bpmodule::system::System;
using bpmodule::system::Atom;
using bpmodule::system::AtomSetUniverse;
using bpmodule::system::MakeGhost;
using bpmodule::system::IsGhost;
using bpmodule::system::SystemMap;
using bpmodule::math::PowerSetItr;
using bpmodule::datastore::OptionMap;
using bpmodule::modulebase::SystemFragmenter;

typedef std::shared_ptr<const AtomSetUniverse> SharedUniverse;
typedef bpmodule::modulemanager::ModulePtr<SystemFragmenter> Fragmenter_t;

namespace bpmethods{
    
   string MakeName(const SN_t& FullSN,const SN_t& SN){
        SN_t NewSN=SN,Ghosts;
        std::set_difference(FullSN.begin(),FullSN.end(),
                            SN.begin(),SN.end(),
                            std::inserter(Ghosts,Ghosts.begin()));
        std::stringstream ss;
        for(size_t i: Ghosts)ss<<"Gh("<<i<<")_";
        size_t counter=0;
        for(size_t i: SN){
            ss<<i;
            if(++counter<SN.size())ss<<"_";
        }
        return ss.str();
}
    
    vector<double> VMFC::Deriv_(size_t Order){
        const System& Mol=*Wfn().system;
        
        AtomSetUniverse NewUniv;
        size_t NAtoms=Mol.Size();
        
        
        //We use two loops so that the real atoms come first and then the ghosts
        //ultimately this will make derivatives easier...
        
        //Because of the index, we need this map
        std::unordered_map<Atom,Atom> Real2Ghost;
        for(const Atom AtomI: Mol)NewUniv<<AtomI;
        for(const Atom AtomI: Mol){
            Real2Ghost.emplace(AtomI,MakeGhost(NAtoms++,AtomI));
            NewUniv<<Real2Ghost.at(AtomI);
        }
        System NewU(NewUniv,true);

        const OptionMap& DaOptions=Options();
        string Fragmentizer=DaOptions.Get<string>("FRAGMENTIZER");
        Fragmenter_t Fragger=CreateChildModule<SystemFragmenter>(Fragmentizer);
        SystemMap NMers=Fragger->Fragmentize(Mol);
        SNList_t SNs=BinNMers(NMers);
        SystemMap AllFrags;
        std::map<string,double> Coeffs;
        for(size_t n=0;n<SNs.size();++n){
            for(const auto& NMerI: SNs[n]){
                const SN_t& FullSN=NMerI.first;
                string FullName=NMerI.second;
                PowerSetItr<SN_t> SubNMers(FullSN,1,n);
                const System& NewUFullSys=NMers.at(FullName);
                System FullNMer=
                   NewU.Partition([&NewUFullSys](const Atom& AtomI){
                   return NewUFullSys.Contains(AtomI);});
                AllFrags.emplace(FullName,FullNMer);
                Coeffs[FullName]=1.0;//Real monomers are always 1
                while(!SubNMers.Done()){//Ghost the sub calculations
                    const SN_t& SN=*SubNMers;
                    System OldUSubSys=NMers.at(SNs[SN.size()-1][SN]);
                    System SubNMerI=
                            NewU.Partition([&OldUSubSys]
                            (const Atom& AtomI){return OldUSubSys.Contains(AtomI);}
                            );
                    System Comp=FullNMer-SubNMerI;
                    for(const Atom AtomI: Comp)
                        SubNMerI.Insert(Real2Ghost.at(AtomI));
                    AllFrags.emplace(MakeName(FullSN,SN),SubNMerI);
                    ++SubNMers;
                }
            }
        }
        
        for(const auto& NMerI: AllFrags)
            std::cout<<NMerI.first<<std::endl<<NMerI.second<<std::endl;
        exit(1);       
       
    }
    
    
    
}//End namespace