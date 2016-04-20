/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <pulsar/system/System.hpp>
#include <pulsar/datastore/OptionMap.hpp>
#include "SystemFragmenters/Bondizer.hpp"

using pulsar::datastore::OptionMap;
using pulsar::system::SystemMap;
using pulsar::system::System;
using pulsar::system::Atom;

typedef std::unordered_map<Atom,std::unordered_set<Atom>> Conn_t;

void Recurse(size_t depth,
             System& NewFrag,
             const Atom& CurrAtom,
             size_t MaxBonds,
             const Conn_t& Conns){
    if(MaxBonds==depth)return;
    for(const Atom& AnAtom: Conns.at(CurrAtom)){//Loop atoms attached to CurrAtm
        //Ensure we didn't find an atom we already knew about
        if(NewFrag.Contains(AnAtom))continue;
        NewFrag.Insert(AnAtom);//It's good
        Recurse(depth+1,NewFrag,AnAtom,MaxBonds,Conns);
    }//End loop over attached atoms
}//End recurse

SystemMap Bondizer::Fragmentize_(const System & mol){
    Conn_t Conns=pulsar::system::GetConns(mol);
    const OptionMap& DaOptions=Options();
    size_t MaxBonds=DaOptions.Get<size_t>("NBONDS");
    SystemMap Frags;
    for(const Atom& AnAtom: mol){
        System NewFrag=mol.Partition(
                [&AnAtom](const Atom& RHS){return AnAtom==RHS;}
                );
        NewFrag.SetMultiplicity(mol.GetMultiplicity());
        Recurse(1,NewFrag,AnAtom,MaxBonds,Conns);
        NewFrag.SetCharge(NewFrag.GetSumCharge());
        NewFrag.SetNElectrons(NewFrag.GetSumNElectrons());
        bool good=true;
        for(const auto& OldFrags: Frags){//Check if unique
            if(OldFrags.second==NewFrag)good=false;
            if(!good)break;
        }
        std::stringstream Name;
        Name<<Frags.size();
        if(good)Frags.emplace(Name.str(),NewFrag);
    }
    return Frags;
}
