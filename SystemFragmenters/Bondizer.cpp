/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <pulsar/system/System.hpp>
#include <pulsar/datastore/OptionMap.hpp>
#include "SystemFragmenters/Bondizer.hpp"

using namespace pulsar::system;
using namespace pulsar::modulebase;

void Recurse(size_t depth,System& NewFrag,const Atom& CurrAtom,
             size_t MaxBonds,const Conn_t& Conns){
    if(MaxBonds==depth)return;
    for(const Atom& AnAtom: Conns.at(CurrAtom)){
        if(NewFrag.Contains(AnAtom))continue;//Already knew about it
        NewFrag.Insert(AnAtom);
        Recurse(depth+1,NewFrag,AnAtom,MaxBonds,Conns);
    }//End loop over attached atoms
}//End recurse

NMerSetType Bondizer::Fragmentize_(const System & mol){
    Conn_t Conns=GetConns(mol);
    size_t MaxBonds=Options().Get<size_t>("NBONDS");
    NMerSetType Frags;
    for(const Atom& AnAtom: mol){
        System NewFrag(mol.GetUniverse(),false);
        NewFrag.Insert(AnAtom);
        Recurse(1,NewFrag,AnAtom,MaxBonds,Conns);
        NewFrag.SetMultiplicity(mol.GetMultiplicity());
        bool good=true;
        for(const auto& OldFrags: Frags)//Check if unique
            if(OldFrags.second.NMer==NewFrag){
                good=false;
                break;
            }
        if(!good)continue;
        NMerInfo Temp;
        Temp.SN.insert(std::to_string(Frags.size()));
        Temp.NMer=NewFrag;
        Frags.emplace(Temp.SN,Temp);
    }
    return MakeNMers(Frags);
}
