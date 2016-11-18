#include "pulsar_modules/system_fragmenters/Atomizer.hpp"
#include<iostream>
using namespace pulsar;

NMerSetType Atomizer::fragmentize_(const System & mol)
{
    NMerSetType ret;
    size_t idx = 0;
    for(const Atom & atom : mol)
    {
        NMerInfo Temp;
        Temp.sn.insert(idx);
        AtomSetUniverse DaAtom(atom);
        Temp.nmer=System(DaAtom,true);
        Temp.weight=1.0;
        ret.emplace(Temp.sn,Temp);
        idx++;
    }   
    
    return ret;
}

