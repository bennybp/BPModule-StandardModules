#include "pulsar_modules/system_fragmenters/Bondizer.hpp"
#include <algorithm>//For find_if

using std::find_if;
using namespace pulsar;
const std::string n_bonds_key="MAX_NBONDS";

void Recurse(size_t depth,System& NewFrag,const Atom& CurrAtom,
             size_t MaxBonds,const Conn_t& Conns){
    if(MaxBonds==depth)return;
    for(const Atom& AnAtom: Conns.at(CurrAtom)){
        if(NewFrag.count(AnAtom))continue;//Already knew about it
        NewFrag.insert(AnAtom);
        Recurse(depth+1,NewFrag,AnAtom,MaxBonds,Conns);
    }//End loop over attached atoms
}//End recurse

NMerSetType Bondizer::fragmentize_(const System & mol){
    const Conn_t Conns=get_connectivity(mol);
    const size_t MaxBonds=options().get<size_t>(n_bonds_key);
    NMerSetType Frags;
    for(const Atom& AnAtom: mol){
        NMerInfo Temp={{Frags.size()},System(mol,false),1.0};
        Temp.nmer.insert(AnAtom);
        Recurse(0,Temp.nmer,AnAtom,MaxBonds,Conns);
        
        //Check if we already found this Frag
        auto OldFrag=find_if(Frags.begin(),Frags.end(),
          [&Temp](const NMerSetType::value_type& Fi){
              return Fi.second.nmer==Temp.nmer;
          }
        );
        if(OldFrag!=Frags.end())continue;
        Frags.emplace(Temp.sn,std::move(Temp));
    }
    return Frags;
}
