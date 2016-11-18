#include "pulsar_modules/system_fragmenters/VMFCGhoster.hpp"
#include <pulsar/math/PowerSetItr.hpp>
using namespace pulsar;

//Options
const std::string fragger_key="SYSTEM_FRAGMENTER_KEY";

NMerSetType VMFCGhoster::fragmentize_(const pulsar::System& mol){
    auto fragger=create_child_from_option<SystemFragmenter>(fragger_key);
    NMerSetType Frags=fragger->fragmentize(mol);
    AtomSetUniverse asu=mol.as_universe();
    for(const Atom& ai: mol)asu.insert(make_ghost_atom(ai));
    System NewMol(asu,false);
    std::map<SNType,System> GhostFrags;
    
    //Get the ghost version of the frags
    for(auto& Fragi: Frags){
        if(Fragi.first.size()>1)continue;
        System Temp(NewMol,false);
        for(const Atom& ai: Fragi.second.nmer)Temp.insert(make_ghost_atom(ai));
        GhostFrags.emplace(Fragi.first,std::move(Temp));
    }
    
    const size_t NFrags=GhostFrags.size();
    
    //Unionize them
    NMerSetType NMers;
    for(auto& Fragi: Frags){
        NMers.insert(Fragi);
        const SNType& SNi=Fragi.first;
        NMers[SNi].weight=1.0;
        PowerSetItr<SNType> real(SNi,1,SNi.size()-1);
        while(real){
            SNType sn=*real;
            System Temp(NewMol,false);
            for(size_t i : Fragi.first){
                const bool is_real=real->count(i);
                Temp+=(is_real?Frags.at({i}).nmer:GhostFrags[{i}]);
                if(!is_real)sn.insert(i+NFrags);
            }
            const bool same_parity=SNi.size()%2==real->size()%2;
            NMers.emplace(sn,NMerInfo({sn,Temp,(same_parity?1.0:-1.0)}));
            ++real;
        }
    }
    return NMers;
}
