#include "pulsar_modules/system_fragmenters/CPGhoster.hpp"

using namespace pulsar;

//Options
const std::string fragger_key="SYSTEM_FRAGMENTER_KEY";

NMerSetType CPGhoster::fragmentize_(const pulsar::System& mol){
    auto fragger=create_child_from_option<SystemFragmenter>(fragger_key);
    NMerSetType Frags=fragger->fragmentize(mol);
    AtomSetUniverse asu=mol.as_universe();
    for(const Atom& ai: mol)asu.insert(make_ghost_atom(ai));
    System NewMol(asu,false);
    for(auto& Fragi: Frags){
        System Temp(NewMol,false);
        for(const Atom& ai: mol)
            Temp.insert(Fragi.second.nmer.count(ai)?ai:make_ghost_atom(ai));
        Fragi.second.nmer=Temp;
    }
    return Frags;
    
    
}
