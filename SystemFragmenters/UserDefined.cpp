/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <pulsar/exception/Exceptions.hpp>
#include "SystemFragmenters/UserDefined.hpp"

using std::vector;
using std::string;

using pulsar::system::System;
using pulsar::system::Atom;
using pulsar::system::SystemMap;
using pulsar::datastore::OptionMap;
using pulsar::exception::GeneralException;

SystemMap UserDefined::Fragmentize_(const System & mol){
    SystemMap NMers;
    const OptionMap& DaOptions=Options();
    vector<string> Names=
            DaOptions.Get<vector<string>>("FRAGMENT_NAMES");
    vector<int> AtomsPerFrag=DaOptions.Get<vector<int>>("ATOMS_PER_FRAG");
    vector<int> Frags=DaOptions.Get<vector<int>>("FRAGMENTS");
    
    if(Names.size()!=AtomsPerFrag.size())
        throw GeneralException("The number of names for your fragments must"
                " match the length of the array of fragment sizes",
                "NNames",Names.size(),
                 "NSizes",AtomsPerFrag.size());
    
    std::vector<Atom> Atoms;
    for(const Atom& AtomI: mol)Atoms.push_back(AtomI);
        
    System Empty=mol.Partition([](const Atom&){return false;});
    
    for(size_t i=0,counter=0;i<Names.size();++i){
        NMers.emplace(Names[i],Empty);
        for(size_t j=0;j<AtomsPerFrag[i];++j)
            NMers.at(Names[i]).Insert(Atoms[Frags[counter++]]);
    }
    return NMers;
}