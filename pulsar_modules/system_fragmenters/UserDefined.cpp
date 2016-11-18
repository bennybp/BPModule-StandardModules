/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <pulsar/exception/PulsarException.hpp>
#include "pulsar_modules/system_fragmenters/UserDefined.hpp"

using std::vector;
using std::string;
using namespace pulsar;

NMerSetType UserDefined::fragmentize_(const System & mol){
    NMerSetType NMers;
    vector<string> Names=options().get<vector<string>>("FRAGMENT_NAMES");
    vector<int> AtomsPerFrag=options().get<vector<int>>("ATOMS_PER_FRAG");
    vector<int> Frags=options().get<vector<int>>("FRAGMENTS");
    if(Names.size()!=AtomsPerFrag.size())
        throw pulsar::PulsarException("The number of names for your fragments must"
                " match the length of the array of fragment sizes",
                "NNames",Names.size(),
                 "NSizes",AtomsPerFrag.size());
    
    vector<Atom> Atoms;
    for(const Atom& AtomI: mol)Atoms.push_back(AtomI);
    System Empty(mol.as_universe(),false);
    for(size_t i=0,counter=0;i<Names.size();++i){
        NMerInfo NMer;
        NMer.sn.insert(i);
        NMer.nmer=Empty;
        for(int j=0;j<AtomsPerFrag[i];++j)
            NMer.nmer.insert(Atoms[Frags[counter++]]);
        NMers.emplace(NMer.sn,NMer);
    }
    //return make_nmers(NMers);
}

