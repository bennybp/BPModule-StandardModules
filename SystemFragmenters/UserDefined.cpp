/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <pulsar/exception/Exceptions.hpp>
#include "SystemFragmenters/UserDefined.hpp"

using std::vector;
using std::string;
using namespace pulsar::system;
using pulsar::exception::GeneralException;
using namespace pulsar::modulebase;

NMerSetType UserDefined::Fragmentize_(const System & mol){
    NMerSetType NMers;
    vector<string> Names=Options().Get<vector<string>>("FRAGMENT_NAMES");
    vector<int> AtomsPerFrag=Options().Get<vector<int>>("ATOMS_PER_FRAG");
    vector<int> Frags=Options().Get<vector<int>>("FRAGMENTS");
    if(Names.size()!=AtomsPerFrag.size())
        throw GeneralException("The number of names for your fragments must"
                " match the length of the array of fragment sizes",
                "NNames",Names.size(),
                 "NSizes",AtomsPerFrag.size());
    
    vector<Atom> Atoms;
    for(const Atom& AtomI: mol)Atoms.push_back(AtomI);
    System Empty(mol.AsUniverse(),false);
    for(size_t i=0,counter=0;i<Names.size();++i){
        NMerInfo NMer;
        NMer.SN.insert(Names[i]);
        NMer.NMer=Empty;
        for(int j=0;j<AtomsPerFrag[i];++j)
            NMer.NMer.Insert(Atoms[Frags[counter++]]);
        NMers.emplace(NMer.SN,NMer);
    }
    return MakeNMers(NMers);
}

