#include "Atomizer.hpp"

using namespace bpmodule::system;

MoleculeMap Atomizer::Fragmentize_(const Molecule & mol)
{
    MoleculeMap ret;

    //! \todo best way to do this? Could we copy
    //        the molecule (universe only) and then
    //        insert?
    for(const auto & atom : mol)
    {
        auto idx = atom.GetIdx();
        std::string idxstr = std::to_string(atom.GetIdx());
        ret.emplace(idxstr, mol.Partition([idx](const Atom & a) { return a.GetIdx() == idx; }));
    }

    return ret;
}

