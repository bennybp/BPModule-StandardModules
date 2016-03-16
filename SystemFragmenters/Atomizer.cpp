#include "Atomizer.hpp"

using namespace bpmodule::system;

SystemMap Atomizer::Fragmentize_(const System & mol)
{
    SystemMap ret;

    //! \todo best way to do this? Could we copy
    //        the System (universe only) and then
    //        insert?
    for(const auto & atom : mol)
    {
        auto idx = atom.GetIdx();
        std::string idxstr = std::to_string(atom.GetIdx());
        ret.emplace(idxstr, mol.Partition([idx](const Atom & a) { return a.GetIdx() == idx; }));
    }

    return ret;
}

