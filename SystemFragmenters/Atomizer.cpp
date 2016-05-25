#include "Atomizer.hpp"

using namespace pulsar::system;

SystemMap Atomizer::Fragmentize_(const System & mol)
{
    SystemMap ret;

    //! \todo best way to do this? Could we copy
    //        the System (universe only) and then
    //        insert?
    size_t idx = 0;
    for(const auto & atom : mol)
    {
        std::string idxstr = std::to_string(idx);
        auto hash = atom.MyHash();
        ret.emplace(idxstr, mol.Partition([hash](const Atom & a) { return a.MyHash() == hash; }));
        idx++;
    }

    return ret;
}

