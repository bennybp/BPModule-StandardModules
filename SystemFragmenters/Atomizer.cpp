#include "Atomizer.hpp"

using namespace pulsar::system;
using namespace pulsar::modulebase;

NMerSetType Atomizer::fragmentize_(const System & mol)
{
    NMerSetType ret;
    size_t idx = 0;
    for(const Atom & atom : mol)
    {
        std::string idxstr = std::to_string(idx);
        NMerInfo Temp;
        Temp.SN.insert(idxstr);
        AtomSetUniverse DaAtom(atom);
        Temp.NMer=System(DaAtom,true);
        ret.emplace(Temp.SN,Temp);
        idx++;
    }

    return make_nmers(ret);
}

