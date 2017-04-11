
#include<pulsar/modulebase/SystemFragmenter.hpp>
#include<pulsar/util/IterTools.hpp>
#include<pulsar/math/Cast.hpp>
#include "pulsar_modules/methods/mbe/MBE.hpp"
#include "pulsar_modules/methods/method_helpers/MethodHelpers.hpp"

//Options
const std::string method_key="METHOD_KEY";
const std::string system_key="SYSTEM_FRAGMENTER_KEY";

using namespace pulsar;
   
DerivReturnType MBE::deriv_(size_t Order,const Wavefunction& Wfn){
    std::vector<std::string> Keys={options().get<std::string>(method_key)};
    const System& Mol=*Wfn.system;
    NMerSetType NMers=
            create_child_from_option<SystemFragmenter>(system_key)->fragmentize(Mol);
    std::vector<double> Cs,Result(numeric_cast<size_t>(std::pow(3*Mol.size(),Order)));
    std::vector<Wavefunction> Wfns;
    for(const typename NMerSetType::value_type& NMerI : NMers){
        Cs.push_back(NMerI.second.weight);
        Wfns.push_back(Wfn);
        Wfns.back().system=std::make_shared<System>(
                                                  std::move(NMerI.second.nmer));
    }
    std::vector<DerivReturnType> Results=
            RunSeriesOfMethods(module_manager(),id(),Keys,Wfns,Order);
    
    std::unordered_map<Atom,size_t> SuperMap;
    
    for(const Atom& AtomI: Mol)SuperMap.insert({AtomI,SuperMap.size()});
    const size_t Offset=SuperMap.size();
    for(const Atom& AtomI: Mol)
        SuperMap.insert({make_ghost_atom(AtomI),SuperMap.size()-Offset});

    for(size_t i: Range<0>(Results.size())){
        std::unordered_map<Atom,size_t> SubMap;
        const System& NMerI=*Results[i].first.system;
        for(const Atom& AtomI: NMerI)
            SubMap.insert({AtomI,SubMap.size()});
        FillDeriv(Result,Results[i].second,Cs[i],NMerI,SuperMap,SubMap,Order);
    }
            
        return {Wfn,Result};
    }
