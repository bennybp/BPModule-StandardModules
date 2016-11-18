/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include<pulsar/modulebase/SystemFragmenter.hpp>
#include<pulsar/modulemanager/ModuleManager.hpp>
#include<pulsar/util/IterTools.hpp>
#include "pulsar_modules/methods/mbe/MBE.hpp"
#include "pulsar_modules/methods/method_helpers/MethodHelpers.hpp"

/*
 * TODO once computations cache results, compute interactions by
 *   "Re-running" the computations and then assembling derivatives
 */

using std::vector;
using std::string;
using namespace pulsar;
using SFer=SystemFragmenter;
using Return_t=DerivReturnType;
using Map_t=std::unordered_map<Atom,size_t>;

namespace pulsarmethods{
    
    Return_t MBE::deriv_(size_t Order,const Wfn_t& Wfn)
    {
        vector<string> Keys={options().get<string>("METHOD")};
        const System& Mol=*Wfn.system;
        NMerSetType NMers=
                create_child_from_option<SFer>("FRAGMENTIZER")->fragmentize(Mol);
        vector<double> Cs;
        vector<Wfn_t> Wfns;
        for(const typename NMerSetType::value_type& NMerI : NMers){
            Cs.push_back(NMerI.second.weight);
            Wfns.push_back(Wfn);
            Wfns.back().system=std::make_shared<System>(NMerI.second.nmer);
        }
        vector<Return_t> Results=
                RunSeriesOfMethods(module_manager(),id(),Keys,Wfns,Order);
        Map_t SuperMap;
        vector<double> Result(std::pow(3*Mol.size(),Order));
        for(const Atom& AtomI: Mol)SuperMap.insert({AtomI,SuperMap.size()});
        const size_t Offset=SuperMap.size();
        for(const Atom& AtomI: Mol)
            SuperMap.insert({make_ghost_atom(AtomI),SuperMap.size()-Offset});
        
        for(size_t i: Range<0>(Results.size())){
            Map_t SubMap;
            const System& NMerI=*Results[i].first.system;
            for(const Atom& AtomI: NMerI)
                SubMap.insert({AtomI,SubMap.size()});
            FillDeriv(Result,Results[i].second,Cs[i],NMerI,SuperMap,SubMap,Order);
        }
            
        return {Wfn,Result};
    }
    
}//End namespace