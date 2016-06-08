/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include<pulsar/modulebase/SystemFragmenter.hpp>
#include<pulsar/modulemanager/ModuleManager.hpp>
#include<pulsar/util/IterTools.hpp>
#include "Methods/MBE/MBE.hpp"
#include "Methods/MethodHelpers/MethodHelpers.hpp"

/*
 * TODO once computations cache results, compute interactions by
 *   "Re-running" the computations and then assembling derivatives
 */

using std::vector;
using std::string;
using Wfn_t=pulsar::datastore::Wavefunction;
using pulsar::util::Range;
using namespace pulsar::modulebase;
using SFer=SystemFragmenter;
using Return_t=DerivReturnType;
using namespace pulsar::system;
using Map_t=std::unordered_map<Atom,size_t>;

namespace pulsarmethods{
    
    Return_t MBE::Deriv_(size_t Order,const Wfn_t& Wfn)
    {
        vector<string> Keys(1,Options().Get<string>("METHOD"));
        const System& Mol=*Wfn.system;
        NMerSetType NMers=
                CreateChildFromOption<SFer>("FRAGMENTIZER")->Fragmentize(Mol);
        vector<double> Cs;
        vector<Wfn_t> Wfns;
        for(const typename NMerSetType::value_type& NMerI : NMers){
            Cs.push_back(NMerI.second.Weight);
            Wfns.push_back(Wfn);
            Wfns.back().system=std::make_shared<System>(NMerI.second.NMer);
        }
        vector<Return_t> Results=
                RunSeriesOfMethods(MManager(),ID(),Keys,Wfns,Cs,Order);
        Map_t SprMap;
        vector<double> Result(std::pow(3*Mol.size(),Order));
        for(const Atom& AtomI: Mol)SprMap.insert({AtomI,SprMap.size()});
        const size_t Offset=SprMap.size();
        for(const Atom& AtomI: Mol)
            SprMap.insert({MakeGhost(AtomI),SprMap.size()-Offset});
        
        for(size_t i: Range<0>(Results.size())){
            Map_t SubMap;
            const System& NMerI=*Results[i].first.system;
            for(const Atom& AtomI: NMerI)
                SubMap.insert({AtomI,SubMap.size()});
            FillDeriv(Result,Results[i].second,Cs[i],NMerI,SprMap,SubMap,Order);
        }
        return {Wfn,Result};
    }
    

}//End namespace
