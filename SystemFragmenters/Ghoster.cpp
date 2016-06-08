/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include<pulsar/exception/Exceptions.hpp>
#include<pulsar/math/PowerSetItr.hpp>
#include "SystemFragmenters/Ghoster.hpp"

using std::string;
using TruncOrder_t=std::map<size_t,size_t>;
using namespace::pulsar::system;
using namespace::pulsar::modulebase;
using namespace::pulsar::math;
using R2G_t=std::unordered_map<Atom,Atom>;
using NMerSet_t=typename NMerSetType::value_type;

const string FraggerKey="SYSTEM_FRAGMENTER_KEY";
const string GhTruncKey="GHOST_TRUNCATION_ORDERS";


/************  The inline functions below are defined to simplify understanding
               of the code.
 */

//Determines N and n given a set of fragments
inline std::pair<size_t,size_t> NMerStats(const NMerSetType& OrigFrags){
    size_t NFrags=0,MaxN=0;
    for(const NMerSet_t& NMerI: OrigFrags){
        const size_t n=NMerI.first.size();
        if(n==1)++NFrags;
        MaxN=std::max(MaxN,n);
    }
    return {NFrags,MaxN};
}

//Wrapper around a semi default NMer initialization
inline NMerSet_t InitializeNMerSet_t(const NMerSet_t& NMerI,const System& NewS){
    NMerSet_t NMerJ;
    NMerJ.second.NMer=System(NewS);
    NMerJ.second.NMer.insert(NMerI.second.NMer.begin(),NMerI.second.NMer.end());
    NMerJ.second.SN.insert(NMerI.first.begin(),NMerI.first.end());
    NMerJ.second.Weight=NMerI.second.Weight;
    return {NMerJ.second.SN,NMerJ.second};
}

//Wrapper around determining the set of active fragments
inline SNType ActiveSN(const SNType& FullSN,const NMerSet_t& NMerI){
    SNType TempSN;
    std::set_difference(FullSN.begin(),FullSN.end(),
                        NMerI.first.begin(),NMerI.first.end(),
                        std::inserter(TempSN,TempSN.end()));
    return TempSN;
}

NMerSetType Ghoster::Fragmentize_(const System& mol){
    NMerSetType OrigFrags=
          CreateChildFromOption<SystemFragmenter>(FraggerKey)->Fragmentize(mol);
    
    //All fragments have same universe so just grab one
    AtomSetUniverse NewU,OldU;
    NewU=OldU=*(OrigFrags.begin()->second.NMer.GetUniverse());
    R2G_t R2G;
    for(const Atom& AtomI : OldU)
        NewU.Insert(R2G.insert({AtomI,MakeGhost(AtomI)}).first->second);
    System NewSys(NewU,false);
    std::pair<size_t,size_t> Stats=NMerStats(OrigFrags);
    const size_t NFrags=Stats.first,MaxN=Stats.second;
    TruncOrder_t TruncOrders=Options().Get<TruncOrder_t>(GhTruncKey);
    bool OnlyFull=true,VMFCn=true;
    
    for(const auto& TO: TruncOrders){
        if(TO.second<NFrags-TO.first)OnlyFull=false;
        if(TO.second!=MaxN-TO.first)VMFCn=false;
    }

    if(!OnlyFull && ! VMFCn)
        throw pulsar::exception::GeneralException(
                "I have not coded up weights for arbitrary "
                "combinations of real and ghost monomers.");
    
    //Supersystem SN
    SNType FullSN;
    for(size_t i=0;i<NFrags;++i)FullSN.insert(std::to_string(i));
    
    //The frags to be returned
    NMerSetType NewFrags;
    for(const NMerSet_t& NMerI: OrigFrags){
        NMerSet_t NMerJ=InitializeNMerSet_t(NMerI,NewSys);
        SNType TempSN=ActiveSN(FullSN,NMerJ);
        size_t Min=0,Max=TempSN.size();
        const size_t r=NMerI.first.size();
        if(OnlyFull)Min=Max;
        else if(TruncOrders.count(r)==1)Max=TruncOrders.at(r);
        else Max=0;
        PowerSetItr<SNType> Ghosts(TempSN,Min,Max);
        while(Ghosts){
            NMerSet_t CopyOfI=NMerJ;
            for(const string& SNI: *Ghosts){
                CopyOfI.second.SN.insert(std::to_string(NFrags+std::atoi(SNI.c_str())));
                for(const Atom& AtomI: OrigFrags.at({SNI}).NMer)
                    CopyOfI.second.NMer.Insert(R2G.at(AtomI));
            }
            if(VMFCn)CopyOfI.second.Weight=(Ghosts->size()%2==0?1.0:-1.0);
            NewFrags.insert({CopyOfI.second.SN,CopyOfI.second});
            ++Ghosts;
        }
    }
    return NewFrags;
}

NMerSetType CommonGuts(const System& Mol,
                       bool UseSuper,
                       const SystemFragmenter* SF,
                       pulsar::modulemanager::ModuleManager& MM){
    using Fragger_t=pulsar::modulemanager::ModulePtr<SystemFragmenter>;
    const string GhostKey=SF->Options().Get<string>("GHOSTER_KEY");
    Fragger_t Fragger=SF->CreateChild<SystemFragmenter>(GhostKey);
    const string OrigKey=Fragger->Options().Get<string>(FraggerKey);
    NMerSetType OrigFrags=
            SF->CreateChild<SystemFragmenter>(OrigKey)->Fragmentize(Mol);
    std::pair<size_t,size_t> Stats=NMerStats(OrigFrags);
    const string NewKey=MM.GenerateUniqueKey();
    MM.DuplicateKey(GhostKey,NewKey);
    TruncOrder_t NewTruncs;
    const size_t Max=(UseSuper?Stats.first:Stats.second);
    for(size_t i=1;i<=Stats.second;++i)NewTruncs.insert({i,Max-i});
    MM.ChangeOption(NewKey,GhTruncKey,NewTruncs);
    return SF->CreateChild<SystemFragmenter>(NewKey)->Fragmentize(Mol);
}


NMerSetType CPGhoster::Fragmentize_(const pulsar::system::System& Mol){
    return CommonGuts(Mol,true,this,MManager());
}

NMerSetType VMFCGhoster::Fragmentize_(const pulsar::system::System& Mol){
    return CommonGuts(Mol,false,this,MManager());
}

