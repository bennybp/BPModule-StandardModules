/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include<pulsar/system/CrystalFunctions.hpp>
#include<pulsar/math/CombItr.hpp>
#include<pulsar/math/PowerSetItr.hpp>
#include<pulsar/math/BLAS.hpp>
#include "SystemFragmenters/CrystalFragger.hpp"

using namespace pulsar::system;
using namespace pulsar::modulebase;
using namespace pulsar::math;
using std::vector;
using std::array;
using std::get;
using std::make_pair;
using Fragger_t=pulsar::modulemanager::ModulePtr<SystemFragmenter>;
using NMerSet_t=typename NMerSetType::value_type;//std::pair<SNType,NMerInfo>

NMerSetType RenumberUC(const NMerSetType& UC,size_t StartFrom){
    NMerSetType UCFrags;
    for(const NMerSet_t& FragI: UC){
        SNType SN;
        for(size_t i=0;i<FragI.first.size();++i)
            SN.insert(std::to_string(StartFrom++));
        NMerInfo NMer(FragI.second);
        NMer.SN=SN;
        UCFrags.insert(make_pair(SN,NMer));
    }
    return UCFrags;
}

vector<SNType> CrysMakeNMers(size_t N,
                      NMerSetType& Frags,
                      const NMerSetType& UCFrags, 
                      const NMerSetType& SCFrags,
                      double DistThresh){
    vector<SNType> FinalSNs;
    for(size_t i=1;i<=std::min(N,UCFrags.size());++i){
        CombItr<NMerSetType> UCTuples(UCFrags,i);
        while(UCTuples){
            NMerInfo NMer;
            for(const NMerSet_t& FragI: *UCTuples){
                if(NMer.SN.size()==0)NMer.NMer=FragI.second.NMer;
                else NMer.NMer+=FragI.second.NMer;
                NMer.SN.insert(FragI.first.begin(),FragI.first.end());
            } 
            ++UCTuples;
            if(i==N){
                Frags.insert(std::make_pair(NMer.SN,NMer));
                FinalSNs.push_back(NMer.SN);
                continue;
            }
            CombItr<NMerSetType> SCTuples(SCFrags,N-i);
            while(SCTuples){            
                NMerInfo Temp(NMer);
                for(const auto& FragJ: *SCTuples){
                    Temp.SN.insert(FragJ.first.begin(),FragJ.first.end());
                    Temp.NMer+=FragJ.second.NMer;
                }
                ++SCTuples;
                FinalSNs.push_back(Temp.SN);
                Frags.insert(make_pair(Temp.SN,Temp));
            }            
        }
    }
    return FinalSNs;
}

NMerSetType GetUniqueNMers(NMerSetType& Frags,
                           const vector<SNType>& FinalSNs){
    std::map<SNType,SVDReturn_t> SVDs;
    for(const NMerSet_t& NMerI: Frags){
        const System& NMerJ=NMerI.second.NMer;
        std::vector<double> Carts=
            ToDoubleStar(NMerJ.Translate(-1.0*NMerJ.CenterOfMass()));       
        SVDs.emplace(NMerI.first,SVD(Carts,Carts.size()/3,3L));
    }
    CombItr<vector<SNType>> NMers(FinalSNs,2);
    while(NMers){
        NMerInfo& NMerI=Frags.at((*NMers)[0]);
        NMerInfo& NMerJ=Frags.at((*NMers)[1]);
        ++NMers;
        if(AreEqual(NMerI.Weight,0.0,1E-6)||AreEqual(NMerJ.Weight,0.0,1E-6))
            continue;
        bool Equal=true;
        for(size_t i=0;i<3;++i){
            if(AreEqual(get<1>(SVDs.at(NMerI.SN))[i],
                        get<1>(SVDs.at(NMerJ.SN))[i],2.0))continue;
            Equal=false;
            break;
        }
        if(Equal){
            NMerJ.Weight+=NMerI.Weight;
            NMerI.Weight=0.0;
        }
    }
    NMerSetType GoodFrags;
    //RMR copy_if is for some reason trying to copy to a const version...
    for(const NMerSet_t FragI: Frags)
        if(!AreEqual(FragI.second.Weight,0.0,1E-6))GoodFrags.insert(FragI);
    return GoodFrags;
}


NMerSetType CrystalFragger::Fragmentize_(const System & mol){
    NMerSetType Frags;
    const array<double,3>& Sides=mol.GetSpace().LatticeSides;
    const vector<size_t> Ls=Options().Get<vector<size_t>>("LATTICE_SIZE");
    const array<size_t,3> LatticeSides={Ls[0],Ls[1],Ls[2]};
    System UC=mol.Translate(Sides);
    System SC(MakeSuperCell(mol.AsUniverse(),LatticeSides,Sides),true);
    UC=SC.Partition([&](const Atom& a){return UC.Contains(a);});
    Fragger_t Fragger=
            CreateChildFromOption<SystemFragmenter>("SYSTEM_FRAGMENTER_KEY");
    SC-=UC;
    NMerSetType SCFrags=Fragger->Fragmentize(SC),
                UCFrags=Fragger->Fragmentize(UC);
    UCFrags=RenumberUC(UCFrags,SCFrags.size());
    const size_t N=Options().Get<size_t>("TRUNCATION_ORDER");
    double Dist=0.0;
    vector<SNType> FinalSNs=CrysMakeNMers(N,Frags,UCFrags,SCFrags,Dist);
    Frags=GetUniqueNMers(Frags,FinalSNs);
    vector<NMerSetType> Ints(N-1);
    vector<vector<SNType>> IntSNs(N-1);
    for(const NMerSet_t& NMerI: Frags){
        PowerSetItr<SNType> NMerJ(NMerI.first,1,N-1);
        while(NMerJ){
            size_t n=NMerJ->size()-1;
            if(Ints[n].count(*NMerJ)!=1){
                NMerInfo NewNMer;
                NewNMer.NMer=System(SC.GetUniverse(),false);
                NewNMer.Weight=0.0;
                NewNMer.SN=*NMerJ;
                for(const std::string& Name: *NMerJ){
                    NMerSetType* DaFrags=&SCFrags;
                    if(UCFrags.count({Name})==1)DaFrags=&UCFrags;
                    NewNMer.NMer+=DaFrags->at({Name}).NMer;
                }
                Ints[n].insert(std::make_pair(NewNMer.SN,NewNMer));
                IntSNs[n].push_back(NewNMer.SN);
            }
            double coef=-1.0;
            if((n+1)%2==N%2)coef=1.0;
            Ints[n][*NMerJ].Weight+=coef*NMerI.second.Weight;
            ++NMerJ;
        }
    }
    for(size_t i=0;i<N-1;++i){    
        NMerSetType Good=GetUniqueNMers(Ints[i],IntSNs[i]);
        Frags.insert(Good.begin(),Good.end());
    }
    return Frags;
}
