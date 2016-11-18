/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include<pulsar/system/CrystalFunctions.hpp>
#include<pulsar/math/CombItr.hpp>
#include<pulsar/math/PowerSetItr.hpp>
#include "pulsar_modules/system_fragmenters/CrystalFragger.hpp"

using namespace pulsar;
using std::vector;
using std::array;
using std::get;
using std::make_pair;
using Fragger_t=pulsar::ModulePtr<SystemFragmenter>;
using NMerSet_t=typename NMerSetType::value_type;//std::pair<SNType,NMerInfo>

NMerSetType RenumberUC(const NMerSetType& UC,size_t StartFrom){
    NMerSetType UCFrags;
    for(const NMerSet_t& FragI: UC){
        SNType SN;
        for(size_t i=0;i<FragI.first.size();++i)
            SN.insert(StartFrom++);
        NMerInfo NMer(FragI.second);
        NMer.sn=SN;
        UCFrags.insert(make_pair(SN,NMer));
    }
    return UCFrags;
}

vector<SNType> Crysmake_nmers(size_t N,
                      NMerSetType& Frags,
                      const NMerSetType& UCFrags, 
                      const NMerSetType& SCFrags,
                      double /*DistThresh*/){
    vector<SNType> FinalSNs;
    for(size_t i=1;i<=std::min(N,UCFrags.size());++i){
        CombItr<NMerSetType> UCTuples(UCFrags,i);
        while(UCTuples){
            NMerInfo NMer;
            for(const NMerSet_t& FragI: *UCTuples){
                if(NMer.sn.size()==0)NMer.nmer=FragI.second.nmer;
                else NMer.nmer+=FragI.second.nmer;
                NMer.sn.insert(FragI.first.begin(),FragI.first.end());
            } 
            ++UCTuples;
            if(i==N){
                Frags.insert(std::make_pair(NMer.sn,NMer));
                FinalSNs.push_back(NMer.sn);
                continue;
            }
            CombItr<NMerSetType> SCTuples(SCFrags,N-i);
            while(SCTuples){            
                NMerInfo Temp(NMer);
                for(const auto& FragJ: *SCTuples){
                    Temp.sn.insert(FragJ.first.begin(),FragJ.first.end());
                    Temp.nmer+=FragJ.second.nmer;
                }
                ++SCTuples;
                FinalSNs.push_back(Temp.sn);
                Frags.insert(make_pair(Temp.sn,Temp));
            }            
        }
    }
    return FinalSNs;
}

NMerSetType GetUniqueNMers(NMerSetType& Frags,
                           const vector<SNType>& FinalSNs){
    /*std::map<SNType,SVDReturn_t> SVDs;
    for(const NMerSet_t& NMerI: Frags){
        const System& NMerJ=NMerI.second.NMer;
        std::vector<double> Carts=
            ToDoubleStar(NMerJ.translate(-1.0*NMerJ.center_of_mass()));       
        SVDs.emplace(NMerI.first,SVD(Carts,Carts.size()/3,3L));
    }
    CombItr<vector<SNType>> NMers(FinalSNs,2);
    while(NMers){
        NMerInfo& NMerI=Frags.at((*NMers)[0]);
        NMerInfo& NMerJ=Frags.at((*NMers)[1]);
        ++NMers;
        if(are_equal(NMerI.Weight,0.0,1E-6)||are_equal(NMerJ.Weight,0.0,1E-6))
            continue;
        bool Equal=true;
        for(size_t i=0;i<3;++i){
            if(are_equal(get<1>(SVDs.at(NMerI.SN))[i],
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
        if(!are_equal(FragI.second.Weight,0.0,1E-6))GoodFrags.insert(FragI);
    return GoodFrags;*/
}


NMerSetType CrystalFragger::fragmentize_(const System & mol){
    NMerSetType Frags;
    const array<double,3>& Sides=mol.space.lattice_sides;
    const vector<size_t> Ls=options().get<vector<size_t>>("LATTICE_SIZE");
    const array<size_t,3> lattice_sides={Ls[0],Ls[1],Ls[2]};
    System UC=pulsar::translate(mol,Sides);
    System SC(MakeSuperCell(mol.as_universe(),lattice_sides,Sides),true);
    UC=SC.partition([&](const Atom& a){return UC.count(a);});
    Fragger_t Fragger=
            create_child_from_option<SystemFragmenter>("SYSTEM_FRAGMENTER_KEY");
    SC-=UC;
    NMerSetType SCFrags=Fragger->fragmentize(SC),
                UCFrags=Fragger->fragmentize(UC);
    UCFrags=RenumberUC(UCFrags,SCFrags.size());
    const size_t N=options().get<size_t>("TRUNCATION_ORDER");
    double Dist=0.0;
    vector<SNType> FinalSNs=Crysmake_nmers(N,Frags,UCFrags,SCFrags,Dist);
    Frags=GetUniqueNMers(Frags,FinalSNs);
    vector<NMerSetType> Ints(N-1);
    vector<vector<SNType>> IntSNs(N-1);
    System AFrag(*SC.get_universe(),false);
    for(const NMerSet_t& NMerI: Frags){
        PowerSetItr<SNType> NMerJ(NMerI.first,1,N-1);
        while(NMerJ){
            size_t n=NMerJ->size()-1;
            if(Ints[n].count(*NMerJ)!=1){
                NMerInfo NewNMer;
                NewNMer.nmer=AFrag;
                NewNMer.weight=0.0;
                NewNMer.sn=*NMerJ;
                for(size_t Name: *NMerJ){
                    NMerSetType* DaFrags=&SCFrags;
                    if(UCFrags.count({Name})==1)DaFrags=&UCFrags;
                    NewNMer.nmer+=DaFrags->at({Name}).nmer;
                }
                Ints[n].insert(std::make_pair(NewNMer.sn,NewNMer));
                IntSNs[n].push_back(NewNMer.sn);
            }
            double coef=-1.0;
            if((n+1)%2==N%2)coef=1.0;
            Ints[n][*NMerJ].weight+=coef*NMerI.second.weight;
            ++NMerJ;
        }
    }
    for(size_t i=0;i<N-1;++i){    
        NMerSetType Good=GetUniqueNMers(Ints[i],IntSNs[i]);
        Frags.insert(Good.begin(),Good.end());
    }
    return Frags;
}
