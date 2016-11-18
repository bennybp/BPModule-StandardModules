#include "pulsar_modules/system_fragmenters/NMerizer.hpp"
#include "pulsar/math/PowerSetItr.hpp"
#include "pulsar/math/Cast.hpp"


using namespace pulsar;
/*************** Options ****************************/
const std::string fragger_key="SYSTEM_FRAGMENTER_KEY";
const std::string trunc_key="TRUNCATION_ORDER";
const std::string dist_key="DISTANCE_THRESHOLDS";


//Recursive function for establishing the weight of each nmer
//Note b/c we are possibly leaving some out we can't use the closed formulas
inline void GetCoef(bool Even,const SNType& SN,NMerSetType& NMers){
    if(NMers.count(SN)==0)return;//We assumed this interaction is negligible
    NMers[SN].weight+=(Even?1.0:-1.0);
    PowerSetItr<SNType> Frags(SN,1,numeric_cast<int>(SN.size())-1);
    while(Frags)GetCoef(!Even,*Frags++,NMers);
}

NMerSetType NMerizer::fragmentize_(const System & mol){
    auto fragger=create_child_from_option<SystemFragmenter>(fragger_key);
    NMerSetType NMers,Frags=fragger->fragmentize(mol);
    const int NEnd=numeric_cast<int>(
        std::min(options().get<size_t>(trunc_key),Frags.size()));
    
    //Make sure we work in two stupid base case scenarios
    if(NEnd==0)return NMers;//Empty set
    if(NEnd==1)return Frags;//The fragments we were given
    

    const auto Dists=options().get<std::map<size_t,double>>(dist_key);

    PowerSetItr<NMerSetType> Comb(Frags,1,NEnd);
    while(Comb){
        NMerInfo DaNMer={{},System(mol,false),0.0};
        for(const auto& Frag : *Comb){
            const NMerInfo& NMer=Frag.second;
            DaNMer.sn.insert(NMer.sn.begin(),NMer.sn.end());
            DaNMer.nmer.insert(NMer.nmer.begin(),NMer.nmer.end());
        }
        const size_t N=Comb->size();
        bool is_bad=Dists.count(N);
        if(is_bad){
            double RHS=1.0;
            const Point CoM=center_of_mass(DaNMer.nmer);
            for(const auto& Frag: *Comb)
                RHS*=center_of_mass(Frag.second.nmer).distance(CoM);
            RHS=std::pow(RHS,1.0/numeric_cast<double>(N));
            if(RHS>Dists.at(N))is_bad=true;
        }
        if(!is_bad)NMers.insert({DaNMer.sn,DaNMer});
        ++Comb;
    }
    for(const auto& NMerI: NMers)GetCoef(true,NMerI.first,NMers);
    return NMers;
}