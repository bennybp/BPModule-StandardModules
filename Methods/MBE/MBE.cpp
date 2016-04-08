/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <string>
#include <bpmodule/system/System.hpp>
#include <bpmodule/modulebase/SystemFragmenter.hpp>
#include <bpmodule/datastore/OptionMap.hpp>
#include <bpmodule/math/PowerSetItr.hpp>
#include <bpmodule/math/Binomial.hpp>
#include "Methods/MBE/MBE.hpp"

using std::string;
using std::vector;
using std::map;
using std::set;
using std::stringstream;

using bpmodule::datastore::OptionMap;
using bpmodule::system::Atom;
using bpmodule::system::System;
using bpmodule::modulemanager::ModuleManager;
using bpmodule::modulebase::SystemFragmenter;
using bpmodule::modulebase::EnergyMethod;
using bpmodule::system::SystemMap;


typedef bpmodule::modulemanager::ModulePtr<SystemFragmenter> Fragmenter_t;
typedef bpmodule::modulemanager::ModulePtr<EnergyMethod> EMethod_t;
typedef vector<double> Return_t;
typedef vector<string> String_t;
typedef map<string,Return_t> DerivMap;
typedef set<size_t> SN_t;//NMer's serial number
typedef vector<map<SN_t,string>> SNList_t;//Type of binned SNs

/*
 * TODO once computations cache results, compute interactions by
 *   "Re-running" the computations and then assembling derivatives
 */

namespace bpmethods{

    
    //Splits the nmer names up
    String_t split(const string &s){
        String_t elems;
        stringstream ss(s);
        string item;
        while(std::getline(ss,item,'_'))elems.push_back(item);
        return elems;
    }
    
    //Computes the MBE coefficients by recursion
    void GetCoef(bool Even,const SN_t& NMer,const SNList_t& SNs,
                 map<string,double>& Coeffs){
        //If we don't have the NMer it's because we're assuming it's n-body
        //interaction is negligible, so don't follow that recursion
        size_t n=NMer.size()-1;
        if(SNs[n].count(NMer)==0)return;
        Coeffs[SNs[n].at(NMer)]+=(Even?1.0:-1.0);
        bpmodule::math::PowerSetItr<SN_t> Frags(NMer,1,n);//no empty sets
        while(!Frags.Done()){
            GetCoef(!Even,*Frags,SNs,Coeffs);
            ++Frags;
        }
    }
    
    Return_t MBE::DerivImpl(size_t Order)const{
        //Load options
        const OptionMap& DaOptions=Options();
        string MethodName=DaOptions.Get<string>("METHOD");
        string Fragmentizer=DaOptions.Get<string>("FRAGMENTIZER");

        //Make N-Mers
        const System& Mol=*(InitialWfn().system);
        Fragmenter_t Fragger=CreateChildModule<SystemFragmenter>(Fragmentizer);
        SystemMap NMers=Fragger->Fragmentize(Mol);
        
        //Parse the names, bin n-mers by size, zero weights
        SNList_t SNs;
        map<string,double> Weights;
        for(const auto& NMerI: NMers){
            const string Tag=NMerI.first;
            String_t StrBuffer=split(Tag);
            size_t N=StrBuffer.size();
            if(N>SNs.size())SNs.resize(N);
            SN_t SN;
            for(const string& Frag: StrBuffer)
                SN.insert(std::stoi(Frag));
            SNs[N-1][SN]=Tag;
            Weights[Tag]=0.0;
        }
        
        /* For the MBE we have the sum of all one-body interactions, the
         * sum of all two-body interactions, etc.  For each n-mer, GetCoef
         * adjusts the coefficients of the systems needed to get the energy
         * of the n-body interaction unique to that n-mer.  To get the
         * total we thus have to loop over all orders
         */
        for(size_t i=0;i<SNs.size();++i)
            for(const auto& NMerI: SNs[i])
                GetCoef(true,NMerI.first,SNs,Weights);        
        
        //Now we need the weights in the same order as the systems
        Return_t SortedWeights;
        for(const auto& NMerI: NMers)
            SortedWeights.push_back(Weights[NMerI.first]);
            
        
        //Pass info to MIM
        MManager().ChangeOption("MIM","METHODS",String_t({MethodName}));
        MManager().ChangeOption("MIM","FRAGMENTIZER",Fragmentizer);
        MManager().ChangeOption("MIM","WEIGHTS",SortedWeights);
        EMethod_t MIM=CreateChildModule<EnergyMethod>("MIM");
        return MIM->Deriv(Order);
    }
    

}//End namespace
