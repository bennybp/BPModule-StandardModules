/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <sstream>
#include <pulsar/system/System.hpp>
#include <pulsar/system/Atom.hpp>
#include <pulsar/math/PowerSetItr.hpp>
#include "Methods/MBE/MBEUtils.hpp"

using std::stringstream;
using std::string;
using std::vector;
using std::map;
using std::unordered_map;

using pulsar::system::System;
using pulsar::system::Atom;
using pulsar::system::SystemMap;

typedef vector<double> Return_t;

namespace pulsarmethods{

//Splits the nmer names up
vector<string> split(const string &s){
    vector<string> elems;
    stringstream ss(s);
    string item;
    while(std::getline(ss,item,'_'))
        elems.push_back(item);
    return elems;
}
        
SNList_t BinNMers(const SystemMap& Frags){
    SNList_t SNs;
    for(const auto& NMerI: Frags){
        const string Tag=NMerI.first;
        vector<string> StrBuffer=split(Tag);
        size_t N=StrBuffer.size();
        if(N>SNs.size())SNs.resize(N);
        SN_t SN;
        for(const string& Frag: StrBuffer)
            SN.insert(std::stoi(Frag));
        SNs[N-1][SN]=Tag;
    }
    return SNs;
}

void FillDeriv(Return_t& Result, 
               const Return_t& SubResult,
               double Coeff,
               const System& Sys, 
               const unordered_map<Atom,size_t>& SuperAtomMap,
               const unordered_map<Atom,size_t>& SubAtomMap,
               size_t Order,
               vector<Atom> Idx,
               vector<size_t> Comp){
   //Handle energy
   if(Order==0)Result[0]+=Coeff*SubResult[0];
   else if(Idx.size()==(Order-1)){//Gradient lands here, Hessian and higher on recursion
      size_t SuperOff=0,SubOff=0;
      for(size_t i=0;i<Order-1;++i){//First offset
         size_t SuperStride=1,SubStride=1;
         for(size_t j=i;j<Order-1;++j){
            SuperStride*=3*SuperAtomMap.size();
            SubStride*=3*SubAtomMap.size();
         }
         SuperOff+=(SuperAtomMap.at(Idx[i])*3+Comp[i])*SuperStride;
         SubOff+=(SubAtomMap.at(Idx[i])*3+Comp[i])*SubStride;
      }
      for(const Atom& AtomI: Sys){//Unrolled loop
         //Second offsets
         size_t SuperOff2=SuperAtomMap.at(AtomI)*3;
         size_t SubOff2=SubAtomMap.at(AtomI)*3;
         for(size_t i=0;i<3;++i)//Actual filling
            Result[SuperOff+SuperOff2+i]+=Coeff*SubResult[SubOff+SubOff2+i];
      }
   }
   else{//Hessian and higher land here
      for(const Atom& AtomI : Sys){
         Idx.push_back(AtomI);
         for(size_t i=0;i<3;++i){
            Comp.push_back(i);
            FillDeriv(Result,SubResult,Coeff,Sys,
                      SuperAtomMap,SubAtomMap,Order,Idx,Comp);
            Comp.pop_back();
         }
         Idx.pop_back();
      }
   }
   
}

    void GetCoef(bool Even,const SN_t& NMer,const SNList_t& SNs,
                 map<string,double>& Coeffs){
        //If we don't have the NMer it's because we're assuming it's n-body
        //interaction is negligible, so don't follow that recursion
        size_t n=NMer.size()-1;
        if(SNs[n].count(NMer)==0)return;
        Coeffs[SNs[n].at(NMer)]+=(Even?1.0:-1.0);
        pulsar::math::PowerSetItr<SN_t> Frags(NMer,1,n);//no empty sets
        while(!Frags.Done()){
            GetCoef(!Even,*Frags,SNs,Coeffs);
            ++Frags;
        }
    }

}//End namespace