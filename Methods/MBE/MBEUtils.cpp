/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <sstream>
#include <bpmodule/system/System.hpp>
#include "Methods/MBE/MBEUtils.hpp"

using std::stringstream;
using std::string;
using std::vector;
using bpmodule::system::SystemMap;



namespace bpmethods{

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
}//End namespace