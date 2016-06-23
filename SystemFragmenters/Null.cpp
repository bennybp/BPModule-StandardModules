/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "SystemFragmenters/Null.hpp"
using namespace pulsar::modulebase;

NMerSetType NullFragmenter::fragmentize_(const pulsar::system::System & mol){
    NMerInfo NMer;
    NMer.SN.insert("SYSTEM");
    NMer.NMer=mol;
    return {{NMer.SN,NMer}};
}

