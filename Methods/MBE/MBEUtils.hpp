/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MBEUtils.hpp
 * Author: richard
 *
 * Created on March 31, 2016, 5:01 PM
 */

#ifndef MBEUTILS_HPP
#define MBEUTILS_HPP
#include<string>
#include<map>
#include<vector>

namespace bpmodule{namespace system{
class System;    
typedef std::map<std::string,System> SystemMap;
}}


namespace bpmethods{
///NMer's serial number (the indices of the fragments that make up the n-mer)
typedef std::set<size_t> SN_t;
///The resulting, binned n-mers
typedef std::vector<std::map<SN_t,std::string>> SNList_t;


/** \brief Bins a set of NMers by size
 * 
 *   In fragment based methods we often would like to have the n-mers
 *   sorted by size.  Furthermore, we often want to do manipulations on
 *   their indices.  This function takes the object that is the current
 *   output of a SystemFragmenter and parses the string keys so that
 *   you get back a list such that the \f$n\f$-th element is a map
 *   between \f$n\f$-mer serial numbers and their corresponding 
 *   key in the input map.  It assumes that keys are of the form: "1_2_3"
 *   for the trimer made by the union of the systems with keys "1", "2", and "3"
 */ 
SNList_t BinNMers(const bpmodule::system::SystemMap& Frags);


}//End namespace

#endif /* MBEUTILS_HPP */

