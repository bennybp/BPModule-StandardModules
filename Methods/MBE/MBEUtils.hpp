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
#include<set>
#include<unordered_map>


namespace bpmodule{namespace system{
class System;
class Atom;
typedef std::map<std::string,System> SystemMap;
}}


namespace bpmethods{

///NMer's serial number (the indices of the fragments that make up the n-mer)
typedef std::set<size_t> SN_t;

///The resulting, binned n-mers, mapped by serial numbers
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

/*   \brief Function for projecting arbitrary order derivatives
 * 
 *   For the purposes of this function "super" refers to the
 *   system we are projecting into and "sub" is the system we are projecting from.
 *   At no point do we actually check if sub is truly a subset of super, nor do
 *   we rely on such a relationship.  This comes from the primary role of this
 *   function for projecting derivatives of subsystems onto those of the
 *   supersystem; however we also use it BSSE methods to project ghost atoms
 *   onto real atoms.
 *  
 *   Here's the plan.  For the Order-th derivative we have to expand a rank 
 *   "Order" tensor into another rank "Order" tensor.  Assume we have two maps: 
 *   SuperAtomMap, and SubAtomMap, which respectively tell us which number an
 *   atom is in the supersystem and the subsystem.  We assume that super
 *   derivative is in the order x, y, z for whatever atom in SuperAtomMap is 0,
 *    x,y,z for the atom in SuperAtomMap that is 1; etc. and similar for the
 *   sub derivative.  Technically the derivatives have symmetry,
 *   but we ignore that for now.  We do the actual filling by recursion.  
 *   At each level we loop over atoms
 *   and components and pass the result to the next level.  
 *   Once we have gone down Order levels, we have a fully
 *   specified index and simply evaluate the index.  We can avoid a 
 *   little-bit of recursion (and possibly pick up
 *   a tiny bit of vectorization if we unroll the last iteration.  
 *   We thus will have two offsets.  The first
 *   offset is the offset of the Order-1 indices.  The second offset is the 
 *   order-th index that we have unrolled.  Note that each order
 *   is sort of two orders because we get an offset from the Cartesian 
 *   components too.  This means that Basically our index looks like:
 *   \f[
 *   Idx=\left[
 *       AtomJ*3+CompJ+
 *       \sum_{i=0}^{Order-1} Atom[i]*(3*NAtoms)^{Order-1-i}+Comp[i]\right]
 *   \f]
 *   where the summation involving \f$Atom[i]\f$ is the first offset and 
 *   \f$AtomJ\f$ is the second offset.
 * 
 *    \param[in/out] Result The final result we are projecting into
 *    \param[in] SubResult The result we are projecting out of
 *    \param[in] Coeff The value our SubResult is getting scaled by
 *    \param[in] SuperAtomMap Map of Atom to index in the supersystem
 *    \param[in] SubAtomMap  Map of Atom to index in the subsystem
 *    \param[in] Order What rank derivative we are doing
 *    \param[in] Idx Buffer used for Atom index during recursion 
 *                    (you can ignore)
 *    \param[in] Comp Buffer used for Cartesian component during recursion
 *                    (you can ignore)
 */   
void FillDeriv(std::vector<double>& Result, 
               const std::vector<double>& SubResult,
               double Coeff,
               const bpmodule::system::System& Sys, 
               const std::unordered_map<bpmodule::system::Atom,size_t>& SuperAtomMap,
               const std::unordered_map<bpmodule::system::Atom,size_t>& SubAtomMap,
               size_t Order,
               std::vector<bpmodule::system::Atom> Idx=
                  std::vector<bpmodule::system::Atom>(),
               std::vector<size_t> Comp=std::vector<size_t>());

//Computes the MBE coefficients by recursion
void GetCoef(bool Even,const SN_t& NMer,
             const SNList_t& SNs,
             std::map<std::string,double>& Coeffs);

}//End namespace

#endif /* MBEUTILS_HPP */

