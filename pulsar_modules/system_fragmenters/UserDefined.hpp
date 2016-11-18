/* 
 * File:   UserDefined.hpp
 * Author: richard
 *
 * Created on March 31, 2016, 4:06 PM
 */

#pragma once

#include <pulsar/modulebase/SystemFragmenter.hpp>


/** \brief A fragmenter that relies on the user to tell it what the fragments
 *  are.
 * 
 *  For the time-being this class is pretty hacky.  For n fragments, it requires
 *  three options:
 *  FRAGMENT_NAMES an n element long list of strings each of which is the name
 *                 of the fragment
 *  ATOMS_PER_FRAG an n element long list of integers where each entry is the
 *                 number of atoms in each fragment
 *  FRAGMENTS      A list of atom indices such that the first set are for frag1
 *                 the second set for atom2, etc.
 * 
 *  This means that usage of this class is somewhat of an expert option as
 *  you have to be sure that this module gets the right system or else the
 *  fragmentation pattern is meaningless.
 */ 
class UserDefined : public pulsar::SystemFragmenter
{
public:
using SystemFragmenter::SystemFragmenter;

virtual pulsar::NMerSetType fragmentize_(const pulsar::System & mol);

};


