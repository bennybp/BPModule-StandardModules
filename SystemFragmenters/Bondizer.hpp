/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Bondizer.hpp
 * Author: richard
 *
 * Created on March 15, 2016, 10:05 AM
 */

#ifndef BONDIZER_HPP
#define BONDIZER_HPP

#include <bpmodule/modulebase/SystemFragmenter.hpp>


/** \brief Breaks a system apart into groups of atoms that are seperated by at
 *         most a pre-defined number of bonds
 * 
 *   This fragmenter first calls bpmodule::system::GetConns() so I suggest 
 *   modifying the covalent radii of the atoms you want to be bonded/not bonded
 *   to accomplish this. For example, say the user wants two atoms to be bonded.
 *   If you compute
 *   the distance between them and then split that distance in the same ratio as
 *   the original radii this gives an effective covalent radius for each atom
 *   consistent with the specified connectivity.  Using the same ratio should
 *   minimally impact the other bonds the atoms are involved in.
 * 
 *   This class recognizes the following options:
 *     - NBONDS The maximum number of bonds that may exist between two atoms
 */ 
class Bondizer : public bpmodule::modulebase::SystemFragmenter
{
public:
    using bpmodule::modulebase::SystemFragmenter::SystemFragmenter;

    virtual bpmodule::system::SystemMap Fragmentize_(const bpmodule::system::System & mol);

};


#endif /* BONDIZER_HPP */

