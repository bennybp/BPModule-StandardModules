/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MBECommonHeader.hpp
 * Author: richard
 *
 * Created on April 4, 2016, 12:55 PM.  
 * 
 * \file Header files, typedefs, and using's that I have used in every file so
 * far.  Now in a common place.
 */

#ifndef MBECOMMONHEADER_HPP
#define MBECOMMONHEADER_HPP

#include <pulsar/system/System.hpp>
#include <pulsar/modulebase/ModuleBase.hpp>
#include <pulsar/modulebase/SystemFragmenter.hpp>
#include <pulsar/modulebase/EnergyMethod.hpp>
#include <pulsar/datastore/OptionMap.hpp>
#include <pulsar/datastore/Wavefunction.hpp>


using pulsar::system::System;
using pulsar::system::Atom;
using pulsar::system::SystemMap;
using pulsar::modulebase::SystemFragmenter;
using pulsar::modulebase::EnergyMethod;
using pulsar::datastore::OptionMap;
using pulsar::datastore::Wavefunction;

typedef pulsar::modulemanager::ModulePtr<SystemFragmenter> Fragmenter_t;
typedef pulsar::modulemanager::ModulePtr<EnergyMethod> EMethod_t;


#endif /* MBECOMMONHEADER_HPP */

