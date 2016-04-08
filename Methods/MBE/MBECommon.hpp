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

#include <bpmodule/system/System.hpp>
#include <bpmodule/modulebase/ModuleBase.hpp>
#include <bpmodule/modulebase/SystemFragmenter.hpp>
#include <bpmodule/modulebase/EnergyMethod.hpp>
#include <bpmodule/datastore/OptionMap.hpp>


using bpmodule::system::System;
using bpmodule::system::Atom;
using bpmodule::system::SystemMap;
using bpmodule::modulebase::SystemFragmenter;
using bpmodule::modulebase::EnergyMethod;
using bpmodule::datastore::OptionMap;

typedef bpmodule::modulemanager::ModulePtr<SystemFragmenter> Fragmenter_t;
typedef bpmodule::modulemanager::ModulePtr<EnergyMethod> EMethod_t;


#endif /* MBECOMMONHEADER_HPP */

