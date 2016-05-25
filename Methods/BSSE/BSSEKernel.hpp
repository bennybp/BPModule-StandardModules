/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BSSEKernel.hpp
 * Author: richard
 *
 * Created on April 4, 2016, 6:11 PM
 * 
 * Functions that are used in several places for BSSE related tasks.
 * 
 * \todo Merge VMFC's kernel with SSFCKernel
 */

#ifndef BSSEKERNEL_HPP
#define BSSEKERNEL_HPP

#include <memory>
#include <unordered_map>
#include <map>
#include <string>
#include <pulsar/system/System.hpp>
#include <pulsar/system/Atom.hpp>
#include <pulsar/modulemanager/ModuleManager.hpp>
#include <pulsar/modulebase/EnergyMethod.hpp>
#include "Methods/MBE/MBEUtils.hpp"

namespace pulsarmethods{


///A struct that contains all the relevant data after ghosting a system
struct RealGhostData{
    ///Map of the Real atom to its corresponding ghost
    std::unordered_map<pulsar::system::Atom,pulsar::system::Atom> Real2Ghost;
    ///Map of any atom (real or ghost) to its index
    std::unordered_map<pulsar::system::Atom,size_t> Atom2Idx;
    ///Map of an atom (real or ghost) to its real index
    std::unordered_map<pulsar::system::Atom,size_t> Atom2RealIdx;
    ///A universe consisting of all the real and ghost atoms
    pulsar::system::AtomSetUniverse NewSystem;
    ///A system comprised of only the real atoms
    std::unique_ptr<pulsar::system::System> RealSystem;
};

///Returns an appropriately structured RealGhostData struct for Mol
RealGhostData GhostTheSystem(const pulsar::system::System& Mol);

///Returns a name consistent with our naming conventions
std::string MakeName(const SN_t& FullSN,const SN_t& SN);


/** \brief The core kernel for creating SSFC-like systems
 * 
 *  This function makes the systems needed for various SSFC-like system
 *  computations as well as figuring out their weights for insertion into
 *  MIM.  I made this pretty general so that someone could at some point
 *  use it to code up many-body counterpoise corrections as well as
 *  many-ghost-many-body corrections by just tweaking the input options.
 *  It is not clear to me that such methods are actually worth coding up
 *  given VMFC, which is why I have punted.
 *
 *  \param[in,out] NMers The set of systems that you need to run
 *  \param[in] Data The struct of info related to the ghosted system
    \param[in] MinRealOrder The minimum number of real monomers in a fragment
 *  \param[in] MaxRealOrder The maximum number of real monomers in a fragment
 *  \param[in] MinGhostOrder The minimum number of ghost monomers in a fragment
 *  \param[in] MaxGhostOrder The maximum number of ghost monomers in a fragment
 *  \param[in] AddReal Should we also include the real-only systems to NMers
 *                     (coefficients will also be added to return)
 *  \return The coefficients corresponding to each system
 */
std::map<std::string,double> SSFCKernel(pulsar::system::SystemMap& NMers,
                                        const RealGhostData& Data, 
                                        size_t MaxRealOrder,size_t MinRealOrder,
                                      size_t MaxGhostOrder,size_t MinGhostOrder
);


/** \brief The code for running any series of BSSE computations
 * 
 *   Basically this code is a wrapper around the creation of a UserDefined
 *   fragmenter, the appropriate setting of that fragmenter's options, and
 *   then an execution of MIM.  
 * 
 *  \param[in] AllFrags All of the real/ghosted systems to run
 *  \param[in] Coeffs The coefficients of the systems, matched by key
 *  \param[in] Data The struct of info related to the new, ghosted system
 *  \param[in] Order What order derivative we are computing
 *  \param[in] ID The parent module's ID
 *  \param[in] MM The module manager of the parent
 *  \param[in] MethodName The key of the method we are going to call
 *  \param[in] MIMName The key for MIM method you want us to call
 *  \return The derivative you were interested in
 */ 
pulsar::modulebase::EnergyMethod::DerivReturnType 
    RunCalcs(const pulsar::system::SystemMap& AllFrags,
        const pulsar::datastore::Wavefunction& Wfn,
                        const std::map<std::string,double>& Coeffs,
                        const RealGhostData& Data,
                        size_t Order,
                        ID_t ID,
                        pulsar::modulemanager::ModuleManager& MM,
                        const std::string& MethodName,
                        const std::string& MIMName);

}//End namespace

#endif /* BSSEKERNEL_HPP */

