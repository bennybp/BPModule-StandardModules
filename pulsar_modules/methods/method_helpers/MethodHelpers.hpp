/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/** \file MethodHelpers.hpp
 *  \brief These are helper functions for doing common tasks
 *  \author Ryan M. Richard
 *  \version 1.0
 *  \date June 2, 2016
 */

#ifndef PULSAR_GUARD_METHODHELPERS_HPP
#define PULSAR_GUARD_METHODHELPERS_HPP

#include <pulsar/modulebase/EnergyMethod.hpp>//For DerivReturnType
#include <pulsar/datastore/Wavefunction.hpp>
#include <pulsar/types.h>


/** \brief Given a series of computations to perform runs them in embarassingly
 *         parallel fashion and combines the results
 * 
 *  Many methods can be thought of as running a linear combination of methods.
 *  This function takes a series of EnergyMethod keys \p Keys, a series of
 *  wavefunctions \p Wfns, and a series of linear combination coeficents \p Cs
 *  and proceeds to run the requested computations in an
 *  embarassingly parallel fashion.  If you only provide one wavefunction we
 *  will assume that we are using the same wavefunction for each computation
 *  and it is the keys that are changing.  The alternative is also possible,
 *  if you provide a single key then we will assume the wavefunctions are
 *  changing.  The number of coefficients will determine the number of tasks
 *  that are to be run and it must be equal to either the number of keys or
 *  the number of wavefunctions (or both).
 * 
 *  \param[in] Keys The key or keys we are to use
 *  \param[in] Wfns The Wfns to run
 *  \param[in] Cs   The weights of each task in the final expansion
 *  \param[in] Deriv What order derivative are we running?
 *  \return A vector of the requested derivatives and their wavefunctions
 */
std::vector<pulsar::DerivReturnType> RunSeriesOfMethods(
    pulsar::ModuleManager& MM,
    ID_t ID,
    const std::vector<std::string>& Keys,
    const std::vector<pulsar::Wavefunction>& Wfns, 
    //const std::vector<double> Cs,
    size_t Deriv
);

/** \brief Maps one atomic tensor quantity to another
 *
 *   \param[in] Result The final result
 *   \param[in] SubResult Some result to map to the Result
 *   \param[in] C The coefficient to scale 
 */
void FillDeriv(std::vector<double>& Result, 
               const std::vector<double>& SubResult,
               double C,
               const pulsar::System& Sys, 
               const std::unordered_map<pulsar::Atom,size_t>& SprAtomMap,
               const std::unordered_map<pulsar::Atom,size_t>& SubAtomMap,
               size_t Order,
               std::vector<pulsar::Atom> Idx=
                  std::vector<pulsar::Atom>(),
               std::vector<size_t> Comp=std::vector<size_t>());

#endif /* PULSAR_GHUARD_METHODHELPERS_HPP */

