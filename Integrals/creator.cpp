#include "Integrals/Overlap.hpp"
#include "Integrals/Dipole.hpp"
#include "Integrals/KineticEnergy.hpp"
#include "Integrals/CoreBuild.hpp"
#include "Integrals/OneElectronPotential.hpp"
#include "Integrals/OneElectronProperty.hpp"
#include "Integrals/ReferenceERI.hpp"
#include "Integrals/OneElectron_Eigen.hpp"
#include "Integrals/NuclearRepulsion.hpp"
#include "Integrals/NuclearDipole.hpp"

using pulsar::modulemanager::ModuleCreationFuncs;


extern "C" {

ModuleCreationFuncs insert_supermodule(void)
{
    ModuleCreationFuncs cf;
    cf.add_cpp_creator<ReferenceERI>("ReferenceERI");
    cf.add_cpp_creator<Overlap>("Overlap");
    cf.add_cpp_creator<Dipole>("Dipole");
    cf.add_cpp_creator<KineticEnergy>("KineticEnergy");
    cf.add_cpp_creator<OneElectronPotential>("OneElectronPotential");
    cf.add_cpp_creator<CoreBuild>("CoreBuild");
    cf.add_cpp_creator<OneElectronProperty>("OneElectronProperty");
    cf.add_cpp_creator<OneElectron_Eigen>("OneElectron_Eigen");
    cf.add_cpp_creator<NuclearRepulsion>("NuclearRepulsion");
    cf.add_cpp_creator<NuclearDipole>("NuclearDipole");
    return cf;
}



}

