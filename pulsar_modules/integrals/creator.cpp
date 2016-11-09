#include "Integrals/OSOverlap.hpp"
#include "Integrals/OSDipole.hpp"
#include "Integrals/OSKineticEnergy.hpp"
#include "Integrals/OSOneElectronPotential.hpp"

#include "Integrals/OneElectronIntegralSum.hpp"
#include "Integrals/OneElectronProperty.hpp"
#include "Integrals/ReferenceERI.hpp"
#include "Integrals/OneElectron_Eigen.hpp"
#include "Integrals/NuclearRepulsion.hpp"
#include "Integrals/NuclearDipole.hpp"

using pulsar::modulemanager::ModuleCreationFuncs;
using namespace psr_modules::integrals;

extern "C" {

ModuleCreationFuncs insert_supermodule(void)
{
    ModuleCreationFuncs cf;
    cf.add_cpp_creator<ReferenceERI>("ReferenceERI");
    cf.add_cpp_creator<OSOverlap>("OSOverlap");
    cf.add_cpp_creator<OSDipole>("OSDipole");
    cf.add_cpp_creator<OSKineticEnergy>("OSKineticEnergy");
    cf.add_cpp_creator<OSOneElectronPotential>("OSOneElectronPotential");
    cf.add_cpp_creator<OneElectronIntegralSum>("OneElectronIntegralSum");
    cf.add_cpp_creator<OneElectronProperty>("OneElectronProperty");
    cf.add_cpp_creator<OneElectron_Eigen>("OneElectron_Eigen");
    cf.add_cpp_creator<NuclearRepulsion>("NuclearRepulsion");
    cf.add_cpp_creator<NuclearDipole>("NuclearDipole");
    return cf;
}



}

