#include "Integrals/Overlap.hpp"
#include "Integrals/Dipole.hpp"
#include "Integrals/KineticEnergy.hpp"
#include "Integrals/CoreBuild.hpp"
#include "Integrals/OneElectronPotential.hpp"
#include "Integrals/OneElectronProperty.hpp"
#include "Integrals/ReferenceERI.hpp"
#include "Integrals/EigenCacher.hpp"
#include "Integrals/NuclearRepulsion.hpp"
#include "Integrals/NuclearDipole.hpp"

using pulsar::modulemanager::ModuleCreationFuncs;


extern "C" {

ModuleCreationFuncs InsertSupermodule(void)
{
    ModuleCreationFuncs cf;
    cf.AddCppCreator<ReferenceERI>("ReferenceERI");
    cf.AddCppCreator<Overlap>("Overlap");
    cf.AddCppCreator<Dipole>("Dipole");
    cf.AddCppCreator<KineticEnergy>("KineticEnergy");
    cf.AddCppCreator<OneElectronPotential>("OneElectronPotential");
    cf.AddCppCreator<CoreBuild>("CoreBuild");
    cf.AddCppCreator<OneElectronProperty>("OneElectronProperty");
    cf.AddCppCreator<EigenCacher>("EigenCacher");
    cf.AddCppCreator<NuclearRepulsion>("NuclearRepulsion");
    cf.AddCppCreator<NuclearDipole>("NuclearDipole");
    return cf;
}



}

