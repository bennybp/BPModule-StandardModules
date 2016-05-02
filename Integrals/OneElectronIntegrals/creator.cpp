#include "Overlap.hpp"
#include "KineticEnergy.hpp"
#include "OneElectronPotential.hpp"

using pulsar::modulemanager::ModuleCreationFuncs;


extern "C" {

ModuleCreationFuncs InsertSupermodule(void)
{
    ModuleCreationFuncs cf;
    cf.AddCppCreator<Overlap>("Overlap");
    cf.AddCppCreator<KineticEnergy>("KineticEnergy");
    cf.AddCppCreator<OneElectronPotential>("OneElectronPotential");
    return cf;
}



}

