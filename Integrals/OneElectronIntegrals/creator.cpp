#include "Overlap.hpp"
#include "KineticEnergy.hpp"

using pulsar::modulemanager::ModuleCreationFuncs;


extern "C" {

ModuleCreationFuncs InsertSupermodule(void)
{
    ModuleCreationFuncs cf;
    cf.AddCppCreator<Overlap>("Overlap");
    cf.AddCppCreator<KineticEnergy>("KineticEnergy");
    return cf;
}



}

