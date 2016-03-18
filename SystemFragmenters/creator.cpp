#include "SystemFragmenters/Atomizer.hpp"
#include "SystemFragmenters/Bondizer.hpp"

using bpmodule::modulemanager::ModuleCreationFuncs;


extern "C" {

ModuleCreationFuncs InsertSupermodule(void)
{
    ModuleCreationFuncs cf;
    cf.AddCppCreator<Atomizer>("Atomizer");
    cf.AddCppCreator<Bondizer>("Bondizer");
    return cf;
}



}

