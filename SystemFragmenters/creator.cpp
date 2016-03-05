#include "Atomizer.hpp"

using bpmodule::modulemanager::ModuleCreationFuncs;


extern "C" {

ModuleCreationFuncs InsertSupermodule(void)
{
    ModuleCreationFuncs cf;
    cf.AddCppCreator<Atomizer>("Atomizer");
    return cf;
}



}

