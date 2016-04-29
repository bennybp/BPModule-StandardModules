#include "ReferenceERI.hpp"

using pulsar::modulemanager::ModuleCreationFuncs;


extern "C" {

ModuleCreationFuncs InsertSupermodule(void)
{
    ModuleCreationFuncs cf;
    cf.AddCppCreator<ReferenceERI>("ReferenceERI");
    return cf;
}



}

