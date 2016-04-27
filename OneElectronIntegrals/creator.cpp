#include "Overlap.hpp"

using pulsar::modulemanager::ModuleCreationFuncs;


extern "C" {

ModuleCreationFuncs InsertSupermodule(void)
{
    ModuleCreationFuncs cf;
    cf.AddCppCreator<Overlap>("Overlap");
    return cf;
}



}

