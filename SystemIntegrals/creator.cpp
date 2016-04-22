#include "NuclearRepulsion.hpp"

using pulsar::modulemanager::ModuleCreationFuncs;


extern "C" {

ModuleCreationFuncs InsertSupermodule(void)
{
    ModuleCreationFuncs cf;
    cf.AddCppCreator<NuclearRepulsion>("NuclearRepulsion");
    return cf;
}



}

