#include "NuclearRepulsion.hpp"
#include "NuclearDipole.hpp"

using pulsar::modulemanager::ModuleCreationFuncs;


extern "C" {

ModuleCreationFuncs InsertSupermodule(void)
{
    ModuleCreationFuncs cf;
    cf.AddCppCreator<NuclearRepulsion>("NuclearRepulsion");
    cf.AddCppCreator<NuclearDipole>("NuclearDipole");
    return cf;
}



}

