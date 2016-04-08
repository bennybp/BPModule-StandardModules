#include "SystemFragmenters/Atomizer.hpp"
#include "SystemFragmenters/Bondizer.hpp"
#include "SystemFragmenters/Null.hpp"
#include "SystemFragmenters/UserDefined.hpp"


using bpmodule::modulemanager::ModuleCreationFuncs;


extern "C" {

ModuleCreationFuncs InsertSupermodule(void)
{
    ModuleCreationFuncs cf;
    cf.AddCppCreator<Atomizer>("Atomizer");
    cf.AddCppCreator<Bondizer>("Bondizer");
    cf.AddCppCreator<NullFragmenter>("NullFragmenter");
    cf.AddCppCreator<UserDefined>("UserDefined");
    return cf;
}



}

