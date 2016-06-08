#include "SystemFragmenters/Atomizer.hpp"
#include "SystemFragmenters/Bondizer.hpp"
#include "SystemFragmenters/CrystalFragger.hpp"
#include "SystemFragmenters/Ghoster.hpp"
#include "SystemFragmenters/Null.hpp"
#include "SystemFragmenters/UserDefined.hpp"

using pulsar::modulemanager::ModuleCreationFuncs;

extern "C" {

ModuleCreationFuncs InsertSupermodule(void)
{
    ModuleCreationFuncs cf;
    cf.AddCppCreator<Atomizer>("Atomizer");
    cf.AddCppCreator<Bondizer>("Bondizer");
    cf.AddCppCreator<CrystalFragger>("CrystalFragger");
    cf.AddCppCreator<Ghoster>("Ghoster");
    cf.AddCppCreator<CPGhoster>("CPGhoster");
    cf.AddCppCreator<VMFCGhoster>("VMFCGhoster");
    cf.AddCppCreator<NullFragmenter>("NullFragmenter");
    cf.AddCppCreator<UserDefined>("UserDefined");
    return cf;
}



}

