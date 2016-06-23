#include "SystemFragmenters/Atomizer.hpp"
#include "SystemFragmenters/Bondizer.hpp"
#include "SystemFragmenters/CrystalFragger.hpp"
#include "SystemFragmenters/Ghoster.hpp"
#include "SystemFragmenters/Null.hpp"
#include "SystemFragmenters/UserDefined.hpp"

using pulsar::modulemanager::ModuleCreationFuncs;

extern "C" {

ModuleCreationFuncs insert_supermodule(void)
{
    ModuleCreationFuncs cf;
    cf.add_cpp_creator<Atomizer>("Atomizer");
    cf.add_cpp_creator<Bondizer>("Bondizer");
    cf.add_cpp_creator<CrystalFragger>("CrystalFragger");
    cf.add_cpp_creator<Ghoster>("Ghoster");
    cf.add_cpp_creator<CPGhoster>("CPGhoster");
    cf.add_cpp_creator<VMFCGhoster>("VMFCGhoster");
    cf.add_cpp_creator<NullFragmenter>("NullFragmenter");
    cf.add_cpp_creator<UserDefined>("UserDefined");
    return cf;
}



}

