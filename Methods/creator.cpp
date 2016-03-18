#include "MBE/MBE.hpp"

using bpmodule::modulemanager::ModuleCreationFuncs;

extern "C" {

ModuleCreationFuncs InsertSupermodule(void){
    ModuleCreationFuncs cf;
    cf.AddCppCreator<bpmethods::MBE>("MBE");
    return cf;
}

}