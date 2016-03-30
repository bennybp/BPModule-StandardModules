#include "MBE/MBE.hpp"
#include "MIM/MIM.hpp"

using bpmodule::modulemanager::ModuleCreationFuncs;

extern "C" {

ModuleCreationFuncs InsertSupermodule(void){
    ModuleCreationFuncs cf;
    cf.AddCppCreator<bpmethods::MBE>("MBE");
    cf.AddCppCreator<bpmethods::MIM>("MIM");
    return cf;
}

}