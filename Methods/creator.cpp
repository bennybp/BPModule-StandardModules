#include "MBE/MBE.hpp"
#include "MIM/MIM.hpp"
#include "VMFC/VMFC.hpp"

using bpmodule::modulemanager::ModuleCreationFuncs;

extern "C" {

ModuleCreationFuncs InsertSupermodule(void){
    ModuleCreationFuncs cf;
    cf.AddCppCreator<bpmethods::MBE>("MBE");
    cf.AddCppCreator<bpmethods::MIM>("MIM");
    cf.AddCppCreator<bpmethods::VMFC>("VMFC");
    return cf;
}

}