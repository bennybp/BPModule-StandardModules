#include "Methods/BSSE/CP.hpp"
#include "Methods/BSSE/VMFC.hpp"
#include "Methods/MBE/MBE.hpp"
#include "Methods/MIM/MIM.hpp"


using bpmodule::modulemanager::ModuleCreationFuncs;

extern "C" {

ModuleCreationFuncs InsertSupermodule(void){
    ModuleCreationFuncs cf;
    cf.AddCppCreator<bpmethods::CP>("CP");
    cf.AddCppCreator<bpmethods::MBE>("MBE");
    cf.AddCppCreator<bpmethods::MIM>("MIM");
    cf.AddCppCreator<bpmethods::VMFC>("VMFC");
    return cf;
}

}