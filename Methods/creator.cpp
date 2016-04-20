#include "Methods/BSSE/CP.hpp"
#include "Methods/BSSE/VMFC.hpp"
#include "Methods/MBE/MBE.hpp"
#include "Methods/MIM/MIM.hpp"


using pulsar::modulemanager::ModuleCreationFuncs;

extern "C" {

ModuleCreationFuncs InsertSupermodule(void){
    ModuleCreationFuncs cf;
    cf.AddCppCreator<pulsarmethods::CP>("CP");
    cf.AddCppCreator<pulsarmethods::MBE>("MBE");
    cf.AddCppCreator<pulsarmethods::MIM>("MIM");
    cf.AddCppCreator<pulsarmethods::VMFC>("VMFC");
    return cf;
}

}