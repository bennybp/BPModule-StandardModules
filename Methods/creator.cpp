#include "Methods/BSSE/CP.hpp"
#include "Methods/BSSE/VMFC.hpp"
#include "Methods/MBE/MBE.hpp"
#include "Methods/MIM/MIM.hpp"
#include "Methods/SCF/Damping.hpp"
#include "Methods/SCF/DIIS.hpp"
#include "Methods/SCF/HFIterate.hpp"
#include "Methods/SCF/CoreGuess.hpp"
#include "Methods/SCF/BasicFockBuild.hpp"


using pulsar::modulemanager::ModuleCreationFuncs;

extern "C" {

ModuleCreationFuncs InsertSupermodule(void){
    ModuleCreationFuncs cf;
    cf.AddCppCreator<pulsarmethods::CP>("CP");
    cf.AddCppCreator<pulsarmethods::MBE>("MBE");
    cf.AddCppCreator<pulsarmethods::MIM>("MIM");
    cf.AddCppCreator<pulsarmethods::VMFC>("VMFC");
    cf.AddCppCreator<pulsarmethods::Damping>("Damping");
    cf.AddCppCreator<pulsarmethods::DIIS>("DIIS");
    cf.AddCppCreator<pulsarmethods::HFIterate>("HFIterate");
    cf.AddCppCreator<pulsarmethods::CoreGuess>("CoreGuess");
    cf.AddCppCreator<pulsarmethods::BasicFockBuild>("BasicFockBuild");
    return cf;
}

}
