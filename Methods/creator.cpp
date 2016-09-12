
#include "Methods/MBE/MBE.hpp"
#include "Methods/SCF/Damping.hpp"
#include "Methods/SCF/DIIS.hpp"
#include "Methods/SCF/HFIterate.hpp"
#include "Methods/SCF/CoreGuess.hpp"
#include "Methods/SCF/BasicFockBuild.hpp"


using pulsar::modulemanager::ModuleCreationFuncs;

extern "C" {

ModuleCreationFuncs insert_supermodule(void){
    ModuleCreationFuncs cf;
    cf.add_cpp_creator<pulsarmethods::MBE>("MBE");
    cf.add_cpp_creator<pulsarmethods::Damping>("Damping");
    cf.add_cpp_creator<pulsarmethods::DIIS>("DIIS");
    cf.add_cpp_creator<pulsarmethods::HFIterate>("HFIterate");
    cf.add_cpp_creator<pulsarmethods::CoreGuess>("CoreGuess");
    cf.add_cpp_creator<pulsarmethods::BasicFockBuild>("BasicFockBuild");
    return cf;
}

}
