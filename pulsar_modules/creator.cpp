#include "pulsar_modules/system_fragmenters/Atomizer.hpp"
#include "pulsar_modules/system_fragmenters/Bondizer.hpp"
#include "pulsar_modules/system_fragmenters/CrystalFragger.hpp"
#include "pulsar_modules/system_fragmenters/CPGhoster.hpp"
#include "pulsar_modules/system_fragmenters/NMerizer.hpp"
#include "pulsar_modules/system_fragmenters/UserDefined.hpp"
#include "pulsar_modules/system_fragmenters/VMFCGhoster.hpp"
#include "pulsar_modules/methods/mbe/MBE.hpp"
//#include "methods/scf/Damping.hpp"
//#include "methods/scf/DIIS.hpp"
//#include "methods/scf/HFIterate.hpp"
//#include "methods/scf/CoreGuess.hpp"
//#include "methods/scf/BasicFockBuild.hpp"


using pulsar::ModuleCreationFuncs;

extern "C" {

ModuleCreationFuncs insert_supermodule(void){
    ModuleCreationFuncs cf;
    cf.add_cpp_creator<MBE>("MBE");
//    cf.add_cpp_creator<pulsarmethods::Damping>("Damping");
//    cf.add_cpp_creator<pulsarmethods::DIIS>("DIIS");
//    cf.add_cpp_creator<pulsarmethods::HFIterate>("HFIterate");
//    cf.add_cpp_creator<pulsarmethods::CoreGuess>("CoreGuess");
//    cf.add_cpp_creator<pulsarmethods::BasicFockBuild>("BasicFockBuild");
    cf.add_cpp_creator<Atomizer>("Atomizer");
    cf.add_cpp_creator<Bondizer>("Bondizer");
    cf.add_cpp_creator<CrystalFragger>("CrystalFragger");
    cf.add_cpp_creator<CPGhoster>("CPGhoster");
    cf.add_cpp_creator<VMFCGhoster>("VMFCGhoster");
    cf.add_cpp_creator<UserDefined>("UserDefined");
    cf.add_cpp_creator<NMerizer>("NMerizer");
    return cf;
}

}
