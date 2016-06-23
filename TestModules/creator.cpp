#include "TestModule1.hpp"
#include "TestExtLib.hpp"

using pulsar::modulemanager::ModuleCreationFuncs;


extern "C" {

ModuleCreationFuncs insert_supermodule(void)
{
    ModuleCreationFuncs cf;
    cf.add_cpp_creator<TestModule1>("TestModule1");
    cf.add_cpp_creator<TestExtLib>("TestExtLib");
    return cf;
}



}

