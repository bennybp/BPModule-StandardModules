#include "TestModule1.hpp"
#include "TestExtLib.hpp"

using bpmodule::modulemanager::ModuleCreationFuncs;


extern "C" {

ModuleCreationFuncs InsertSupermodule(void)
{
    ModuleCreationFuncs cf;
    cf.AddCppCreator<TestModule1>("TestModule1");
    cf.AddCppCreator<TestExtLib>("TestExtLib");
    return cf;
}



}

