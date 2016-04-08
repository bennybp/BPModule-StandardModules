#include <bpmodule/output/Output.hpp>
#include <bpmodule/modulemanager/ModuleManager.hpp>

#include "TestExtLib.hpp"

// From this module
#include "staticlib.hpp"
#include "dynlib.hpp"

using namespace bpmodule::output;
using namespace bpmodule::modulemanager;
using namespace bpmodule::exception;



TestExtLib::TestExtLib(ID_t id)
               : Test_Base(id)
{
}



TestExtLib::~TestExtLib()
{
}



void TestExtLib::RunTest_(void)
{
    out.Output("+++ In TestExtLib: RunTest. Info: (%?) %? %? v%?\n", ID(), Key(), Name(), Version());
    out.Output("*** Testing calling other libraries\n");
    std::string sstr = Static_GetString();
    std::string dstr = Dynamic_GetString();

    out.Output("       From static: %?\n", sstr);
    out.Output("      From dynamic: %?\n", dstr);
}



void TestExtLib::CallRunTest_(const std::string & other)
{
}

void TestExtLib::CallRunTest2_(const std::string & other1, const std::string & other2)
{
}



void TestExtLib::TestThrow_(void)
{
}



void TestExtLib::CallThrow_(const std::string & other)
{
}

