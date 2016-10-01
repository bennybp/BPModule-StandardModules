#include <pulsar/output/Output.hpp>
#include <pulsar/modulemanager/ModuleManager.hpp>

#include "TestExtLib.hpp"

// From this module
#include "staticlib.hpp"
#include "dynlib.hpp"

using namespace pulsar::output;
using namespace pulsar::modulemanager;
using namespace pulsar::exception;



TestExtLib::TestExtLib(ID_t id)
               : Test_Base(id)
{
}



TestExtLib::~TestExtLib()
{
}



void TestExtLib::run_test_(void)
{
    out.output("+++ In TestExtLib: run_test. Info: (%?) %? %? v%?\n", id(), key(), name(), version());
    out.output("*** Testing calling other libraries\n");
    std::string sstr = Static_GetString();
    std::string dstr = Dynamic_GetString();

    out.output("       From static: %?\n", sstr);
    out.output("      From dynamic: %?\n", dstr);
}



void TestExtLib::call_run_test_(const std::string & other)
{
}

void TestExtLib::call_run_test2_(const std::string & other1, const std::string & other2)
{
}



void TestExtLib::test_throw_(void)
{
}



void TestExtLib::call_throw_(const std::string & other)
{
}

void TestExtLib::add_to_cache_(const std::string & key, unsigned int policy)
{
}

void TestExtLib::get_from_cache_(const std::string & key)
{
}

