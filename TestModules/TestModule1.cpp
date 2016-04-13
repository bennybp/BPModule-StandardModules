#include <bpmodule/output/Output.hpp>
#include <bpmodule/modulemanager/ModuleManager.hpp>
#include "TestModule1.hpp"


using namespace bpmodule::modulemanager;
using namespace bpmodule::exception;



TestModule1::TestModule1(ID_t id)
    : Test_Base(id)
{
}



TestModule1::~TestModule1()
{
}



void TestModule1::RunTest_(void)
{
    out.Output("+++ In TestModule1: RunTest. Info: (%?) %? %? v%?\n", ID(), Key(), Name(), Version());

    out.Output("    Wavefunction: %?\n", InitialWfn().MyHash().String());
    out.Output("   Cache entries: %?\n", Cache().Size());
    for(const auto & it : Cache().GetKeys())
        out.Output("                  > %?\n", it);

    out.Output("   double_opt_def:    %?\n", Options().Get<double>("double_opt_def"));
    out.Output("      int_opt_def:    %?\n", Options().Get<int>("int_opt_def"));
    out.Output("     bool_opt_def:    %?\n", Options().Get<bool>("bool_opt_def"));
    out.Output("      str_opt_def:    %?\n", Options().Get<std::string>("str_opt_def"));
    out.Output("\n");

    if(Options().Has("double_opt"))
        out.Output("       double_opt:    %?\n", Options().Get<double>("double_opt"));

    if(Options().Has("int_opt"))
        out.Output("          int_opt:    %?\n", Options().Get<int>("int_opt"));

    if(Options().Has("bool_opt"))
        out.Output("         bool_opt:    %?\n", Options().Get<bool>("bool_opt"));

    if(Options().Has("str_opt"))
        out.Output("          str_opt:    %?\n", Options().Get<std::string>("str_opt"));

    Cache().Set( "Element 1", std::string("Something in the python cache") );
    Cache().Set( "Element 2", 42);
    Cache().Set( "Element 3", 42.0 );
    Cache().Set( "Element 4", std::vector<int>{ 1, 2, 3, 4} );
}



void TestModule1::CallRunTest_(const std::string & other)
{
    out.Output("+++ In TestModule1: CallRunTest with %?\n", other);

    ModulePtr<Test_Base> tb2 = CreateChildModule<Test_Base>(other);
    out.Output("  + Obtained scoped module ID %?\n", tb2->ID());
    tb2->RunTest();
    out.Output("  + Finished with scoped module %?. Deleting automatically\n", tb2->ID());

    out.Output("+++Done\n");
}


void TestModule1::CallRunTest2_(const std::string & other1, const std::string & other2)
{
    out.Output("+++ In TestModule1: CallRunTest with %? %?\n", other1, other2);

    ModulePtr<Test_Base> tb2 = CreateChildModule<Test_Base>(other1);
    out.Output("  + Obtained scoped module ID %?\n", tb2->ID());
    tb2->CallRunTest(other2);
    out.Output("  + Finished with scoped module %?. Deleting automatically\n", tb2->ID());

    out.Output("+++Done\n");
}



void TestModule1::TestThrow_(void)
{
    out.Warning("+++ In TestModule1: Throwing an exception!\n");
    throw GeneralException("This is a test exception",
                           "Data1", "Hi",
                           "Data 2", "Hello");
}



void TestModule1::CallThrow_(const std::string & other)
{
    out.Output("+++ In TestModule1: CallThrowTest with %?\n", other);

    ModulePtr<Test_Base> tb2 = CreateChildModule<Test_Base>(other);
    out.Output("  + Obtained scoped module ID %?\n", tb2->ID());
    tb2->TestThrow();

    // shouldn't be called
    out.Output("+++Done\n");
}

