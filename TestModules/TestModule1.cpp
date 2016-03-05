#include <bpmodule/output/Output.hpp>
#include <bpmodule/modulemanager/ModuleManager.hpp>
#include "TestModule1.hpp"


using namespace bpmodule::output;
using namespace bpmodule::modulemanager;
using namespace bpmodule::exception;



TestModule1::TestModule1(unsigned long id)
    : Test_Base(id)
{
}



TestModule1::~TestModule1()
{
}



void TestModule1::RunTest_(void)
{
    Output("+++ In TestModule1: RunTest. Info: (%1%) %2% %3% v%4%\n", ID(), Key(), Name(), Version());

    Output("    Wavefunction: %1%\n", Wfn().UniqueString());
    Output("   Cache entries: %1%\n", Cache().Size());
    for(const auto & it : Cache().GetKeys())
        Output("                  > %1%\n", it);

    Output("   double_opt_def:    %1%\n", Options().Get<double>("double_opt_def"));
    Output("      int_opt_def:    %1%\n", Options().Get<int>("int_opt_def"));
    Output("     bool_opt_def:    %1%\n", Options().Get<bool>("bool_opt_def"));
    Output("      str_opt_def:    %1%\n", Options().Get<std::string>("str_opt_def"));
    Output("\n");

    if(Options().Has("double_opt"))
        Output("       double_opt:    %1%\n", Options().Get<double>("double_opt"));

    if(Options().Has("int_opt"))
        Output("          int_opt:    %1%\n", Options().Get<int>("int_opt"));

    if(Options().Has("bool_opt"))
        Output("         bool_opt:    %1%\n", Options().Get<bool>("bool_opt"));

    if(Options().Has("str_opt"))
        Output("          str_opt:    %1%\n", Options().Get<std::string>("str_opt"));

    Cache().Set( "Element 1", std::string("Something in the python cache") );
    Cache().Set( "Element 2", 42);
    Cache().Set( "Element 3", 42.0 );
    Cache().Set( "Element 4", std::vector<int>{ 1, 2, 3, 4} );
}



void TestModule1::CallRunTest_(const std::string & other)
{
    Output("+++ In TestModule1: CallRunTest with %1%\n", other);

    ModulePtr<Test_Base> tb2 = CreateChildModule<Test_Base>(other);
    Output("  + Obtained scoped module ID %1%\n", tb2->ID());
    tb2->RunTest();
    Output("  + Finished with scoped module %1%. Deleting automatically\n", tb2->ID());

    Output("+++Done\n");
}



void TestModule1::TestThrow_(void)
{
    Warning("+++ In TestModule1: Throwing an exception!\n");
    throw GeneralException("This is a test exception",
                           "Data1", "Hi",
                           "Data 2", "Hello");
}



void TestModule1::CallThrow_(const std::string & other)
{
    Output("+++ In TestModule1: CallThrowTest with %1%\n", other);

    ModulePtr<Test_Base> tb2 = CreateChildModule<Test_Base>(other);
    Output("  + Obtained scoped module ID %1%\n", tb2->ID());
    tb2->TestThrow();

    // shouldn't be called
    Output("+++Done\n");
}

