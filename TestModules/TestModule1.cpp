#include <pulsar/output/Output.hpp>
#include <pulsar/modulemanager/ModuleManager.hpp>
#include "TestModule1.hpp"


using namespace pulsar;

struct SomeStruct
{
    int i;
};

struct SomeSerializableStruct
{
    int i;

    template<typename Archive>
    void serialize(Archive & ar)
    {
        ar(i);
    }
};

TestModule1::TestModule1(ID_t id)
    : Test_Base(id)
{
}



TestModule1::~TestModule1()
{
}



void TestModule1::run_test_(void)
{
    out.output("+++ In TestModule1: run_test. Info: (%?) %? %? v%?\n", id(), key(), name(), version());

    out.output("   double_opt_def:    %?\n", options().get<double>("double_opt_def"));
    out.output("      int_opt_def:    %?\n", options().get<int>("int_opt_def"));
    out.output("     bool_opt_def:    %?\n", options().get<bool>("bool_opt_def"));
    out.output("      str_opt_def:    %?\n", options().get<std::string>("str_opt_def"));
    out.output("\n");

    if(options().has("double_opt"))
        out.output("       double_opt:    %?\n", options().get<double>("double_opt"));

    if(options().has("int_opt"))
        out.output("          int_opt:    %?\n", options().get<int>("int_opt"));

    if(options().has("bool_opt"))
        out.output("         bool_opt:    %?\n", options().get<bool>("bool_opt"));

    if(options().has("str_opt"))
        out.output("          str_opt:    %?\n", options().get<std::string>("str_opt"));
}



void TestModule1::call_run_test_(const std::string & other)
{
    out.output("+++ In TestModule1: call_run_test with %?\n", other);

    ModulePtr<Test_Base> tb2 = create_child<Test_Base>(other);
    out.output("  + Obtained scoped module ID %?\n", tb2->id());
    tb2->run_test();
    out.output("  + Finished with scoped module %?. Deleting automatically\n", tb2->id());

    out.output("+++Done\n");
}


void TestModule1::call_run_test2_(const std::string & other1, const std::string & other2)
{
    out.output("+++ In TestModule1: call_run_test with %? %?\n", other1, other2);

    ModulePtr<Test_Base> tb2 = create_child<Test_Base>(other1);
    out.output("  + Obtained scoped module ID %?\n", tb2->id());
    tb2->call_run_test(other2);
    out.output("  + Finished with scoped module %?. Deleting automatically\n", tb2->id());

    out.output("+++Done\n");
}



void TestModule1::test_throw_(void)
{
    out.warning("+++ In TestModule1: Throwing an exception!\n");
    throw PulsarException("This is a test exception",
                           "Data1", "Hi",
                           "Data 2", "Hello");
}



void TestModule1::call_throw_(const std::string & other)
{
    out.output("+++ In TestModule1: call_throwTest with %?\n", other);

    ModulePtr<Test_Base> tb2 = create_child<Test_Base>(other);
    out.output("  + Obtained scoped module ID %?\n", tb2->id());
    tb2->test_throw();

    // shouldn't be called
    out.output("+++Done\n");
}


void TestModule1::add_to_cache_(const std::string & key, unsigned int policy)
{
    std::string val("TestModule1::DataInCacheForKey:");
    val += key;
    cache().set(key, val, policy);
}


void TestModule1::get_from_cache_(const std::string & key)
{
    std::string val = *(cache().get<std::string>(key));
    out.output("From cache: Key = %? Value = %?\n", key, val);
}


