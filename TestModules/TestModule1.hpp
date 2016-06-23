#ifndef _GUARD_TESTMODULE1_HPP_
#define _GUARD_TESTMODULE1_HPP_

#include <pulsar/modulebase/Test_Base.hpp>

class TestModule1 : public pulsar::modulebase::Test_Base
{
public:
    TestModule1(ID_t id);

    virtual void run_test_(void);

    virtual void call_run_test_(const std::string & other);

    virtual void call_run_test2_(const std::string & other1, const std::string & other2);

    virtual void test_throw_(void);

    virtual void call_throw_(const std::string & other);

    virtual ~TestModule1();

};


#endif
