#ifndef _GUARD_TESTMODULE1_HPP_
#define _GUARD_TESTMODULE1_HPP_

#include <pulsar/modulebase/Test_Base.hpp>

class TestModule1 : public pulsar::modulebase::Test_Base
{
public:
    TestModule1(ID_t id);

    virtual void RunTest_(void);

    virtual void CallRunTest_(const std::string & other);

    virtual void CallRunTest2_(const std::string & other1, const std::string & other2);

    virtual void TestThrow_(void);

    virtual void CallThrow_(const std::string & other);

    virtual ~TestModule1();

};


#endif
