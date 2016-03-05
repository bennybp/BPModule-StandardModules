#ifndef _GUARD_TESTMODULE1_HPP_
#define _GUARD_TESTMODULE1_HPP_

#include <bpmodule/modulebase/Test_Base.hpp>

class TestModule1 : public bpmodule::modulebase::Test_Base
{
public:
    TestModule1(unsigned long id);

    virtual void RunTest_(void);

    virtual void CallRunTest_(const std::string & other);

    virtual void TestThrow_(void);

    virtual void CallThrow_(const std::string & other);

    virtual ~TestModule1();

};


#endif
