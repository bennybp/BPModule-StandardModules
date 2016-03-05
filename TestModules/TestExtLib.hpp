#ifndef _GUARD_TESTEXTLIB_HPP_
#define _GUARD_TESTEXTLIB_HPP_

#include <bpmodule/modulebase/Test_Base.hpp>

class TestExtLib : public bpmodule::modulebase::Test_Base
{
public:
    TestExtLib(unsigned long id);

    virtual void RunTest_(void);

    virtual void CallRunTest_(const std::string & other);

    virtual void TestThrow_(void);

    virtual void CallThrow_(const std::string & other);

    virtual ~TestExtLib();

};


#endif
