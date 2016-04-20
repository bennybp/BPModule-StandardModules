#ifndef _GUARD_TESTEXTLIB_HPP_
#define _GUARD_TESTEXTLIB_HPP_

#include <pulsar/modulebase/Test_Base.hpp>

class TestExtLib : public pulsar::modulebase::Test_Base
{
public:
    TestExtLib(ID_t id);

    virtual void RunTest_(void);

    virtual void CallRunTest_(const std::string & other);

    virtual void CallRunTest2_(const std::string & other1, const std::string & other2);

    virtual void TestThrow_(void);

    virtual void CallThrow_(const std::string & other);

    virtual ~TestExtLib();

};


#endif
