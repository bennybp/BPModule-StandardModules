#ifndef _GUARD_TESTEXTLIB_HPP_
#define _GUARD_TESTEXTLIB_HPP_

#include <pulsar/modulebase/Test_Base.hpp>

class TestExtLib : public pulsar::Test_Base
{
public:
    TestExtLib(ID_t id);

    virtual void run_test_(void);

    virtual void call_run_test_(const std::string & other);

    virtual void call_run_test2_(const std::string & other1, const std::string & other2);

    virtual void test_throw_(void);

    virtual void call_throw_(const std::string & other);

    virtual void add_to_cache_(const std::string & key, unsigned int policy);

    virtual void get_from_cache_(const std::string & key);

    virtual ~TestExtLib();

};


#endif
