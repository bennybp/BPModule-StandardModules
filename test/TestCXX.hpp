/*! \file
 *
 * \brief  Contains the macro definition of a dummy class for testing Pulsar
 *         C++ side
 */

#pragma once

#include<pulsar/testing/Tester.hpp>
#include<pulsar/modulebase/Test_Base.hpp>

using pulsar::Tester;
using pulsar::Test_CXX_Base;
using pulsar::ModuleCreationFuncs;
#define TEST_CLASS(test_name)\
    class test_name : public Test_CXX_Base {\
protected:\
        virtual void run_test_(void);\
public:\
    using Test_CXX_Base::Test_CXX_Base;\
};\
extern "C" {\
ModuleCreationFuncs insert_supermodule(void){\
    ModuleCreationFuncs cf;\
    cf.add_cpp_creator<test_name>(#test_name);\
    return cf;\
}}\
void test_name::run_test_()

