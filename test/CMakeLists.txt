cmake_minimum_required(VERSION 3.0)
include(CTest)

SET(CTEST_BINARY_DIRECTORY "${CMAKE_INSTALL_DIR}")
#Should find the Pulsar in the staging area
find_package(pulsar REQUIRED)

#Macro for defining a python test
function(pulsar_py_test test_name)
add_test(NAME ${test_name}_PY
         COMMAND ${PYTHON_EXECUTABLE} ${test_name}.py
         WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/test
)
install(FILES ${test_name}.py
        DESTINATION test
)
endfunction()

#Macro for defining a C++ test
function(pulsar_cxx_test test_name)
add_library(${test_name} MODULE ${test_name}.cpp)
target_include_directories(${test_name} PRIVATE ${PROJECT_SOURCE_DIR} pulsar)
target_link_libraries(${test_name} PRIVATE pulsar)
add_test(NAME ${test_name}_CXX
         COMMAND ${PYTHON_EXECUTABLE} TestCXX.py ${test_name}
         WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/test
)
install(TARGETS ${test_name}
        DESTINATION test
)
endfunction()

#Macro for defining both a Python and C++ test
function(pulsar_test test_name)
   pulsar_py_test(${test_name})
   pulsar_cxx_test(${test_name})
endfunction()

#This file will allow us to run ctest in the top-level build dir
#Basically it just defers to the actual top-level CTestTestfile.cmake in the
#build directory for this project
file(WRITE ${CMAKE_INSTALL_PREFIX}/CTestTestfile.cmake "subdirs(test)")

add_subdirectory(system_fragmenters)


configure_file("TestFxns.py.in" "TestFxns.py" @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/TestFxns.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/test)
install(FILES TestCXX.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/test)