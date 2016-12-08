cmake_minimum_required(VERSION 3.0)
include(CTest)

SET(CTEST_BINARY_DIRECTORY "${CMAKE_INSTALL_DIR}")

#Should find the Pulsar in the staging area
#This sets PULSAR_RUNTEST variable
find_package(pulsar REQUIRED)


#Macro for running a test via Pulsar's RunTest script
function(pulsar_runtest test_name test_path)
  add_test(NAME ${test_name}
        COMMAND ${PYTHON_EXECUTABLE} ${PULSAR_RUNTEST} ${test_path} ${STAGE_DIR}/modules
        WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}
        )

endfunction()

#Macro for defining a python test
function(pulsar_py_test dir test_name)
  install(FILES ${test_name}.py DESTINATION ${dir})
  set(${test_name}_PY_PATH ${CMAKE_INSTALL_PREFIX}/${dir}/${test_name}.py)
  pulsar_runtest(${test_name}_PY ${CMAKE_INSTALL_PREFIX}/${dir}/${test_name}.py)
endfunction()


#Macro for defining a C++ test
function(pulsar_cxx_test dir test_name)
  add_library(${test_name} MODULE ${test_name}.cpp)
  target_include_directories(${test_name} PRIVATE ${PROJECT_SOURCE_DIR})
  target_link_libraries(${test_name} PRIVATE pulsar)
  install(TARGETS ${test_name} DESTINATION ${dir})
  pulsar_runtest(${test_name}_CPP $<TARGET_FILE:${test_name}>) 
endfunction()

#Macro for defining both a Python and C++ test
function(pulsar_test dir test_name)
   pulsar_py_test(${dir} ${test_name})
   pulsar_cxx_test(${dir} ${test_name})
endfunction()

########################
# Core pulsar module tests
########################
foreach(dir methods system_fragmenters)
  add_subdirectory(${dir})
  install(FILES ${CMAKE_BINARY_DIR}/${dir}/CTestTestfile.cmake DESTINATION ${dir})
endforeach()

# This file was created from the various add_test commands
install(FILES ${CMAKE_BINARY_DIR}/CTestTestfile.cmake DESTINATION .)
