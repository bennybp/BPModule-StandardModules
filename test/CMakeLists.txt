cmake_minimum_required(VERSION 3.0)
#Should find the Pulsar in the staging area
#This sets PULSAR_RUNTEST variable
find_package(pulsar REQUIRED)
include(pulsar/CTestMacros)
########################
# Core pulsar module tests
########################
foreach(dir methods system_fragmenters)
  add_subdirectory(${dir})
  install(FILES ${CMAKE_BINARY_DIR}/${dir}/CTestTestfile.cmake DESTINATION ${dir})
endforeach()

add_subdirectory(testmodules)

# This file was created from the various add_test commands
install(FILES ${CMAKE_BINARY_DIR}/CTestTestfile.cmake DESTINATION .)
