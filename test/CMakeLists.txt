cmake_minimum_required(VERSION 3.0)
find_package(pulsar REQUIRED)
include(CTestMacros)
########################
# Core pulsar module tests
########################

#Macro for defining a python test
function(pulsar_sm_py_test dir test_name)
  pulsar_py_test(${dir} ${test_name} ${STAGE_DIR}${CMAKE_INSTALL_PREFIX})
endfunction()

foreach(dir analysis methods system_fragmenters)
  add_subdirectory(${dir})
  install(FILES ${CMAKE_BINARY_DIR}/${dir}/CTestTestfile.cmake
          DESTINATION ${dir})
endforeach()

add_subdirectory(testmodules)

# This file was created from the various add_test commands
install(FILES ${CMAKE_BINARY_DIR}/CTestTestfile.cmake DESTINATION .)
