cmake_minimum_required(VERSION 3.0)

list(APPEND CMAKE_PREFIX_PATH "${PULSAR_PATH}")
find_package(pulsar REQUIRED)

include_directories(${PULSAR_CXX_INCLUDES})
include_directories(${CMAKE_CURRENT_LIST_DIR})
add_subdirectory(SystemFragmenters)
add_subdirectory(Methods)
add_subdirectory(SystemIntegrals)
add_subdirectory(TestModules)

install(FILES "StandardModules.py" DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES "CallPsi4.py" DESTINATION ${CMAKE_INSTALL_PREFIX})
