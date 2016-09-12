#Each method should append their source files to this variable
set(BPMETHODS_SRC creator.cpp export.cpp)

#Each method should append their include files to this variable
set(BPMETHODS_INC "")

add_subdirectory(CompositeMethods)
add_subdirectory(MBE)
add_subdirectory(MethodHelpers)
add_subdirectory(SCF)


add_library(Methods MODULE ${BPMETHODS_SRC})
target_compile_options(Methods PRIVATE ${PULSAR_CXX_STRICT_FLAGS})
target_link_libraries(Methods pulsar Common)
target_include_directories(Methods SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIR})
set_target_properties(Methods PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/Common")


install(TARGETS Methods
        LIBRARY DESTINATION Methods
        RUNTIME DESTINATION Methods)

install(FILES "__init__.py" DESTINATION Methods)
install(FILES "modinfo.py"  DESTINATION Methods)
install(FILES "creator.py"  DESTINATION Methods)
