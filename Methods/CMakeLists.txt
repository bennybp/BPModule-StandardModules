#Each method should append their source files to this variable
set(BPMETHODS_SRC creator.cpp)

#Each method should append their include files to this variable
set(BPMETHODS_INC "")

add_subdirectory(CC)
add_subdirectory(CompositeMethods)
add_subdirectory(MBE)
add_subdirectory(MethodHelpers)
add_subdirectory(MP2)
add_subdirectory(SCF)


add_library(Methods MODULE ${BPMETHODS_SRC})
target_compile_options(Methods PRIVATE ${PULSAR_CXX_STRICT_FLAGS})
target_link_libraries(Methods pulsar)
target_include_directories(Methods SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIR})


install(TARGETS Methods
        LIBRARY DESTINATION Methods
        RUNTIME DESTINATION Methods)

install(FILES "__init__.py" DESTINATION Methods)
install(FILES "modinfo.py"  DESTINATION Methods)
install(FILES "creator.py"  DESTINATION Methods)
