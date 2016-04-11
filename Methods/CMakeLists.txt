#Each method should append their source files to this variable
set(BPMETHODS_SRC creator.cpp)

#Each method should append their include files to this variable
set(BPMETHODS_INC "")

add_subdirectory(BSSE)
add_subdirectory(CBS)
add_subdirectory(MBE)
add_subdirectory(MIM)
add_subdirectory(MP2)
add_subdirectory(SCF)


add_library(Methods MODULE ${BPMETHODS_SRC})

target_link_libraries(Methods bpmodule)


install(TARGETS Methods
        LIBRARY DESTINATION Methods
        RUNTIME DESTINATION Methods)

install(FILES "__init__.py" DESTINATION Methods)
install(FILES "modinfo.py"  DESTINATION Methods)
install(FILES "creator.py"  DESTINATION Methods)