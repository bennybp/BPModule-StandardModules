# Create test libraries
add_library(dynlib SHARED dynlib.cpp)
add_library(staticlib STATIC staticlib.cpp)
set_target_properties(staticlib PROPERTIES COMPILE_FLAGS "-fPIC")


set(TESTMODULES_FILES creator.cpp
                      TestModule1.cpp
                      TestExtLib.cpp
   )

set (TESTMODULES_PYFILES creator.py
                         TestOptions.py
                         TestPyModule1.py
    )



add_library(TestModules MODULE ${TESTMODULES_FILES})

target_link_libraries(TestModules dynlib staticlib pulsar)

install(TARGETS TestModules LIBRARY DESTINATION TestModules
                            RUNTIME DESTINATION TestModules)

install(FILES "__init__.py" DESTINATION TestModules)
install(FILES "modinfo.py"  DESTINATION TestModules)
install(FILES ${TESTMODULES_PYFILES} DESTINATION TestModules)
                            
set_target_properties(TestModules PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/TestModules)

# install the dynamic library
install(TARGETS dynlib LIBRARY DESTINATION TestModules
                       RUNTIME DESTINATION TestModules)
