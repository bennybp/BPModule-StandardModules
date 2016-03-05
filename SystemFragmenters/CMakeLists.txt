set(SYSTEMFRAG_INSTALL_PATH "modules/SystemFragmenters")

set(SYSTEMFRAG_FILES creator.cpp
                     Atomizer.cpp
   )

add_library(SystemFragmenters MODULE ${SYSTEMFRAG_FILES})

target_link_libraries(SystemFragmenters bpmodule)


install(TARGETS SystemFragmenters
        LIBRARY DESTINATION ${SYSTEMFRAG_INSTALL_PATH}
        RUNTIME DESTINATION ${SYSTEMFRAG_INSTALL_PATH})

install(FILES "__init__.py" DESTINATION ${SYSTEMFRAG_INSTALL_PATH})
install(FILES "modinfo.py"  DESTINATION ${SYSTEMFRAG_INSTALL_PATH})
