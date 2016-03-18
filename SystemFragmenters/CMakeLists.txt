set(SYSTEMFRAG_FILES creator.cpp
                     Atomizer.cpp
                     Bondizer.cpp
   )

add_library(SystemFragmenters MODULE ${SYSTEMFRAG_FILES})

target_link_libraries(SystemFragmenters bpmodule)


install(TARGETS SystemFragmenters
        LIBRARY DESTINATION SystemFragmenters
        RUNTIME DESTINATION SystemFragmenters)

install(FILES "__init__.py" DESTINATION SystemFragmenters)
install(FILES "modinfo.py"  DESTINATION SystemFragmenters)
