set(SYSTEMFRAG_FILES creator.cpp
                     Atomizer.cpp
                     Bondizer.cpp
                     CrystalFragger.cpp
                     Ghoster.cpp
                     UserDefined.cpp
                     Null.cpp
   )

add_library(SystemFragmenters MODULE ${SYSTEMFRAG_FILES})
target_compile_options(SystemFragmenters PRIVATE ${PULSAR_CXX_STRICT_FLAGS})

target_link_libraries(SystemFragmenters pulsar)


install(TARGETS SystemFragmenters
        LIBRARY DESTINATION SystemFragmenters
        RUNTIME DESTINATION SystemFragmenters)

install(FILES "__init__.py" DESTINATION SystemFragmenters)
install(FILES "modinfo.py"  DESTINATION SystemFragmenters)
