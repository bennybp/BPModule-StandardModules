set(SYSTEMINTEGRALS_FILES creator.cpp
                          NuclearRepulsion.cpp
                          NuclearDipole.cpp
   )

set (SYSTEMINTEGRALS_PYFILES
    )


add_library(SystemIntegrals MODULE ${SYSTEMINTEGRALS_FILES})

target_link_libraries(SystemIntegrals pulsar)
target_compile_options(SystemIntegrals PRIVATE ${PULSAR_CXX_STRICT_FLAGS})

install(TARGETS SystemIntegrals LIBRARY DESTINATION SystemIntegrals
                                RUNTIME DESTINATION SystemIntegrals)


install(FILES "__init__.py" DESTINATION SystemIntegrals)
install(FILES "modinfo.py"  DESTINATION SystemIntegrals)
install(FILES ${SYSTEMINTEGRALS_PYFILES} DESTINATION SystemIntegrals)
