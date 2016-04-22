set(SYSTEMINTEGRALS_FILES creator.cpp
                          NuclearRepulsion.cpp
   )

set (SYSTEMINTEGRALS_PYFILES
    )


add_library(SystemIntegrals MODULE ${SYSTEMINTEGRALS_FILES})

target_link_libraries(SystemIntegrals pulsar)

install(TARGETS SystemIntegrals LIBRARY DESTINATION SystemIntegrals
                                RUNTIME DESTINATION SystemIntegrals)


install(FILES "__init__.py" DESTINATION SystemIntegrals)
install(FILES "modinfo.py"  DESTINATION SystemIntegrals)
install(FILES ${SYSTEMINTEGRALS_PYFILES} DESTINATION SystemIntegrals)
