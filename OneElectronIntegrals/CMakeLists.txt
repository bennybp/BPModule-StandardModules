set(ONEELECTRONINTEGRALS_FILES creator.cpp
                          Overlap.cpp
   )

set (ONEELECTRONINTEGRALS_PYFILES
    )


add_library(OneElectronIntegrals MODULE ${ONEELECTRONINTEGRALS_FILES})

target_link_libraries(OneElectronIntegrals pulsar)

install(TARGETS OneElectronIntegrals LIBRARY DESTINATION OneElectronIntegrals
                                RUNTIME DESTINATION OneElectronIntegrals)


install(FILES "__init__.py" DESTINATION OneElectronIntegrals)
install(FILES "modinfo.py"  DESTINATION OneElectronIntegrals)
install(FILES ${ONEELECTRONINTEGRALS_PYFILES} DESTINATION OneElectronIntegrals)
