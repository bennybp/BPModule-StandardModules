set(TWOELECTRONINTEGRALS_FILES creator.cpp
                          ReferenceERI.cpp
                          ValeevRef.cpp
   )

set (TWOELECTRONINTEGRALS_PYFILES
    )


add_library(TwoElectronIntegrals MODULE ${TWOELECTRONINTEGRALS_FILES})

target_link_libraries(TwoElectronIntegrals pulsar Integral_Common)
target_compile_options(TwoElectronIntegrals PRIVATE ${PULSAR_CXX_STRICT_FLAGS})

install(TARGETS TwoElectronIntegrals LIBRARY DESTINATION TwoElectronIntegrals
                                     RUNTIME DESTINATION TwoElectronIntegrals)


install(FILES "__init__.py" DESTINATION TwoElectronIntegrals)
install(FILES "modinfo.py"  DESTINATION TwoElectronIntegrals)
install(FILES ${TWOELECTRONINTEGRALS_PYFILES} DESTINATION TwoElectronIntegrals)

set_target_properties(TwoElectronIntegrals PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/StandardModules)
