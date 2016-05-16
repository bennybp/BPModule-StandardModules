set(ONEELECTRONINTEGRALS_FILES creator.cpp
                          Overlap.cpp
                          KineticEnergy.cpp
                          OneElectronPotential.cpp
                          OneElectronPotential_LUT.cpp
                          CoreBuild.cpp
   )

set (ONEELECTRONINTEGRALS_PYFILES
    )


add_library(OneElectronIntegrals MODULE ${ONEELECTRONINTEGRALS_FILES})

target_link_libraries(OneElectronIntegrals pulsar Integral_Common)
target_compile_options(OneElectronIntegrals PRIVATE ${PULSAR_CXX_STRICT_FLAGS})

install(TARGETS OneElectronIntegrals LIBRARY DESTINATION OneElectronIntegrals
                                     RUNTIME DESTINATION OneElectronIntegrals)


install(FILES "__init__.py" DESTINATION OneElectronIntegrals)
install(FILES "modinfo.py"  DESTINATION OneElectronIntegrals)
install(FILES ${ONEELECTRONINTEGRALS_PYFILES} DESTINATION OneElectronIntegrals)

set_target_properties(OneElectronIntegrals PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/StandardModules)
