set(INTEGRALS_FILES creator.cpp
                    OSOverlapTerms.cpp
                    OSOverlap.cpp
                    OSKineticEnergy.cpp
                    OSDipole.cpp
                    OSOneElectronPotential.cpp
                    OSOneElectronPotential_LUT.cpp

                    OneElectronProperty.cpp
                    OneElectron_Eigen.cpp
                    OneElectronIntegralSum.cpp
                    ReferenceERI.cpp
                    ValeevRef.cpp
                    NuclearRepulsion.cpp
                    NuclearDipole.cpp

                    boys/Boys_shortgrid.cpp
                    boys/Boys_longfac.cpp
   )

add_library(Integrals MODULE ${INTEGRALS_FILES})

target_link_libraries(Integrals pulsar Common)
target_compile_options(Integrals PRIVATE ${PULSAR_CXX_STRICT_FLAGS})
target_include_directories(Integrals SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIR})
set_target_properties(Integrals PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/Common")

install(TARGETS Integrals LIBRARY DESTINATION Integrals
                          RUNTIME DESTINATION Integrals)

install(FILES "__init__.py" DESTINATION Integrals)
install(FILES "modinfo.py"  DESTINATION Integrals)
