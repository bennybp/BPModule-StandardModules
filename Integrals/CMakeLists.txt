set(INTEGRAL_COMMON_FILES
                 Common.cpp
                 Boys.cpp
   )


#target_link_libraries(Integral_Common pulsar)
#add_library(Integral_Common OBJECT ${INTEGRAL_COMMON_FILES})


# See CMAKE bug https://cmake.org/Bug/view.php?id=14778
add_library(Integral_Common SHARED ${INTEGRAL_COMMON_FILES})

target_link_libraries(Integral_Common pulsar)



add_subdirectory(OneElectronIntegrals)
add_subdirectory(TwoElectronIntegrals)


install(TARGETS Integral_Common LIBRARY DESTINATION StandardModules
                                RUNTIME DESTINATION StandardModules)
