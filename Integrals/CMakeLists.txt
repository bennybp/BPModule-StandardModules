set(INTEGRAL_COMMON_FILES
                 Common.cpp
   )


#target_link_libraries(Integral_Common pulsar)
#add_library(Integral_Common OBJECT ${INTEGRAL_COMMON_FILES})


# See CMAKE bug https://cmake.org/Bug/view.php?id=14778
add_library(Integral_Common STATIC ${INTEGRAL_COMMON_FILES})

#TODO - a better way to do this? shared library?
target_compile_options(Integral_Common PRIVATE "-fPIC")
target_link_libraries(Integral_Common pulsar)



add_subdirectory(OneElectronIntegrals)
add_subdirectory(TwoElectronIntegrals)


# Nothing from the common library should be installed
