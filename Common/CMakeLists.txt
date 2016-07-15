set(STANDARDMODULES_COMMON_FILES 
        BasisSetCommon.cpp
 )


add_library(Common SHARED ${STANDARDMODULES_COMMON_FILES})
target_compile_options(Common PRIVATE ${PULSAR_CXX_STRICT_FLAGS})
target_link_libraries(Common pulsar)
target_include_directories(Common SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIR})

install(TARGETS Common
        LIBRARY DESTINATION Common
        RUNTIME DESTINATION Common)

