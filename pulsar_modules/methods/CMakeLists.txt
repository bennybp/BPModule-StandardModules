#Each method should append their source files to this variable
set(PULSAR_METHODS_SRC export.cpp)

#Each method should append their include files to this variable
set(PULSAR_METHODS_INC "")

#At the moment there is nothing to build in composite_methods
foreach(sub_dir mbe method_helpers optimizer scf)
    add_subdirectory(${sub_dir})
endforeach()

set(PULSAR_METHODS_SRC ${PULSAR_METHODS_SRC} PARENT_SCOPE)
