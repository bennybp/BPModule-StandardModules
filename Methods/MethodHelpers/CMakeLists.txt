set(MBE_FILES 
        ${CMAKE_CURRENT_SOURCE_DIR}/MethodHelpers.cpp
 )

set(BPMETHODS_SRC
       ${BPMETHODS_SRC}
       ${MBE_FILES}
    PARENT_SCOPE
)
