set(MBE_FILES 
        ${CMAKE_CURRENT_SOURCE_DIR}/MBE.cpp
 )

set(BPMETHODS_SRC
       ${BPMETHODS_SRC}
       ${MBE_FILES}
    PARENT_SCOPE
)

install(FILES "MBEHelper.py"  DESTINATION Methods/MBE)