install(FILES "CP.py" DESTINATION Methods/CompositeMethods)
install(FILES "FPA.py"  DESTINATION Methods/CompositeMethods)
install(FILES "FellerCBS.py"  DESTINATION Methods/CompositeMethods)
install(FILES "HelgakerCBS.py"  DESTINATION Methods/CompositeMethods)
set(MIM_FILES 
       ${CMAKE_CURRENT_SOURCE_DIR}/MIM.cpp
 )

set(BPMETHODS_SRC
       ${BPMETHODS_SRC}
       ${MIM_FILES}
    PARENT_SCOPE
)
