set(BPTEST_FILES 
        ${CMAKE_CURRENT_SOURCE_DIR}/Damping.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/DIIS.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SCFCommon.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/CoreGuess.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/HFIterate.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/BasicFockBuild.cpp
 )

set(BPMETHODS_SRC
       ${BPMETHODS_SRC}
       ${BPTEST_FILES}
    PARENT_SCOPE
)

