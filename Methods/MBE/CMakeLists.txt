set(MBE_FILES creator.cpp
              MBE.cpp
   )

add_library(MBE MODULE ${MBE_FILES})

target_link_libraries(MBE bpmodule)

#set_target_properties(MBE PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)


install(TARGETS MBE
        LIBRARY DESTINATION MBE
        RUNTIME DESTINATION MBE)

install(FILES "__init__.py" DESTINATION MBE)
install(FILES "modinfo.py"  DESTINATION MBE)
