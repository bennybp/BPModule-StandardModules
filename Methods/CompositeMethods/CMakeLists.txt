foreach(da_file "CP.py" "FPA.py" "FellerCBS.py" "HelgakerCBS.py" 
                "MyCrazyCompositeMethod.py")
    install(FILES ${da_file} DESTINATION Methods/CompositeMethods)
endforeach()