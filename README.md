# Pulsar-StandardModules
This is the repo for the reference implementation of Pulsar's modules.  For the time being there aren't any other options aside from this reference implementation, so this is a requirement for using Pulsar.  

## Installation
The repo has no dependencies aside from those of Pulsar-Core itself.  This means that aside from the standard CMake Variables, the only thing you should have to specify is where you installed the core and where you want to install these modules.  By convention, all Pulsar modules are being installed into a sub-directory of the main core instillation called "modules" and you should set `CMAKE_INSTALL_PREFIX` accordingly.

TODO: If the last statement is our official policy do we want to make that the default in the scripts and then print a warning if the user changes it?

    cmake -Bbuild -H. \
          -DCMAKE_CXX_COMPILER=/Your/CXX_Compiler \
          -DCMAKE_C_COMPILER=/Your/C_Compiler \
          -DMPI_C_COMPILER=/Your/MPI_C_Compiler \
          -DMPI_CXX_COMPILER=/Your/MPI_CXX_COMPILER \
          -DCMAKE_INSTALL_PREFIX=/pulsar/install/path/modules/ \
          -DPULSAR_PATH=/where/you/installed/pulsar/core/

    make
    make install
    
Note: As part of finding Pulsar-Core CMake will have to find some system libraries like Python and Eigen.  If you specified these manually while installing the core you will need to specify them again when installing the Standard Modules to ensure you use the same versions.
