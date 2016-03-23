# BPModule-StandardModules

Quick and dirty build:

    mkdir build
    cd build
    CC=compiler CXX=compiler cmake ../ \
                -DCMAKE_INSTALL_PREFIX=/bpmodule/install/path/modules/ \
                -DBPMODULE_PATH=/bpmodule/install/path/

    make install
