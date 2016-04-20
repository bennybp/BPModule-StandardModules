# Pulsar-StandardModules

Quick and dirty build:

    mkdir build
    cd build
    CC=compiler CXX=compiler cmake ../ \
                -DCMAKE_INSTALL_PREFIX=/pulsar/install/path/modules/ \
                -DPULSAR_PATH=/pulsar/install/path/

    make install
