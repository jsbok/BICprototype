cmake -DCMAKE_INSTALL_PREFIX=./install -S . -B build
cd build && make -j8 && make install
