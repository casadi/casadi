#export TARGET=manylinux2014-x64
export TARGET=manylinux1-x64
#export TARGET=windows-shared-x64-posix
#docker run --rm dockcross/$TARGET:latest > dockcross.$TARGET
#chmod +x dockcross.$TARGET

docker run --rm ghcr.io/jgillis/$TARGET:production > dockcross.$TARGET
chmod +x dockcross.$TARGET

export FLAGS="-DWITH_COMMON=OFF -DWITH_BUILD_REQUIRED=ON -DWITH_BUILD_BONMIN=OFF -DWITH_BONMIN=OFF -DWITH_IPOPT=ON -DWITH_BUILD_LAPACK=ON -DWITH_LAPACK=ON -DWITH_MUMPS=ON -DWITH_CLP=OFF -DWITH_CBC=OFF -DWITH_THREAD=ON -DWITH_QPOASES=ON -DWITH_HPIPM=ON -DWITH_BLASFEO=ON -DWITH_BUILD_HPIPM=ON -DWITH_BUILD_BLASFEO=ON -DWITH_HIGHS=OFF -DWITH_BUILD_HIGHS=OFF -DWITH_BUILD_SPRAL=ON -DWITH_SPRAL=ON"
export FLAGS=""

docker run --rm -v`pwd`:/local ghcr.io/casadi/ci-swig:latest /bin/bash -c "mkdir build-temp && cd build-temp && cmake -DWITH_SELFCONTAINED=ON -DWITH_PYTHON=ON -DSWIG_EXPORT=ON -DWITH_COMMON=OFF .. && make python_source && cd .. && rm -rf build-temp"

#docker run --rm -v`pwd`:/local ghcr.io/casadi/ci-swig:latest /bin/bash -c "mkdir build-temp && cd build-temp && cmake -DWITH_SELFCONTAINED=ON -DWITH_PYTHON3=ON -DWITH_PYTHON=ON -DSWIG_EXPORT=ON -DWITH_COMMON=OFF .. && make python_source && cd .. && rm -rf build-temp"

#docker run --rm -v`pwd`:/local ghcr.io/casadi/ci-swig:latest /bin/bash -c "mkdir build-temp && cd build-temp && cmake -DWITH_SELFCONTAINED=ON -DWITH_MATLAB=ON -DSWIG_EXPORT=ON -DWITH_COMMON=OFF .. && make matlab_source && cd .. && rm -rf build-temp"

export PYTHON_INCLUDE_DIR=/opt/python/cp27-cp27m/include/python2.7/
mkdir -p build-local-$TARGET
./dockcross.$TARGET cmake -Bbuild-local-$TARGET -DWITH_SELFCONTAINED=ON $FLAGS -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=/work/install-$TARGET -DSWIG_IMPORT=ON -DWITH_PYTHON=ON -DPYTHON_LIBRARIES=$PYTHON_LIBRARIES -DPYTHON_INCLUDE_DIR=$PYTHON_INCLUDE_DIR -H.
./dockcross.$TARGET --args "--env PKG_CONFIG_PATH_x86_64_w64_mingw32_shared_posix=/work/build/external_projects/lib64/pkgconfig:/work/build/external_projects/lib/pkgconfig:/work/build/external_projects/share/pkgconfig" -- cmake --build build-local-$TARGET --parallel 12 -v
./dockcross.$TARGET cmake --build build-local-$TARGET --target install -v
