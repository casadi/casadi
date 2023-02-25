#!/bin/bash

set -o pipefail
set -e

export TARGET=manylinux2014-x64
#export TARGET=manylinux1-x64
#export TARGET=windows-shared-x64-posix
#docker run --rm dockcross/$TARGET:latest > dockcross.$TARGET
#chmod +x dockcross.$TARGET

echo "TARGET: $TARGET"

#docker run --rm ghcr.io/jgillis/$TARGET:production > dockcross.$TARGET
#chmod +x dockcross.$TARGET

INT_TARGET=python
#INT_TARGET=matlab
#INT_TARGET=octave
#INT_TARGET=cpp

#INT_TARGET=octave

BUILD_TYPE=Debug

export FLAGS="-DWITH_COMMON=OFF -DWITH_BUILD_REQUIRED=ON -DWITH_BUILD_BONMIN=OFF -DWITH_BONMIN=OFF -DWITH_IPOPT=ON -DWITH_BUILD_LAPACK=ON -DWITH_LAPACK=ON -DWITH_MUMPS=ON -DWITH_CLP=OFF -DWITH_CBC=OFF -DWITH_THREAD=ON -DWITH_QPOASES=ON -DWITH_HPIPM=ON -DWITH_BLASFEO=ON -DWITH_BUILD_HPIPM=ON -DWITH_BUILD_BLASFEO=ON -DWITH_HIGHS=OFF -DWITH_BUILD_HIGHS=OFF -DWITH_BUILD_SPRAL=OFF -DWITH_SPRAL=OFF  -DWITH_PROXQP=OFF -DWITH_BUILD_PROXQP=OFF -DWITH_BUILD_EIGEN3=ON -DWITH_BUILD_SIMDE=ON -DWITH_BUILD_SUNDIALS=OFF -DWITH_SUNDIALS=OFF -DWITH_OSQP=ON -DWITH_BUILD_OSQP=ON -DWITH_SUPERSCS=ON -DWITH_BUILD_SUPERSCS=ON -DWITH_KNITRO=ON -DWITH_CPLEX=ON"


BUILD_TYPE_MOCKUPS=Release

curl -OL https://github.com/casadi/mockups/releases/download/v45/mockups_${TARGET}_${BUILD_TYPE_MOCKUPS}.zip

rm -rf mockups-${TARGET}-${BUILD_TYPE_MOCKUPS}
unzip mockups_${TARGET}_${BUILD_TYPE_MOCKUPS}.zip -d mockups-${TARGET}-${BUILD_TYPE_MOCKUPS}


rm -f  build-local-$TARGET/CMakeCache.txt
mkdir -p build-local-$TARGET

echo "INT_TARGET: $INT_TARGET"

case $INT_TARGET in
  python)
    export PYTHON_INCLUDE_DIR=/opt/python/cp38-cp38/include/python3.8/
    if [[ "$PYTHON_INCLUDE_DIR" == *py2* ]]; then
        docker run --rm -v`pwd`:/local ghcr.io/casadi/ci-swig:latest /bin/bash -c "mkdir build-temp && cd build-temp && cmake -DWITH_SELFCONTAINED=ON -DWITH_PYTHON=ON -DSWIG_EXPORT=ON -DWITH_COMMON=OFF .. && make python_source && cd .. && rm -rf build-temp"
    else
        docker run --rm -v`pwd`:/local ghcr.io/casadi/ci-swig:latest /bin/bash -c "mkdir build-temp && cd build-temp && cmake -DWITH_SELFCONTAINED=ON -DWITH_PYTHON3=ON -DWITH_PYTHON=ON -DSWIG_EXPORT=ON -DWITH_COMMON=OFF .. && make python_source && cd .. && rm -rf build-temp"
    fi
    ./dockcross.$TARGET cmake -Bbuild-local-$TARGET -DSKIP_CONFIG_H_GENERATION=ON -DWITH_SELFCONTAINED=ON $FLAGS -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_INSTALL_PREFIX=/work/install-$INT_TARGET-$TARGET -DSWIG_IMPORT=ON -DWITH_PYTHON=ON -DCMAKE_PREFIX_PATH=/work/mockups-${TARGET}-${BUILD_TYPE_MOCKUPS}/cmake -DPYTHON_LIBRARIES=$PYTHON_LIBRARIES -DPYTHON_INCLUDE_DIR=$PYTHON_INCLUDE_DIR -H.
    ;;

  matlab)
    docker run --rm -v`pwd`:/local ghcr.io/casadi/ci-swig:latest /bin/bash -c "mkdir build-temp && cd build-temp && cmake -DWITH_SELFCONTAINED=ON -DWITH_MATLAB=ON -DSWIG_EXPORT=ON -DWITH_COMMON=OFF .. && make matlab_source && cd .. && rm -rf build-temp"
    ./dockcross.$TARGET cmake -Bbuild-local-$TARGET -DSKIP_CONFIG_H_GENERATION=ON -DWITH_SELFCONTAINED=ON $FLAGS -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_INSTALL_PREFIX=/work/install-$INT_TARGET-$TARGET -DSWIG_IMPORT=ON -DWITH_MATLAB=ON -DCMAKE_PREFIX_PATH=/work/mockups-${TARGET}-${BUILD_TYPE_MOCKUPS}/cmake -H.
    ;;

  octave)
    docker run --rm -v`pwd`:/local ghcr.io/casadi/ci-swig:latest /bin/bash -c "mkdir build-temp && cd build-temp && cmake -DWITH_SELFCONTAINED=ON -DWITH_MATLAB=ON -DSWIG_EXPORT=ON -DWITH_COMMON=OFF .. && make matlab_source && cd .. && rm -rf build-temp"
    ./dockcross.$TARGET cmake -Bbuild-local-$TARGET -DSKIP_CONFIG_H_GENERATION=ON -DWITH_SELFCONTAINED=ON $FLAGS -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_INSTALL_PREFIX=/work/install-$INT_TARGET-$TARGET -DSWIG_IMPORT=ON -DWITH_OCTAVE=ON -DWITH_OCTAVE_IMPORT=ON -DCMAKE_PREFIX_PATH=/work/mockups-$TARGET-${BUILD_TYPE_MOCKUPS}/cmake -H.
    ;;

  cpp)
    ./dockcross.$TARGET cmake -Bbuild-local-$TARGET -DWITH_SELFCONTAINED=ON $FLAGS -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_PREFIX_PATH=/work/mockups-$TARGET/cmake -DCMAKE_INSTALL_PREFIX=/work/install-$INT_TARGET-$TARGET -H.
    ;;

  *)
    echo "Not supported" && exit 1
    ;;
esac

cp mockups-${TARGET}-${BUILD_TYPE_MOCKUPS}/cplex/lib/* install-$INT_TARGET-$TARGET/casadi

./dockcross.$TARGET --args "--env PKG_CONFIG_PATH_x86_64_w64_mingw32_shared_posix=/work/build/external_projects/lib64/pkgconfig:/work/build/external_projects/lib/pkgconfig:/work/build/external_projects/share/pkgconfig" -- .github/workflows/patch_toolchain cmake --build build-local-$TARGET --parallel 12 -v
./dockcross.$TARGET cmake --build build-local-$TARGET --target install -v
