#!/usr/bin/env bash

set -ex

versions=( $1 )

update-alternatives --remove "gcc"        "/usr/bin/gcc"
update-alternatives --remove "g++"        "/usr/bin/g++"
update-alternatives --remove "gcov"       "/usr/bin/gcov"
update-alternatives --remove "gcc-ar"     "/usr/bin/gcc-ar"
update-alternatives --remove "gcc-ranlib" "/usr/bin/gcc-ranlib"
update-alternatives --remove "cpp"        "/usr/bin/cpp"
update-alternatives --remove "gfortran"   "/usr/bin/gfortran"

for v in "${versions[@]}"; do
    update-alternatives \
        --install "/usr/bin/gcc"      "gcc" "$(which gcc-$v)" $(($v * 10)) \
        --slave "/usr/bin/g++"        "g++" "$(which g++-$v)" \
        --slave "/usr/bin/gcov"       "gcov" "$(which gcov-$v)" \
        --slave "/usr/bin/gcc-ar"     "gcc-ar" "$(which gcc-ar-$v)" \
        --slave "/usr/bin/gcc-ranlib" "gcc-ranlib" "$(which gcc-ranlib-$v)" \
        --slave "/usr/bin/cpp"        "cpp" "$(which cpp-$v)"
    
    if which gfortran-$v; then 
        update-alternatives \
            --install "/usr/bin/gfortran" \
                "gfortran" "$(which gfortran-$v)" $(($v * 10))
    fi
done
