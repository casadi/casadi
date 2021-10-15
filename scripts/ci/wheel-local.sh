#!/usr/bin/env bash

cd "$(dirname "${BASH_SOURCE[0]}")"
set -ex

image=panocpy-wheel-local-img6
container=panocpy-wheel-local-cnt6

if [ -z $(docker image ls -q $image) ]; then
    docker create --interactive \
        -v"$PWD/../../:/mnt" \
        --name $container \
        tttapa/panoc-alm-build-python-gcc:3.9-11 bash
    docker start $container
    docker exec -i $container bash << 'EOF'
set -e
apt-get update -y
apt-get install -y autotools-dev automake libpcre3-dev libopenblas-dev
apt-get clean autoclean
apt-get autoremove -y
rm -rf /var/lib/apt/lists/*
export CXXFLAGS="-march=skylake -static-libstdc++ -static-libgcc"
export LDFLAGS="-static-libstdc++ -static-libgcc"
export CFLAGS="-march=skylake -static-libgcc"
cd /mnt
python3 -m venv /tmp/py-venv
. /tmp/py-venv/bin/activate
python -m pip install -r scripts/requirements.txt -r scripts/requirements-wheel.txt
FC=gfortran \
    bash scripts/ci/install-dependencies-wheel.sh
EOF
    docker commit $container $image
else
    if [ -z $(docker ps -aqf name="$container") ]; then 
        docker create --interactive \
            -v"$PWD/../../:/mnt" \
            --name $container \
            $image bash
    fi
    docker start $container
fi

docker exec -i $container bash << 'EOF'
set -e
if [ ! -e /usr/local/bin/patchelf ]; then
    pushd /tmp
    git clone https://github.com/NixOS/patchelf --branch 0.13
    pushd patchelf
    ./bootstrap.sh
    ./configure
    make
    make check
    make install
    popd
    popd
fi
cd /mnt
. /tmp/py-venv/bin/activate
rm -rf _skbuild/
FC=gfortran \
CXXFLAGS="-march=skylake -static-libstdc++ -static-libgcc" \
CFLAGS="-march=skylake -static-libgcc" \
LDFLAGS="-static-libstdc++ -static-libgcc" \
    python setup.py bdist_wheel --build-type RelWithDebInfo -j$(nproc) --generator Ninja --skip-generator-test
cpv=$(echo $PYTHON_VERSION | awk -F. '{print $1 $2}')
LD_LIBRARY_PATH=$VIRTUAL_ENV/lib \
    auditwheel repair --plat manylinux_2_27_x86_64 \
        dist/panocpy-0.0.2-cp$cpv-cp$cpv-linux_x86_64.whl
EOF
docker stop $container