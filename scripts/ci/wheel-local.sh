#!/usr/bin/env bash

cd "$(dirname "${BASH_SOURCE[0]}")"
set -ex

image=alpaqa-wheel-local-img6
container=alpaqa-wheel-local-cnt6

if [ -z $(docker image ls -q $image) ]; then
    docker create --interactive \
        -v"$PWD/../../:/mnt" \
        --name $container \
        tttapa/alpaqa-build-python-gcc:3.9-11 bash
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
python -m pip install -U build
rm -rf /tmp/dist
mkdir /tmp/dist
FC=gfortran \
CXXFLAGS="-march=skylake -static-libstdc++ -static-libgcc" \
CFLAGS="-march=skylake -static-libgcc" \
LDFLAGS="-static-libstdc++ -static-libgcc" \
    python -m build --outdir /tmp/dist
cpv=$(echo $PYTHON_VERSION | awk -F. '{print $1 $2}')
LD_LIBRARY_PATH=$VIRTUAL_ENV/lib \
    auditwheel repair --plat manylinux_2_27_x86_64 \
        /tmp/dist/alpaqa-*.whl
cp /tmp/dist/alpaqa-*.tar.gz wheelhouse
EOF
docker stop $container