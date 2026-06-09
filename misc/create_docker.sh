#docker build -t ghcr.io/casadi/ci-swig:latest -f Dockerfile.swig .

# wasm-js build toolchain (flang->wasm + emsdk); used by binaries.yml `wasm` job
docker build -t ghcr.io/casadi/build-wasm:latest -f Dockerfile.wasm .
docker push ghcr.io/casadi/build-wasm:latest

#docker build -t ghcr.io/casadi/coinbuild:latest -f Dockerfile.coinbuild .
#cp ../docs/users_guide/requirements.txt doc_requirements.txt
#docker build -t ghcr.io/casadi/ci-doc:latest -f Dockerfile.doc .

docker build -t ghcr.io/casadi/web:latest -f Dockerfile.web .
# docker push

docker build -t ghcr.io/jgillis/windows-shared-x64-posix:production -f Dockerfile.windows .
docker push ghcr.io/jgillis/windows-shared-x64-posix:production

docker tag dockcross/manylinux2014-aarch64:latest ghcr.io/jgillis/manylinux2014-aarch64:production
