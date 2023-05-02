#docker build -t ghcr.io/casadi/ci-swig:latest -f Dockerfile.swig .
#docker build -t ghcr.io/casadi/coinbuild:latest -f Dockerfile.coinbuild .
#cp ../docs/users_guide/requirements.txt doc_requirements.txt
#docker build -t ghcr.io/casadi/ci-doc:latest -f Dockerfile.doc .

docker build -t ghcr.io/casadi/web:latest -f Dockerfile.web .
# docker push

docker build -t ghcr.io/jgillis/windows-shared-x64-posix:production -f Dockerfile.windows .
docker push ghcr.io/jgillis/windows-shared-x64-posix:production
