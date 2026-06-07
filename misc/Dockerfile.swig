FROM ubuntu:24.04

# Circumvent bug in SWIG
RUN apt-get update && DEBIAN_FRONTEND=noninteractive  apt install -y g++ gcc libpcre3-dev automake yodl byacc autoconf git cmake python3-dev libpcre2-dev bison

# matlab-customdoc2 carries the python + matlab + wasm-js SWIG backends.
# The `generate wasm-js` step in binaries.yml needs -wasm-js in this swig.
# (Pin to a commit once the wasm-js fixes are pushed to this branch.)
RUN git clone https://github.com/jaeandersson/swig.git && \
    cd swig && git checkout matlab-customdoc2 && ./autogen.sh && ./configure && make && make install && ./configure --prefix=/opt/swig && make && make install

WORKDIR /local
