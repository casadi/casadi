FROM ubuntu:20.04

# Circumvent bug in SWIG
RUN apt-get update && DEBIAN_FRONTEND=noninteractive  apt install -y g++-7 gcc-7 libpcre3-dev automake yodl byacc autoconf git cmake python-dev python3-dev
ENV CXX=g++-7
ENV CC=gcc-7

RUN git clone https://github.com/jaeandersson/swig.git && \
    cd swig && git checkout matlab-customdoc && ./autogen.sh && ./configure && make && make install

WORKDIR /local
