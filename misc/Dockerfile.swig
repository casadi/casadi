FROM ubuntu:24.04

# Circumvent bug in SWIG
RUN apt-get update && DEBIAN_FRONTEND=noninteractive  apt install -y g++ gcc libpcre3-dev automake yodl byacc autoconf git cmake python3-dev libpcre2-dev bison

RUN git clone https://github.com/jaeandersson/swig.git && \
    cd swig && git checkout 38d80f608d5cb74138e445b1507a5bf40293f628 && ./autogen.sh && ./configure && make && make install && ./configure --prefix=/opt/swig && make && make install

WORKDIR /local
