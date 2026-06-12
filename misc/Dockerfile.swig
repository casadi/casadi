FROM ubuntu:24.04

# Circumvent bug in SWIG
RUN apt-get update && DEBIAN_FRONTEND=noninteractive  apt install -y g++ gcc libpcre3-dev automake yodl byacc autoconf git cmake python3-dev libpcre2-dev bison

RUN git clone https://github.com/jaeandersson/swig.git && \
    cd swig && git checkout 42149b079f17df42112e79279c0c00e530cce273 && ./autogen.sh && ./configure && make && make install && ./configure --prefix=/opt/swig && make && make install

WORKDIR /local
