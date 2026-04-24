FROM ubuntu:24.04

# Circumvent bug in SWIG
RUN apt-get update && DEBIAN_FRONTEND=noninteractive  apt install -y g++ gcc libpcre3-dev automake yodl byacc autoconf git cmake python3-dev libpcre2-dev bison

RUN git clone https://github.com/jaeandersson/swig.git && \
    cd swig && git checkout 223525dd46d810aa5c4c830805f5b80df1ba70b3 && ./autogen.sh && ./configure && make && make install && ./configure --prefix=/opt/swig && make && make install

WORKDIR /local
