ARG GCC_VERSION
ARG PYTHON_VERSION 

FROM ubuntu:bionic as gcc-install

ARG GCC_VERSION

RUN export DEBIAN_FRONTEND=noninteractive && \
    apt-get update -y && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    apt-get update -y && \
    add-apt-repository ppa:git-core/ppa && \
    apt-get update -y && \
    apt-get install -y \
        gcc-${GCC_VERSION} g++-${GCC_VERSION} gfortran-${GCC_VERSION} \
        autotools-dev automake libpcre3-dev libopenblas-dev \
        git zip unzip wget && \
    apt-get clean autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

COPY gcc.sh .
RUN chmod +x gcc.sh && ./gcc.sh ${GCC_VERSION}
RUN gcc --version

RUN wget https://github.com/ninja-build/ninja/releases/download/v1.11.0/ninja-linux.zip && \
    unzip ninja-linux.zip && \
    rm -f ninja-linux.zip && \
    mv ninja /usr/local/bin
RUN ninja --version

RUN wget https://github.com/Kitware/CMake/releases/download/v3.23.1/cmake-3.23.1-linux-x86_64.sh && \
    bash cmake-3.23.1-linux-x86_64.sh --skip-license --prefix=/usr/local && \
    rm -f cmake-3.23.1-linux-x86_64.sh
RUN cmake --version

ENV GCC_VERSION=${GCC_VERSION}

FROM gcc-install as casadi-install

COPY install-casadi.sh .
RUN chmod +x install-casadi.sh && ./install-casadi.sh /tmp/casadi-install

FROM gcc-install as casadi-install-static

COPY install-casadi-static.sh .
RUN chmod +x install-casadi-static.sh && ./install-casadi-static.sh /tmp/casadi-install

FROM gcc-install as eigen-install

COPY install-eigen.sh .
RUN chmod +x install-eigen.sh && ./install-eigen.sh /tmp/eigen-install

FROM gcc-install as gtest-install

COPY install-gtest.sh .
RUN chmod +x install-gtest.sh && ./install-gtest.sh /tmp/gtest-install

FROM gcc-install as alpaqa-build

COPY --from=casadi-install /tmp/casadi-install/ /usr/local/
COPY --from=eigen-install /tmp/eigen-install/ /usr/local/
COPY --from=gtest-install /tmp/gtest-install/ /usr/local/
RUN ldconfig

FROM ubuntu:impish as gcc-install-test
RUN export DEBIAN_FRONTEND=noninteractive && \
    apt-get update -y && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    apt-get update -y && \
    add-apt-repository ppa:git-core/ppa && \
    apt-get update -y && \
    apt-get install -y \
        gcc g++ gfortran ninja-build \
        python3-dev python3-pip python3-venv \
        autotools-dev automake libpcre3-dev libopenblas-dev \
        git zip unzip wget && \
    apt-get clean autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/Kitware/CMake/releases/download/v3.23.1/cmake-3.23.1-linux-x86_64.sh && \
    bash cmake-3.23.1-linux-x86_64.sh --skip-license --prefix=/usr/local && \
    rm -f cmake-3.23.1-linux-x86_64.sh
RUN cmake --version

FROM gcc-install-test as alpaqa-test

COPY --from=casadi-install /tmp/casadi-install/ /usr/local/
COPY --from=eigen-install /tmp/eigen-install/ /usr/local/
COPY --from=gtest-install /tmp/gtest-install/ /usr/local/
RUN ldconfig

FROM ghcr.io/tttapa/docker-python-build:gcc${GCC_VERSION}-py${PYTHON_VERSION} as alpaqa-python-build

RUN export DEBIAN_FRONTEND=noninteractive && \
    apt-get update -y && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:git-core/ppa && \
    apt-get update -y && \
    apt-get install -y \
        git zip unzip && \
    apt-get clean autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/ninja-build/ninja/releases/download/v1.11.0/ninja-linux.zip && \
    unzip ninja-linux.zip && \
    rm -f ninja-linux.zip && \
    mv ninja /usr/local/bin
RUN ninja --version

RUN wget https://github.com/Kitware/CMake/releases/download/v3.23.1/cmake-3.23.1-linux-x86_64.sh && \
    bash cmake-3.23.1-linux-x86_64.sh --skip-license --prefix=/usr/local && \
    rm -f cmake-3.23.1-linux-x86_64.sh
RUN cmake --version

COPY --from=casadi-install-static /tmp/casadi-install/ /usr/local/
COPY --from=eigen-install /tmp/eigen-install/ /usr/local/
COPY --from=gtest-install /tmp/gtest-install/ /usr/local/
RUN ldconfig

RUN python3 -m pip install -U build auditwheel patchelf