FROM tttapa/alpaqa-python-build:py3.9 AS alpaqa-docs

RUN export DEBIAN_FRONTEND=noninteractive && \
    apt-get update -y && \
    apt-get install -y \
        bison flex graphviz libjson-perl libperlio-gzip-perl perl && \
    apt-get clean autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

COPY install-lcov.sh .
RUN bash ./install-lcov.sh /usr/local

COPY install-doxygen.sh .
RUN bash ./install-doxygen.sh /usr/local
