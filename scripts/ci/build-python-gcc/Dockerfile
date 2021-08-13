ARG PYTHON_VERSION
FROM python:${PYTHON_VERSION}-buster

RUN echo 'deb http://deb.debian.org/debian testing main' >> /etc/apt/sources.list && \
    apt-get update -y && \
    apt-get install -y gcc-10 g++-10 gfortran-10 cmake bison patchelf && \
    apt-get clean autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*
