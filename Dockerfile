# syntax=docker/dockerfile:1

# This file contains three stages: 'build', 'install', 'test'
# https://docs.docker.com/develop/develop-images/dockerfile_best-practices/#use-multi-stage-builds
#
# 'build' builds sdpb binaries. This image is heavy, because it contains all sources etc.
# 'install' contains only sdpb binaries (in /usr/local/bin/) + necessary dynamic libraries
# 'test' is based on 'install', but also copies sdpb/build/ and sdpb/test/ folders to /home/testuser/sdpb/
# - this allows to in order to run unit and integration tests by calling sdpb/test/run_all_tests.sh
# Final (lightweight) image is made from 'install' target.
#
# How to check that SDPB image works:
#
# docker build . -t sdpb
# docker run sdpb mpirun --allow-run-as-root sdpb --help
#
# How to run unit/integration tests to check that SDPB works correctly:
#
# docker build . -t sdpb-test --target test
# docker run sdpb-test ./test/run_all_tests.sh mpirun --oversubscribe
#
# Note: 'mpirun --oversubscribe' is necessary only if your environment has less than 6 CPUs available

# Elemental binaries built from https://gitlab.com/bootstrapcollaboration/elemental.git
# based on the same alpine:3.24 build as below

# Built from https://github.com/vasdommes/flint/tree/docker-main
FROM vasdommes/flint:main AS flint

FROM bootstrapcollaboration/elemental:master AS build

RUN apk add \
    autoconf \
    automake \
    binutils \
    byacc \
    cmake \
    flex \
    g++ \
    git \
    make \
    python3 \
    unzip  \
    boost-dev \
    gmp-dev \
    libarchive-dev \
    libtool \
    libxml2-dev \
    mpfr-dev \
    openblas-dev \
    openmpi \
    openmpi-dev \
    rapidjson-dev
# Include FLINT
COPY --from=flint /usr/local /usr/local
COPY --from=flint /usr/local/lib /usr/local/lib
COPY --from=flint /usr/local/include /usr/local/include

WORKDIR /usr/local/src
# Build MPSolve (takes ~1 minute)
# TODO: create separate MPSolve image?
# -Wno-error=implicit-function-declaration: MPSolve calls undeclared yyparse(), which is an error since GCC 14.
RUN git clone https://github.com/robol/MPSolve.git --depth=1 && \
    cd MPSolve && \
    ./autogen.sh && \
    CC=mpicc CXX=mpicxx CFLAGS="-g -O2 -Wno-error=implicit-function-declaration" \
        ./configure --disable-dependency-tracking --disable-examples --disable-ui --disable-graphical-debugger --disable-documentation && \
    make && \
    make check && \
    make install && \
    echo "simple mpsolve test:" && \
    mpsolve -p "x^4 - 6*x^9 + 6/7*x + 5"

WORKDIR /usr/local/src/sdpb
# Build SDPB from current sources, print build/config.log in configuration failed
COPY . .
RUN (./waf configure --elemental-dir=/usr/local --flint-dir=/usr/local --prefix=/usr/local \
        || (cat build/config.log && exit 1) \
    ) && \
    python3 ./waf && \
    python3 ./waf install

# Take only sdpb binaries + load necessary dynamic libraries
# libboost_stacktrace_addr2line is not packaged separately, only in boost1.84-libs (all Boost libraries).
# binutils provides addr2line, which prints function names and source locations in stacktraces.
FROM alpine:3.24 AS install
RUN apk add \
    binutils \
    boost1.84-libs \
    gmp \
    libarchive \
    libgmpxx \
    libstdc++ \
    libxml2 \
    mpfr \
    openblas \
    openmpi \
    openssh
COPY --from=build /usr/local/bin /usr/local/bin
COPY --from=build /usr/local/lib /usr/local/lib

# Separate test target, see
# https://docs.docker.com/language/java/run-tests/#multi-stage-dockerfile-for-testing
# Contains /home/testuser/sdpb/build and /home/testuser/sdpb/test folders,
# which is sufficient to run tests as shown below:
#
# docker build . -t sdpb-test --target test
# docker run sdpb-test ./test/run_all_tests.sh mpirun --oversubscribe
FROM install AS test
# Create testuser to run Docker non-root (and avoid 'mpirun --allow-run-as-root' warning)
RUN addgroup --gid 10000 testgroup && \
    adduser --disabled-password --uid 10000 --ingroup testgroup --shell /bin/sh testuser
WORKDIR /home/testuser/sdpb
# Tests write to build/test/out
COPY --from=build --chown=testuser:testgroup /usr/local/src/sdpb/build build
COPY --from=build --chown=testuser:testgroup /usr/local/src/sdpb/test test
USER testuser:testgroup

# Resulting image
FROM install
# TODO best practices suggest to run containter as non-root.
# https://github.com/dnaprawa/dockerfile-best-practices#run-as-a-non-root-user
# But this requires some extra work with permissions when mounting folders for docker run.