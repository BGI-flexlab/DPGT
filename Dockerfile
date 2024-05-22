FROM archieyoung/spark:2.4.5

RUN apt-get -y update
RUN apt-get install -y \
    autoconf automake make gcc g++ \
    perl zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-gnutls-dev libssl-dev \
    cmake libboost-all-dev libjemalloc-dev

RUN mkdir /opt/dpgt-src
COPY src /opt/dpgt-src/src
# COPY test /opt/dpgt-src/test

WORKDIR /opt/dpgt-src
# Clean the vendor library to prevent the contamination of building DPGT in this image.
RUN cd /opt/dpgt-src/src/main/native/vendor/htslib && make clean
RUN cd /opt/dpgt-src/src/main/native/vendor/libdeflate && make clean

RUN mkdir build && \
    cd build && \
    cmake -DJAVA_INCLUDE_PATH=/usr/local/openjdk-8/include ../src/main/native && \
    make -j
RUN mkdir /opt/dpgt && cp -r build/lib /opt/dpgt/lib

# remove src
RUN rm -r /opt/dpgt-src

WORKDIR /
# copy dpgt jar
COPY target/dpgt-1.3.2.0.jar /opt/dpgt/dpgt-1.3.2.0.jar

# set environment variable to load dpgt cpp library and jemalloc library
ENV LD_LIBRARY_PATH=/opt/dpgt/lib:$LD_LIBRARY_PATH
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so

# create input and output dir
RUN mkdir /input
RUN mkdir /output

RUN mkdir /workdir
RUN chmod 777 /workdir

WORKDIR /workdir
