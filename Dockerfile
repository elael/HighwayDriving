FROM ubuntu:18.04

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

# Some C++ dev tools
RUN apt-get update \
    && apt-get -y install build-essential cppcheck valgrind \
                  git libuv1-dev wget zlib1g-dev tar libssl-dev gcc g++ cmake make \
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# Install uWebSockets
RUN git clone https://github.com/uWebSockets/uWebSockets \
    && cd uWebSockets \
    && git checkout e94b6e1 \
    && mkdir build && cd build && cmake .. \
    && make && make install \
    && cd ../.. \
    && ln -s /usr/lib64/libuWS.so /usr/lib/libuWS.so \
    && rm -r uWebSockets

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog
