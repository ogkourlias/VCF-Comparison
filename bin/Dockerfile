################## BASE IMAGE ######################

FROM ubuntu:23.04

################## METADATA ######################
LABEL base_image="ubuntu:23.04"
LABEL version="1.0.0"
LABEL software="VCF-Comparison Pipeline"

################## MAINTAINER ######################

LABEL maintainer="Orfeas Gkourlias <o.gkourlias@umcg.nl>"


################## INSTALLATION ######################
ADD . /tmp/init
WORKDIR /tmp/init

ENV SHELL=/bin/bash
ENV LC_ALL=C
ENV PIP_BREAK_SYSTEM_PACKAGES=1
ENV LANG=C.UTF-8
ENV TZ=Europe
ENV DEBIAN_FRONTEND=noninteractive

# Getting base apps & languages.
RUN apt-get update -y \
    && apt-get upgrade -y \
    && apt-get install -y \
        # Getters & VSC.
        wget=1.21.3-1ubuntu1 \
        curl \
        # Languages.
        # C++
        g++=4:12.2.0-3ubuntu1 \
        # Python
        python-is-python3 \
        pip \
        # Other compilers, builders & compressors & libs.
        libbz2-dev \
        liblzma-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        make=4.3-4.1build1
        # Tools
        # Other

# Non-apt Tools

# Python Packages
RUN pip install pytabix pandas argparse numpy scipy
    
# HTSlib 1.18 Install
RUN cd /tmp/init \
    && wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 \
    && tar -xvjf htslib-1.18.tar.bz2 \
    && cd htslib-1.18 \
    && ./configure --prefix=/usr/ --enable-libcurl \
    && make \
    && make install


################## CLEANUP ######################
# Apt cleanup.
RUN apt-get clean \
    && apt-get autoremove -y

# Build files cleanup
RUN cd /
RUN rm -rf /tmp/init