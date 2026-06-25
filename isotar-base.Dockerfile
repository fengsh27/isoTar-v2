FROM ubuntu:16.04

LABEL maintainer="rosario.distefano.ict@gmail.com"
ENV DEBIAN_FRONTEND=noninteractive

LABEL edu.osumc.dept="Department of Cancer Biology and Genetics - The Ohio State University" \
      edu.osumc.version="1.2" \
      edu.osumc.is-final="" \
      edu.osumc.released="March 27, 2020"

# a few minor docker-specific tweaks
# see https://github.com/docker/docker/blob/9a9fc01af8fb5d98b8eec0740716226fadb3735c/contrib/mkimage/debootstrap
RUN set -xe \
	\
# https://github.com/docker/docker/blob/9a9fc01af8fb5d98b8eec0740716226fadb3735c/contrib/mkimage/debootstrap#L40-L48
	&& echo '#!/bin/sh' > /usr/sbin/policy-rc.d \
	&& echo 'exit 101' >> /usr/sbin/policy-rc.d \
	&& chmod +x /usr/sbin/policy-rc.d \
	\
# https://github.com/docker/docker/blob/9a9fc01af8fb5d98b8eec0740716226fadb3735c/contrib/mkimage/debootstrap#L54-L56
	&& dpkg-divert --local --rename --add /sbin/initctl \
	&& cp -a /usr/sbin/policy-rc.d /sbin/initctl \
	&& sed -i 's/^exit.*/exit 0/' /sbin/initctl \
	\
# https://github.com/docker/docker/blob/9a9fc01af8fb5d98b8eec0740716226fadb3735c/contrib/mkimage/debootstrap#L71-L78
	&& echo 'force-unsafe-io' > /etc/dpkg/dpkg.cfg.d/docker-apt-speedup \
	\
# https://github.com/docker/docker/blob/9a9fc01af8fb5d98b8eec0740716226fadb3735c/contrib/mkimage/debootstrap#L85-L105
	&& echo 'DPkg::Post-Invoke { "rm -f /var/cache/apt/archives/*.deb /var/cache/apt/archives/partial/*.deb /var/cache/apt/*.bin || true"; };' > /etc/apt/apt.conf.d/docker-clean \
	&& echo 'APT::Update::Post-Invoke { "rm -f /var/cache/apt/archives/*.deb /var/cache/apt/archives/partial/*.deb /var/cache/apt/*.bin || true"; };' >> /etc/apt/apt.conf.d/docker-clean \
	&& echo 'Dir::Cache::pkgcache ""; Dir::Cache::srcpkgcache "";' >> /etc/apt/apt.conf.d/docker-clean \
	\
# https://github.com/docker/docker/blob/9a9fc01af8fb5d98b8eec0740716226fadb3735c/contrib/mkimage/debootstrap#L109-L115
	&& echo 'Acquire::Languages "none";' > /etc/apt/apt.conf.d/docker-no-languages \
	\
# https://github.com/docker/docker/blob/9a9fc01af8fb5d98b8eec0740716226fadb3735c/contrib/mkimage/debootstrap#L118-L130
	&& echo 'Acquire::GzipIndexes "true"; Acquire::CompressionTypes::Order:: "gz";' > /etc/apt/apt.conf.d/docker-gzip-indexes \
	\
# https://github.com/docker/docker/blob/9a9fc01af8fb5d98b8eec0740716226fadb3735c/contrib/mkimage/debootstrap#L134-L151
	&& echo 'Apt::AutoRemove::SuggestsImportant "false";' > /etc/apt/apt.conf.d/docker-autoremove-suggests

# delete all the apt list files since they're big and get stale quickly
RUN rm -rf /var/lib/apt/lists/*
# this forces "apt-get update" in dependent images, which is also good
# (see also https://bugs.launchpad.net/cloud-images/+bug/1699913)

# enable the universe
RUN sed -i 's/^#\s*\(deb.*universe\)$/\1/g' /etc/apt/sources.list

# make systemd-detect-virt return "docker"
# See: https://github.com/systemd/systemd/blob/aa0c34279ee40bce2f9681b496922dedbadfca19/src/basic/virt.c#L434
RUN mkdir -p /run/systemd && echo 'docker' > /run/systemd/container


RUN apt-get clean && apt-get update && apt-get install -y --no-install-recommends \
	apt-transport-https \
	apt-utils \
	autoconf \
	automake \
	autotools-dev \
	bioperl \
	build-essential \
	checkinstall \
	cmake \
	dh-autoreconf \
	g++ \
	gcc \
	graphviz \
	gridengine-drmaa1.0 \
	gsl-bin \
	gunicorn \
	hdf5-tools \
	libbz2-dev \
	libc6-dev \
	libg2-dev \
	libgdbm-dev \
	libgraphviz-dev \
	libhdf5-serial-dev \
	libncursesw5-dev \
	libreadline-dev \
	libsqlite3-dev \
	libssl-dev \
	libmysqlclient-dev \
	libtclap-dev \
	libxml2 \
	libxml2-dev \
	nano \
	nginx \
	pkg-config \
	python3-dev \
	python3-mysqldb \
	python3-pip \
	python3-setuptools \
	python3-venv \
	python-dev \
	python-pip \
	rabbitmq-server \
	software-properties-common \
	supervisor \
	tk-dev \
	wget \
	curl \
	bzip2 \
	ca-certificates \
	&& rm -rf /var/lib/apt/lists/*

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
	libssl-dev \
	zlib1g-dev \
	libbz2-dev \
	libreadline-dev \
	libsqlite3-dev \
	libffi-dev \
	libncursesw5-dev \
	libgdbm-dev \
	liblzma-dev \
	libuuid1 \
	libuuidm-ocaml-dev \
	ca-certificates \
	make \
	gcc \
	wget \
	&& rm -rf /var/lib/apt/lists/*

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
	curl \
	libcairo2-dev \
	libcurl4-openssl-dev \
	librsvg2-dev \
	libssl-dev \
	libv8-3.14-dev \
	libwebp-dev \
	r-base \
	&& rm -rf /var/lib/apt/lists/*

RUN wget -q https://www.python.org/ftp/python/3.6.13/Python-3.6.13.tgz -O /tmp/Python-3.6.13.tgz \
	&& tar -xzf /tmp/Python-3.6.13.tgz -C /tmp \
	&& cd /tmp/Python-3.6.13 \
	&& ./configure --with-ensurepip=install \
	&& make -j$(nproc) \
	&& make altinstall \
	&& cd / \
	&& rm -rf /tmp/Python-3.6.13 /tmp/Python-3.6.13.tgz

# --- miRmap 2 toolchain: Python 3.11 (needs OpenSSL >= 1.1.1, which ubuntu:16.04
# does not ship). Built with altinstall so it sits alongside python2.7/3.5/3.6
# without clobbering the system python3. Invoke explicitly as `python3.11`. The
# DMISO (python3.6) and miRmap-1 (python2.7) paths are untouched; this is purely
# additive. Required build headers (libffi/zlib/bz2/sqlite/readline/lzma/ssl)
# are already installed by the apt blocks above. ---

# OpenSSL 1.1.1 into a private prefix (16.04's system OpenSSL is 1.0.2, too old
# for Python 3.10+ to build a working `ssl` module -> pip can't reach PyPI).
RUN wget -q https://www.openssl.org/source/openssl-1.1.1w.tar.gz -O /tmp/openssl-1.1.1w.tar.gz \
	&& tar -xzf /tmp/openssl-1.1.1w.tar.gz -C /tmp \
	&& cd /tmp/openssl-1.1.1w \
	&& ./config --prefix=/usr/local/openssl-1.1.1 --openssldir=/usr/local/openssl-1.1.1 shared zlib \
	&& make -j$(nproc) \
	&& make install_sw \
	&& cd / \
	&& rm -rf /tmp/openssl-1.1.1w /tmp/openssl-1.1.1w.tar.gz

# Python 3.11 built against that OpenSSL. --with-openssl-rpath=auto bakes the lib
# path into the binary so no LD_LIBRARY_PATH is needed at runtime.
RUN wget -q https://www.python.org/ftp/python/3.11.9/Python-3.11.9.tgz -O /tmp/Python-3.11.9.tgz \
	&& tar -xzf /tmp/Python-3.11.9.tgz -C /tmp \
	&& cd /tmp/Python-3.11.9 \
	&& ./configure --with-openssl=/usr/local/openssl-1.1.1 \
	               --with-openssl-rpath=auto \
	               --with-ensurepip=install \
	&& make -j$(nproc) \
	&& make altinstall \
	&& cd / \
	&& rm -rf /tmp/Python-3.11.9 /tmp/Python-3.11.9.tgz

# Verify ssl actually compiled (fails the build early if not), then miRmap 2.
# miRmap declares bare `ViennaRNA`/`dendropy` (no version pins), so pip would
# otherwise grab the newest ViennaRNA, whose wheels target a glibc newer than
# 16.04's (glibc 2.23) -> no matching wheel -> source build -> fail. Pin to a
# release that still ships a manylinux_2_17 (glibc 2.17) cp311 wheel and force a
# binary install (--only-binary :all: fails loudly instead of compiling). This
# also satisfies miRmap's unpinned dep so `pip install mirmap` keeps this version.
RUN python3.11 -m pip install --no-cache-dir --upgrade pip \
	&& python3.11 -c "import ssl; print('python3.11 OpenSSL:', ssl.OPENSSL_VERSION)" \
	&& python3.11 -m pip install --no-cache-dir --only-binary :all: "ViennaRNA==2.7.0" dendropy \
	&& python3.11 -m pip install --no-cache-dir mirmap \
	&& python3.11 -c "import mirmap, RNA, dendropy; print('miRmap 2 import OK')"
