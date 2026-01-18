FROM ubuntu:16.04

MAINTAINER rosario.distefano.ict@gmail.com
ENV DEBIAN_FRONTEND noninteractive

LABEl edu.osumc.dept="Department of Cancer Biology and Genetics - The Ohio State University" \
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

ADD tools /opt/
ADD v2 /opt/

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
	libreadline-gplv2-dev \
	libsqlite3-dev \
	libssl-dev \
	libtclap-dev \
	libxml2 \
	libxml2-dev \
	nano \
	nginx \
	pkg-config \
	python-dev \
	python-mysqldb \
	python-pip \
	python-setuptools \
	python-virtualenv \
	python2.7 \
	rabbitmq-server \
	software-properties-common \
	supervisor \
	tk-dev \
	wget \
	curl \
	bzip2 \
	ca-certificates \
	&& rm -rf /var/lib/apt/lists/*

# RUN pip install --upgrade pip
RUN pip install "numpy==1.16.6"

# Setup Miniconda for DMISO (Python 3.6)
ENV CONDA_DIR=/opt/conda
RUN curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh -o /tmp/miniconda.sh \
	&& bash /tmp/miniconda.sh -b -p ${CONDA_DIR} \
	&& rm /tmp/miniconda.sh \
	&& ${CONDA_DIR}/bin/conda clean -a -y
ENV PATH="${CONDA_DIR}/bin:${PATH}"
RUN conda create -y -n dmiso python=3.6 \
	&& conda run -n dmiso pip install --no-cache-dir "tensorflow==1.15.0" "keras==2.3.1" "numpy" \
	&& conda run -n dmiso python -c "import tensorflow, keras; print('dmiso env ok')" \
	&& conda clean -a -y
ENV DMISO_HOME=/opt/DMISO/DMISO-main
ENV PATH="${PATH}:${DMISO_HOME}"
RUN printf '%s\n' "#!/bin/sh" "exec conda run -n dmiso python ${DMISO_HOME}/dmiso.py \"$@\"" > /usr/local/bin/dmiso \
	&& chmod +x /usr/local/bin/dmiso

#########################
## R                   ##
#########################

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 \
	&& add-apt-repository 'deb [arch=amd64,i386] https://ftp.ussg.iu.edu/CRAN/bin/linux/ubuntu xenial/' \
	&& apt-get clean \
	&& apt-get update \
	&& apt-get install -y --no-install-recommends \
	curl \
	libcairo2-dev \
	libcurl4-openssl-dev \
	librsvg2-dev \
	libssl-dev \
	libv8-3.14-dev \
	libwebp-dev \
	r-base \
	r-base-dev \
	r-cran-bitops \
	r-cran-igraph \
	r-cran-rcurl \
	r-cran-xml \
	r-omegahat-xmlrpc

RUN cd /opt/R \
	&& /usr/bin/Rscript install_r_packages.R 

# Setup flask application
RUN mkdir -p /app
RUN /usr/bin/pip install -r /opt/requirements.txt

# Setup Vienna-rna
RUN cd /opt \
	&& dpkg -i viennarna_2.4.11-1_amd64.deb


# Setup Spatt
RUN cd /opt/spatt \
	&& mkdir build \
	&& cd build \
	&& cmake .. \
	&& make \
	&& make install

# Setup Phast and CLAPACK
RUN cd /opt/CLAPACK-3.2.1 && cp make.inc.example make.inc \
	&& cd /opt/CLAPACK-3.2.1 && make f2clib \
	&& cd /opt/CLAPACK-3.2.1 && make blaslib \
	&& cd /opt/CLAPACK-3.2.1 && make lib
	
RUN cd /opt/ \
	&& dpkg -i phast.v1_4.x86_64.deb

# miRanda setup
RUN cd /opt/miRanda \
	&& ./configure \
	&& make \
	&& make install

# Setup miRmap
RUN cd /opt/miRmap/libs \
	&& mv lib-archlinux-x86_64/ default
	
# RNAhybrid-2.1.2 setup
RUN cd /opt/RNAhybrid \
	&& make clean \
	&& ./configure \
        && make \
        && make install

###################################
## miRNA target prediction tools ##
###################################

#########################
## miRanda v3.3        ##
#########################
ENV PATH="${PATH}:/opt/miRanda/bin"

#########################
## miRmap v1.1         ##
#########################
ENV PATH="${PATH}:/opt/miRmap/scripts"
ENV PYTHONPATH="${PYTHONPATH}:/opt/miRmap/src/"
ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/opt/miRmap/libs/default"

###############################
## Setting PERL5LIB for PITA ##
###############################
ENV PERL5LIB="${PERL5LIB}:/opt/PITA64bit/lib"

#########################
## RNAhybrid 2.1.2     ##
#########################
ENV PATH="${PATH}:/opt/RNAhybrid/src"

##########################
## PITA (v.6 31-Aug-08) ##
##########################
ENV PATH="${PATH}:/opt/PITA64bit"

#########################
## Spatt 2.0           ##
#########################
ENV PATH="${PATH}:/opt/spatt/build/src"

#########################
## TargetScan          ##
#########################
ENV PATH="${PATH}:/opt/TargetScan/TargetScan_70"
ENV PATH="${PATH}:/opt/TargetScan/TargetScan7_BL_PCT"
ENV PATH="${PATH}:/opt/TargetScan/TargetScan7_context_scores"

WORKDIR /app

COPY app /app

RUN mkdir -p /opt/TargetScan/Datasets \
	&& mkdir /input


# Setup nginx
RUN rm /etc/nginx/nginx.conf
COPY nginx.conf /etc/nginx/nginx.conf
RUN rm /etc/nginx/sites-enabled/default
COPY app.conf /etc/nginx/sites-available/
RUN ln -s /etc/nginx/sites-available/app.conf /etc/nginx/sites-enabled/app.conf

# Setup supervisord
RUN mkdir -p /var/log/supervisor
COPY app_supervisord.conf /etc/supervisor/conf.d/app_supervisord.conf
COPY app_gunicorn.conf /etc/supervisor/conf.d/app_gunicorn.conf

ADD rabbitmq.sh /opt/
RUN chmod a+x /opt/rabbitmq.sh

ADD kill_zombies.sh /opt/
RUN chmod a+x /opt/kill_zombies.sh

# Start processes
CMD ["/usr/bin/supervisord"]
