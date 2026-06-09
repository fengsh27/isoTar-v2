FROM frankfeng78/isotar-v2-base:0.2.2

LABEL maintainer="rosario.distefano.ict@gmail.com"
ENV DEBIAN_FRONTEND=noninteractive

LABEL edu.osumc.dept="Department of Cancer Biology and Genetics - The Ohio State University" \
      edu.osumc.version="1.2" \
      edu.osumc.is-final="" \
      edu.osumc.released="March 27, 2020"

# Application version. Pass at build time: `--build-arg VERSION=$(cat VERSION)`.
# Exposed at runtime via ISOTAR_VERSION (read by app_v1/version.py and the
# /api/v1/version endpoint) and as a standard OCI image label.
ARG VERSION=unknown
ENV ISOTAR_VERSION=${VERSION}
LABEL org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.title="isotar-v2-backend"

ADD tools /opt/
COPY v2/*.py /opt/v2/
COPY v2/opt/reference_files /opt/reference_files
# Per-species miRNA metadata: the human file is named without a prefix
# (mature_pre_mirna_ext.json); every other species ships as
# <species>_mature_pre_mirna_ext.json (e.g. mmu_, dre_, cel_, ...).
# v2/mirna_processing.py picks the right file at runtime from the miRBase
# prefix in the miRNA id, looking under ISOTAR_RESOURCES_DIR.
COPY v2/opt/resources /opt/resources
ENV ISOTAR_RESOURCES_DIR=/opt/resources
# Fail the build if the per-species files didn't make it into the image
# (e.g. accidental rename / removal under v2/opt/resources). The human file
# is required; we additionally require at least one <species>_*.json so a
# regression that silently drops every non-human file is caught here, not
# at job time.
RUN test -f /opt/resources/mature_pre_mirna_ext.json \
    && ls /opt/resources/*_mature_pre_mirna_ext.json >/dev/null \
    && echo "miRNA metadata present:" \
    && ls /opt/resources/*mature_pre_mirna_ext.json


RUN python3.6 -m pip install --no-cache-dir --upgrade "pip==21.3.1"

RUN ln -sf /usr/local/bin/python3.6 /usr/local/bin/python3

ENV DMISO_HOME=/opt/DMISO/DMISO-main
ENV PATH="${PATH}:${DMISO_HOME}"
RUN python3.6 -m pip install --no-cache-dir "tensorflow==1.15.0" "keras==2.3.1" "numpy" "h5py==2.10.0" \
	&& python3.6 -c "import tensorflow, keras; print('dmiso env ok')" \
	&& printf '%s\n' "#!/bin/sh" "exec /usr/local/bin/python3.6 ${DMISO_HOME}/dmiso.py \"$@\"" > /usr/local/bin/dmiso \
	&& chmod +x /usr/local/bin/dmiso

#########################
## R                   ##
#########################

RUN printf '%s\n' \
	"CXX=g++ -std=gnu++11" \
	"CXXFLAGS=-std=gnu++11" \
	"CXX11=g++ -std=gnu++11" \
	"CXX11FLAGS=-std=gnu++11" \
	> /etc/R/Makevars

RUN cd /opt/R \
	&& /usr/bin/Rscript -e "options(repos='https://cloud.r-project.org'); source('install_r_packages.R')" 

# Setup flask application
RUN mkdir -p /app
RUN python2.7 -m pip install --no-cache-dir --upgrade "pip==20.3.4" "setuptools==44.1.1" \
	&& python2.7 -m pip install --no-cache-dir "numpy==1.16.6" \
	&& python2.7 -m pip install --no-cache-dir "dendropy==4.3.0" \
	&& python2.7 -m pip install --no-cache-dir -r /opt/requirements.txt

# Setup app_v1
COPY app_v1 /app_v1
# Bundle VERSION as a fallback for /api/v1/version when the build arg is omitted.
COPY VERSION /app_v1/VERSION
RUN python3.6 -m pip install --no-cache-dir -r /app_v1/requirements.txt \
    && mkdir -p /app/logs/celery /app/logs/gunicorn

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
ENV PYTHONPATH=/opt/miRmap/src/
ENV LD_LIBRARY_PATH=/opt/miRmap/libs/default

###############################
## Setting PERL5LIB for PITA ##
###############################
ENV PERL5LIB=/opt/PITA64bit/lib

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
