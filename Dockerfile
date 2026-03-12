FROM frankfeng78/isotar-v2-base:0.2.2

MAINTAINER rosario.distefano.ict@gmail.com
ENV DEBIAN_FRONTEND noninteractive

LABEl edu.osumc.dept="Department of Cancer Biology and Genetics - The Ohio State University" \
      edu.osumc.version="1.2" \
      edu.osumc.is-final="" \
      edu.osumc.released="March 27, 2020"

ADD tools /opt/
COPY v2/*.py /opt/v2/
COPY v2/opt/human /opt/human
COPY v2/opt/resources /opt/resources


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
RUN python3.6 -m pip install --no-cache-dir -r /app_v1/requirements.txt \
    && mkdir -p /app/logs/celery

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
ENV PYTHONPATH="${PYTHONPATH:-}:/opt/miRmap/src/"
ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:/opt/miRmap/libs/default"

###############################
## Setting PERL5LIB for PITA ##
###############################
ENV PERL5LIB="${PERL5LIB:-}:/opt/PITA64bit/lib"

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
