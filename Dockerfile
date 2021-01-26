# NOTE:
# This Dockerfile is stored under "rplatform" directory but it will be
# moved to top level project directory during image build by RP so
# all local paths should be relative to your project's top level
# directory.
#
# NOTE2:
# Base images with tag 3.5.1_rp0.0.75 and later are built as user 'rstudio'
# which belongs to sudoers, so every command requiring root permissions should be
# preceded with 'sudo' - otherwise add layer 'USER root' to run all commands as root.
# It is recommended to not to change 'rstudio' user due to permissions issues
# within Docker container, because container's RStudio Server is run as 'rstudio'.

FROM registry.rplatform.org:5000/rocker-rstudio-tst:4.0.0_rp0.0.79

# ------ Be aware that any changes in following may cause issue with RPlatform and CBS ---------------------------

LABEL MAINTAINER="bartosz.czech@contractors.roche.com"
LABEL NAME=gdrutils
LABEL CACHE_IMAGE="registry.rplatform.org:5000/githubroche/gdrplatform/gdrutils"

# ----------------------------------------------------------------------------------------------------------------
# install system dependencies
RUN sudo apt-get update && sudo apt-get install -y \
    libssl-dev \
    libsasl2-dev \
    libxml2-dev \
    libicu-dev \
    bzip2 \
    liblzma-dev \
    libbz2-dev \
    subversion \
    curl \
    libmariadbclient-dev \
    libv8-dev \
    procps \
    systemd \
    libmagick++-dev \
    libssh2-1-dev \
    ssh \
    openssl \
    supervisor \
    passwd \
    vim

#================= copy ssh keys
COPY rplatform/ssh_keys/id_rsa /home/rstudio/.ssh/id_rsa
COPY rplatform/ssh_keys/id_rsa.pub /home/rstudio/.ssh/id_rsa.pub

#================= Add Roche certs
RUN sudo wget -O /usr/local/share/ca-certificates/Roche_G3_Root_CA.crt  http://certinfo.roche.com/rootcerts/Roche%20G3%20Root%20CA.crt 
RUN sudo update-ca-certificates 

#================= Remove openssl settings (two short keys)
#TODO: contact auth team regarding this issue
RUN sudo grep -v "^CipherString = DEFAULT@SECLEVEL=2" /etc/ssl/openssl.cnf > /tmp/openssl.fixed.cnf
RUN sudo mv /tmp/openssl.fixed.cnf /etc/ssl/openssl.cnf 

## Define your system dependencies in this Dockerfile

COPY rplatform/DESCRIPTION_dependencies.yaml /mnt/vol/rplatform/DESCRIPTION_dependencies.yaml
COPY rplatform/install_dependencies.R /mnt/vol/rplatform/install_dependencies.R
RUN R -f /mnt/vol/rplatform/install_dependencies.R

#================= install packages from specific sources
COPY rplatform/git_dependencies.yml /mnt/vol/rplatform/git_dependencies.yml
COPY rplatform/install_from_source.R /mnt/vol/rplatform/install_from_source.R
RUN R -f /mnt/vol/rplatform/install_from_source.R

RUN sudo rm -rf /mnt/vol/*

LABEL GENERATE_SINGULARITY_IMAGE=false
LABEL production=false
LABEL VERSION=0.0.0.9000
