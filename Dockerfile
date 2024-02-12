ARG BASE_IMAGE=bioconductor/bioconductor_docker:devel
FROM ${BASE_IMAGE}

# ------ Be aware that any changes in following may cause issue with RPlatform and CBS

LABEL MAINTAINER="bartosz.czech@contractors.roche.com"
LABEL NAME=gdrutils
LABEL GENERATE_SINGULARITY_IMAGE=false
LABEL production=false
LABEL VERSION=0.0.0.9100
#LABEL CACHE_IMAGE="registry.rplatform.org:5000/githubroche/gdrplatform/gdrutils"

# temporary fix
# GitHub token for downloading private dependencies
# Need to be defined after FROM as it flushes ARGs
ARG GITHUB_TOKEN

#================= Install dependencies
RUN apt-get update && apt-get -y install git
RUN mkdir -p /mnt/vol
COPY rplatform/dependencies.yaml rplatform/.github_access_token.txt* /mnt/vol
RUN echo "$GITHUB_TOKEN" >> /mnt/vol/.github_access_token.txt
RUN Rscript -e 'BiocManager::install(c("gDRstyle", "gDRtestData", "gDRutils", "BiocStyle", "qs"))'
RUN Rscript -e "gDRstyle::installAllDeps()"

#================= Check & build package
COPY ./ /tmp/gDRutils/
RUN Rscript -e "gDRstyle::installLocalPackage('/tmp/gDRutils')"

#================= Clean up
#RUN sudo rm -rf /mnt/vol/* /tmp/gDRutils/
