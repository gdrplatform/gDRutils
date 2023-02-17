ARG BASE_IMAGE=marcinkam/gdrshiny:0.11

FROM ${BASE_IMAGE}

# Need to be defined after FROM as it flushes ARGs
ARG GITHUB_TOKEN

# ------ Be aware that any changes in following may cause issue with RPlatform and CBS

LABEL MAINTAINER="bartosz.czech@contractors.roche.com"
LABEL NAME=gdrutils
LABEL GENERATE_SINGULARITY_IMAGE=false
LABEL production=false
LABEL VERSION=0.0.0.9100
LABEL CACHE_IMAGE="registry.rplatform.org:5000/githubroche/gdrplatform/gdrutils"

#================= Install dependencies
RUN mkdir -p /mnt/vol
COPY rplatform/dependencies.yaml rplatform/.github_access_token.txt* /mnt/vol
RUN echo "$GITHUB_TOKEN" >> /mnt/vol/.github_access_token.txt
RUN Rscript -e "gDRstyle::installAllDeps()"

#================= Check & build package
COPY gDRutils/ /mnt/vol/gDRutils/
RUN Rscript -e "gDRstyle::installLocalPackage('/mnt/vol/gDRutils')"

#================= Clean up
RUN sudo rm -rf /mnt/vol/rplatform/.github_access_token.txt
