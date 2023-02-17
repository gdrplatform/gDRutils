ARG BASE_IMAGE=marcinkam/gdrshiny:0.11
FROM ${BASE_IMAGE}

# ------ Be aware that any changes in following may cause issue with RPlatform and CBS

LABEL MAINTAINER="bartosz.czech@contractors.roche.com"
LABEL NAME=gdrutils
LABEL GENERATE_SINGULARITY_IMAGE=false
LABEL production=false
LABEL VERSION=0.0.0.9100
LABEL CACHE_IMAGE="registry.rplatform.org:5000/githubroche/gdrplatform/gdrutils"

# temporary fix
# GitHub token for downloading private dependencies
ARG GITHUB_TOKEN

#================= Install dependencies
RUN mkdir -p /mnt/vol
RUN if [[ ! -f "rplatform/.github_access_token.txt" ]] ; then echo "$GITHUB_TOKEN" > rplatform/.github_access_token.txt ; fi
COPY rplatform/dependencies.yaml rplatform/.github_access_token.txt* /mnt/vol
RUN Rscript -e "gDRstyle::installAllDeps()"

#================= Check & build package
COPY gDRutils/ /tmp/gDRutils/
RUN Rscript -e "gDRstyle::installLocalPackage('/tmp/gDRutils')"

#================= Clean up
RUN sudo rm -rf /mnt/vol/* /tmp/gDRutils/
