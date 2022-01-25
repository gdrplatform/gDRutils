FROM arkadiuszgladki/gdr_shiny:0.08

# ------ Be aware that any changes in following may cause issue with RPlatform and CBS

LABEL MAINTAINER="bartosz.czech@contractors.roche.com"
LABEL NAME=gdrutils-gh
LABEL GENERATE_SINGULARITY_IMAGE=false
LABEL production=false
LABEL VERSION=0.0.0.9100
LABEL CACHE_IMAGE="registry.rplatform.org:5000/github/gdrplatform/gdrutils-gh"

# temporary fix
# GitHub token for downloading private dependencies
ARG GITHUB_TOKEN

#================= Install dependencies
RUN mkdir -p /mnt/vol
COPY rplatform/dependencies.yaml rplatform/.github_access_token.txt* /mnt/vol
COPY rplatform/install_all_deps.R /mnt/vol/install_all_deps.R
RUN R -f /mnt/vol/install_all_deps.R

#================= Check & build package
COPY gDRutils/ /tmp/gDRutils/
COPY rplatform/install_repo.R /mnt/vol
RUN R -f /mnt/vol/install_repo.R 

#================= Clean up
RUN sudo rm -rf /mnt/vol/* /tmp/gDRutils/
