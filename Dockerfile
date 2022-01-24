FROM dockerreg.kubemeainfra.science.roche.com/cbs/githubroche/gdrplatform/gdr_shiny:master-c274a040-build-8

# ------ Be aware that any changes in following may cause issue with RPlatform and CBS

LABEL MAINTAINER="bartosz.czech@contractors.roche.com"
LABEL NAME=gdrutils-gh
LABEL GENERATE_SINGULARITY_IMAGE=false
LABEL production=false
LABEL VERSION=0.0.0.9100
LABEL CACHE_IMAGE="registry.rplatform.org:5000/github/gdrplatform/gdrutils-gh"

#================= Install dependencies
COPY rplatform/dependencies.yaml rplatform/.github_access_token.txt* /mnt/vol/
COPY rplatform/install_all_deps.R /mnt/vol/install_all_deps.R
RUN R -f /mnt/vol/install_all_deps.R

#================= Check & build package
COPY gDRutils/ /tmp/gDRutils/
COPY rplatform/install_repo.R /mnt/vol/
RUN R -f /mnt/vol/install_repo.R 

#================= Clean up
RUN sudo rm -rf /mnt/vol/* /tmp/gDRutils/
