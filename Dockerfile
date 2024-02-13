ARG BASE_IMAGE=bioconductor/bioconductor_docker:devel
FROM ${BASE_IMAGE}

# temporary fix
# GitHub token for downloading private dependencies
# Need to be defined after FROM as it flushes ARGs
ARG GITHUB_TOKEN

#================= Install dependencies
RUN mkdir -p /mnt/vol
COPY rplatform/dependencies.yaml rplatform/.github_access_token.txt* /mnt/vol
RUN echo "$GITHUB_TOKEN" >> /mnt/vol/.github_access_token.txt
RUN Rscript -e 'BiocManager::install(c("gDRstyle", "gDRtestData", "gDRutils", "BiocStyle", "qs"))'
RUN Rscript -e "gDRstyle::installAllDeps()"

#================= Check & build package
COPY ./ /tmp/gDRutils/
RUN Rscript -e "gDRstyle::installLocalPackage('/tmp/gDRutils')"

#================= Clean up
RUN sudo rm -rf /mnt/vol/* /tmp/gDRutils/
