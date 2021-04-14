FROM rocker/verse:latest

RUN apt-get update
RUN apt-get -y install libudunits2-dev

RUN Rscript -e "update.packages(ask=FALSE)"

RUN Rscript -e "install.packages('reshape2', quiet = TRUE)"
RUN Rscript -e "install.packages('circlize', quiet = TRUE)"
RUN Rscript -e "install.packages('BiocManager', quiet = TRUE)"

RUN Rscript -e "options(repos = BiocManager::repositories())"

RUN Rscript -e "BiocManager::install('GenomicRanges', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('ChIPseeker', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('SummarizedExperiment', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('bamsignals', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg18.knownGene', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('Gviz', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('ComplexHeatmap', quiet = TRUE)"
