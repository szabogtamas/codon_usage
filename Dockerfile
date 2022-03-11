FROM szabogtamas/jupy_rocker

RUN sudo apt-get update -y && \
    sudo apt-get install -y libxml2-dev && \
    sudo apt-get install -y libxt-dev && \
    sudo apt-get install -y libx11-dev && \
    sudo apt-get install -y libz-dev && \
    sudo apt-get install -y libbz2-dev && \
    sudo apt-get install -y liblzma-dev && \
    sudo apt-get install -y zlib1g-dev && \
    sudo apt-get install -y libglpk-dev && \
    sudo apt-get install -y libcairo2-dev

RUN pip3 install jupytext && \
    pip3 install numpy && \
    pip3 install pandas && \
    pip3 install matplotlib && \
    pip3 install seaborn && \
    pip3 install biopython

RUN install2.r --error \
    --deps TRUE \
    devtools \
    plotly \
    heatmaply \
    RColorBrewer \
    ggsci \
    ggridges \
    readxl \
    msigdbr

RUN R -e "BiocManager::install('EnsDb.Hsapiens.v86')"
RUN R -e "BiocManager::install('EnsDb.Mmusculus.v79')"
RUN R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')"
RUN R -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')"
RUN R -e "BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')"
RUN R -e "BiocManager::install('TxDb.Mmusculus.UCSC.mm10.knownGene')"
RUN R -e "BiocManager::install('BSgenome.Drerio.UCSC.danRer10')"
RUN R -e "BiocManager::install('TxDb.Drerio.UCSC.danRer10.refGene')"
RUN R -e "BiocManager::install('BSgenome.Dmelanogaster.UCSC.dm6')"
RUN R -e "BiocManager::install('TxDb.Dmelanogaster.UCSC.dm6.ensGene')"
RUN R -e "BiocManager::install('org.Dm.eg.db')"
RUN R -e "BiocManager::install('BSgenome.Scerevisiae.UCSC.sacCer3')"
RUN R -e "BiocManager::install('TxDb.Scerevisiae.UCSC.sacCer3.sgdGene')"
RUN R -e "BiocManager::install('org.Sc.sgd.db')"
RUN R -e "BiocManager::install('GenomicFeatures')"
RUN R -e "BiocManager::install('AnnotationHub')"
RUN R -e "BiocManager::install('clusterProfiler')"

RUN chmod a+rwx -R /home/rstudio

ADD ./configs/rstudio-prefs.json /home/rstudio/.config/rstudio/rstudio-prefs.json
