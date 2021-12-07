FROM szabogtamas/jupy_rocker

RUN sudo apt-get update -y && \
    sudo apt-get install -y libxml2-dev && \
    sudo apt-get install -y libxt-dev && \
    sudo apt-get install -y libx11-dev && \
    sudo apt-get install -y libz-dev && \
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
    readxl

RUN R -e "BiocManager::install('clusterProfiler')"

RUN chmod a+rwx -R /home/rstudio

ADD ./configs/rstudio-prefs.json /home/rstudio/.config/rstudio/rstudio-prefs.json