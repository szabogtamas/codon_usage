FROM szabogtamas/jupy_rocker

RUN sudo apt-get update -y && \
    sudo apt-get install -y libxt-dev && \
    sudo apt-get install -y libx11-dev

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
    readxl

RUN chmod a+rwx -R /home/rstudio

ADD ./configs/rstudio-prefs.json /home/rstudio/.config/rstudio/rstudio-prefs.json