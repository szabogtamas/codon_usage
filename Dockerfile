FROM szabogtamas/jupy_rocker

RUN install2.r --error \
    --deps TRUE \
    devtools \
    readxl

RUN chmod a+rwx -R /home/rstudio

ADD ./configs/rstudio-prefs.json /home/rstudio/.config/rstudio/rstudio-prefs.json