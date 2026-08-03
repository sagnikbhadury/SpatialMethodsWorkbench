FROM rocker/shiny:4.5.2

RUN apt-get update && apt-get install -y --no-install-recommends \
    git libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev \
    libjpeg-dev libglpk-dev libgsl-dev && \
    rm -rf /var/lib/apt/lists/*

RUN R -q -e "install.packages(c('shiny','bslib','DT','ggplot2','jsonlite','markdown','zip','testthat','remotes'), repos='https://cloud.r-project.org')"

COPY . /srv/shiny-server/spatial-methods-workbench

ARG INSTALL_ADVANCED_ENGINES=false
RUN if [ "$INSTALL_ADVANCED_ENGINES" = "true" ]; then \
      R -q -e "remotes::install_local('/srv/shiny-server/spatial-methods-workbench/vendor/engines/GP-GHS', dependencies=TRUE, upgrade='never'); remotes::install_local('/srv/shiny-server/spatial-methods-workbench/vendor/engines/ISPAT-3D', dependencies=TRUE, upgrade='never'); remotes::install_local('/srv/shiny-server/spatial-methods-workbench/vendor/engines/ISPAT', dependencies=TRUE, upgrade='never')"; \
    fi

RUN chown -R shiny:shiny /srv/shiny-server/spatial-methods-workbench

EXPOSE 3838
CMD ["/usr/bin/shiny-server"]
