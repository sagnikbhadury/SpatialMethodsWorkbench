FROM rocker/shiny:4.5.2

RUN apt-get update && apt-get install -y --no-install-recommends \
    git libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev \
    libjpeg-dev libglpk-dev libgsl-dev cmake && \
    rm -rf /var/lib/apt/lists/*

RUN R -q -e "install.packages(c('shiny','bslib','DT','ggplot2','igraph','jsonlite','markdown','zip','nnet','testthat','remotes'), repos='https://cloud.r-project.org')"

COPY . /srv/shiny-server/spatial-methods-workbench

# The bundled engine sources are the exact versions the Workbench adapters target:
# ISPAT3D 0.3.1 (CRAN; provides both the planar ispat3d_fit_2d and volumetric ispat3d_fit paths)
# and GPGHS 1.0.1 (GitHub sagnikbhadury/GP-GHS). gpboost is ISPAT3D's compiled dependency.
ARG INSTALL_ADVANCED_ENGINES=false
RUN if [ "$INSTALL_ADVANCED_ENGINES" = "true" ]; then \
      R -q -e "install.packages('gpboost', repos='https://cloud.r-project.org'); remotes::install_local('/srv/shiny-server/spatial-methods-workbench/vendor/engines/GP-GHS', dependencies=FALSE, upgrade='never', build=FALSE); remotes::install_local('/srv/shiny-server/spatial-methods-workbench/vendor/engines/ISPAT-3D', dependencies=FALSE, upgrade='never', build=FALSE); stopifnot(packageVersion('ISPAT3D') >= '0.3.0', packageVersion('GPGHS') >= '1.0.0')"; \
    fi

RUN chown -R shiny:shiny /srv/shiny-server/spatial-methods-workbench

EXPOSE 3838
CMD ["/usr/bin/shiny-server"]
