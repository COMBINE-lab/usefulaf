# image: COMBINE-lab/dockeraf
FROM condaforge/mambaforge:4.10.3-7
MAINTAINER salmon.maintainer@gmail.com

RUN mamba install -c anaconda git
RUN mamba install -c bioconda salmon alevin-fry gffread
RUN mamba install -c conda-forge cxx-compiler
RUN mamba install -c conda-forge r-base r-essentials
RUN mamba install -c bioconda bioconductor-eisar bioconductor-biostrings bioconductor-bsgenome bioconductor-genomicfeatures bioconductor-singlecellexperiment
RUN mamba install -c conda-forge r-stringr r-argparser r-rjson r-matrix
RUN mamba install -c conda-forge time
RUN mamba install -c conda-forge r-devtools
RUN R -e "devtools::install_github('COMBINE-lab/roe')"

RUN touch /root/.bashrc \
       && echo "export ALEVIN_FRY_HOME=/workdir/.afhome" >> /root/.bashrc \
       && echo "export TIME_BIN=/opt/conda/bin/time" >> /root/.bashrc

RUN touch /root/.Rprofile 

RUN echo 'local({r <- getOption("repos") \n\
      r["CRAN"] <- "http://cran.r-project.org" \n\
      options(repos=r) \n\
      })' >> /root/.Rprofile

RUN git clone https://github.com/COMBINE-lab/usefulaf.git

ENV ALEVIN_FRY_HOME=/workdir/.afhome
ENV TIME_BIN=/opt/conda/bin/time
