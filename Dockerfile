FROM nfcore/base@sha256:2043dea2e3215a32576e2e9fa957d8d41f439d209abbf1f858fd02829d2b7d64
#FROM nfcore/base:1.10.2

LABEL authors="Christina Chatzipantsiou" \
      description="Docker image containing all software requirements for the siteqc pipeline"


# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/relate-1.0dev/bin:$PATH

# Install stringi R package and the ones that depend on it.
# (Issue with stringi package from conda that it depends on libicu64 that
# is not available for Debian 10)

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/stringi/stringi_1.4.5.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/stringr/stringr_1.3.1.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/plyr/plyr_1.8.5.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/reshape2/reshape2_1.4.3.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/tidyr/tidyr_1.0.2.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/BH/BH_1.72.0-2.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_0.8.5.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.12.8.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/magrittr_1.5.tar.gz', repos=NULL, type='source')"
RUN R -e "withr::with_makevars(c(PKG_LIBS = '-liconv'), install.packages('tidyverse',repos='http://cran.us.r-project.org'), assignment = '+=')"
RUN R -e "install.packages('randomForest',repos='http://cran.us.r-project.org')"

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-relate-1.0dev > nf-core-relate-1.0dev.yml
# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron


# Install GAWK
RUN apt-get update && \
    apt-get install -y \
                   gawk \
                   tabix \
