FROM ubuntu:16.04

RUN apt-get update && apt-get install -y \
   libxml2-dev \
   libcurl4-openssl-dev \
   pandoc\
   xorg \
   libx11-dev \
   libglu1-mesa-dev \
   libssl-dev \
   gcc \
   g++

RUN apt-get -y install software-properties-common
RUN add-apt-repository ppa:webupd8team/java
RUN apt-get update

# Install python libraries
RUN apt-get -y install python-pip
RUN pip install docopt==0.6.2
RUN pip install beautifulsoup4
RUN pip install lxml
RUN pip install path.py

# R install https://www.r-bloggers.com/how-to-install-r-on-linux-ubuntu-16-04-xenial-xerus/
RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 | apt-key add -
RUN apt-get update
RUN apt-get install -y \
   r-base-core=3.4.2-2xenial2 \
   r-base-dev=3.4.2-2xenial2\
   r-cran-rmysql

# set the repos and make sure that the Bioconductor release matches the chosen version of R
# https://bioconductor.org/about/release-announcements/
# MRAN snapshot date should be set to a date before the next R version was uploaded here:
# https://cran.r-project.org/bin/linux/ubuntu/trusty/
# ...or the present date if the latest R is used
RUN echo 'source("https://bioconductor.org/biocLite.R")' > /opt/packages.R
RUN echo 'options(repos = c(' >> /opt/packages.R
RUN echo 'CRAN = "http://mran.revolutionanalytics.com/snapshot/2017-11-09/",' >> /opt/packages.R
RUN echo 'BioCsoft = "http://bioconductor.statistik.tu-dortmund.de/packages/3.6/bioc",' >> /opt/packages.R
RUN echo 'BioCann = "http://bioconductor.statistik.tu-dortmund.de/packages/3.6/data/annotation",' >> /opt/packages.R
RUN echo 'BioCexp = "http://bioconductor.statistik.tu-dortmund.de/packages/3.6/data/experiment",' >> /opt/packages.R
RUN echo 'BioCextra = "http://bioconductor.statistik.tu-dortmund.de/packages/3.6/extra"' >> /opt/packages.R
RUN echo '))' >> /opt/packages.R

# installing required R packages from Bioconductor
RUN echo 'biocLite("tximport")' >> /opt/packages.R
RUN echo 'biocLite("GenomicFeatures")' >> /opt/packages.R

# installing required R packages from CRAN
#RUN echo 'options(repos = "http://mran.revolutionanalytics.com/snapshot/2017-11-11/")' > /opt/packages.R
RUN echo 'install.packages("sp")' >> /opt/packages.R
RUN echo 'install.packages("pixmap")' >> /opt/packages.R
RUN echo 'install.packages("snowfall")' >> /opt/packages.R
RUN echo 'install.packages("VGAM")' >> /opt/packages.R
RUN echo 'install.packages("mclust")' >> /opt/packages.R
RUN echo 'install.packages("logcondens")' >> /opt/packages.R
RUN echo 'install.packages("Iso")' >> /opt/packages.R
RUN echo 'install.packages("XML")' >> /opt/packages.R
RUN echo 'install.packages("rgl")' >> /opt/packages.R
RUN echo 'install.packages("plyr")' >> /opt/packages.R
RUN echo 'install.packages("git2r")' >> /opt/packages.R
RUN echo 'install.packages("httr")' >> /opt/packages.R
RUN echo 'install.packages("gplots")' >> /opt/packages.R
RUN echo 'install.packages("rjson")' >> /opt/packages.R
RUN echo 'install.packages("knitr")' >> /opt/packages.R
RUN echo 'install.packages("rmarkdown")' >> /opt/packages.R
RUN echo 'install.packages("ggplot2")' >> /opt/packages.R
RUN echo 'install.packages("devtools")' >> /opt/packages.R
RUN echo 'source("http://www.math.ntnu.no/inla/givemeINLA.R")' >> /opt/packages.R
RUN echo 'library(devtools)' >> /opt/packages.R
RUN echo 'install_github("markvdwiel/ShrinkBayes")' >> /opt/packages.R
RUN Rscript /opt/packages.R

COPY Dockerfile /opt/

MAINTAINER Viktorian Miok, Seven Bridges, <viktorian.miok@sbgenomics.com>
