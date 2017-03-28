FROM r-base
MAINTAINER ome-devel@lists.openmicroscopy.org.uk

ADD r-base-docker.R /
RUN Rscript /r-base-docker.R
RUN useradd -m user -d /scratch
WORKDIR /scratch
ENTRYPOINT ["/usr/bin/Rscript"]
