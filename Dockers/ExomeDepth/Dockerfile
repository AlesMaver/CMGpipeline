# Base Image
FROM bioconductor/bioconductor_docker

# R dependencies
RUN Rscript -e "install.packages('BiocManager')" && \
	Rscript -e "BiocManager::install(c('Biostrings', 'IRanges', 'Rsamtools', 'GenomicRanges', 'GenomicAlignments', 'rtracklayer'), update=F)" && \
	Rscript -e "install.packages('aod')" && \
	Rscript -e "install.packages('VGAM')" && \
	Rscript -e "install.packages('dplyr')" && \
	Rscript -e "install.packages('tidyverse')" && \
	Rscript -e "install.packages('optparse')" && \
	Rscript -e "install.packages('doParallel')" && \
	Rscript -e "install.packages('RcppRoll')" && \
	Rscript -e "install.packages(pkgs = 'https://cran.r-project.org/package=ExomeDepth&version=1.1.16', repos = NULL)"

# Create a directory /home/exome_depth and copy ExomeDepth.R script into it
RUN mkdir /home/exome_depth
COPY ExomeDepth.R /home/exome_depth

CMD ["/bin/bash"]