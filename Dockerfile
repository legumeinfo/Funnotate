# Base R Shiny image
FROM rocker/shiny

# Bioconductor
FROM bioconductor/bioconductor_docker:devel

# Install R packages
RUN R -e "install.packages(c('shiny', 'BiocManager', 'cicerone', 'DT', 'future', 'httr', 'jsonlite', 'promises', 'shinyjs', 'stringi', 'stringr', 'yaml'))"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('Biostrings', 'InterMineR'))"

# Make required directories in the container
RUN mkdir /home/Funnotate
RUN mkdir /home/Funnotate/static
RUN mkdir /home/Funnotate/static/css
RUN mkdir /home/Funnotate/static/js

# Copy the Shiny application code
COPY app.R backend.R server.R settings.yml tour.R ui.R /home/Funnotate/

# Copy static files
COPY static/funnotate-process.png static/gene-family-help.html static/lis-6044923.png static/tools-512.png /home/Funnotate/static/
COPY static/css/phylogram.css /home/Funnotate/static/css/
COPY static/js/biojs-io-newick.min.js static/js/im.js static/js/phylogram.js /home/Funnotate/static/js/

# Expose the application's port
EXPOSE 8181

# Set user to shiny (instead of root)
#USER shiny

# Run the Shiny application
#WORKDIR /home/Funnotate
CMD Rscript /home/Funnotate/app.R

