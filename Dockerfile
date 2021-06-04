FROM ubuntu:21.04 AS ahrd

RUN apt update && apt install -y --no-install-recommends \
  ant \
  default-jdk \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src/
ARG AHRD_VERSION=3.3.3
ADD https://github.com/groupschoof/AHRD/archive/v${AHRD_VERSION}.tar.gz .
RUN tar -xzf v${AHRD_VERSION}.tar.gz \
  && cd AHRD-${AHRD_VERSION} \
  && ant dist \
  && mv dist/ahrd.jar /usr/src \
  && cd .. \
  && rm -rf AHRD-${AHRD_VERSION}

FROM ubuntu:21.04

RUN apt update \
  && DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends \
  estscan \
  default-jre-headless \
  hmmer \
  ncbi-blast+ \
  r-bioc-biostrings \
  r-cran-dt \
  r-cran-future \
  r-cran-promises \
  r-cran-shiny \
  r-cran-shinyjs \
  r-cran-stringi \
  r-cran-yaml \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /srv/shiny-server
COPY --from=ahrd /usr/src/ahrd.jar .

RUN Rscript -e 'download.file(url="https://legumeinfo.org/data/v2/LEGUMES/Fabaceae/genefamilies/legume.genefam.fam1.M65K/legume.genefam.fam1.M65K.hmm.tar.gz", destfile="legume.genefam.fam1.M65K.hmm.tar.gz")' \
  && tar -xzf legume.genefam.fam1.M65K.hmm.tar.gz \
  && rm legume.genefam.fam1.M65K.hmm.tar.gz

RUN mkdir blastdb
RUN Rscript -e 'download.file(url="https://www.arabidopsis.org/download_files/Proteins/TAIR10_protein_lists/TAIR10_pep_20101214", destfile="blastdb/TAIR10_pep_20101214")'
RUN cd blastdb && makeblastdb -parse_seqids -dbtype prot -taxid 3702 -in TAIR10_pep_20101214

RUN Rscript -e 'download.file(url="https://de.cyverse.org/anon-files/iplant/home/mtruncatula/public/Mt4.0/Annotation/Mt4.0v2/Mt4.0v2_GenesProteinSeq_20140818_1100.fasta", destfile="blastdb/Mt4.0v2_GenesProteinSeq_20140818_1100.fasta")'
RUN cd blastdb && makeblastdb -parse_seqids -dbtype prot -taxid 3880 -in Mt4.0v2_GenesProteinSeq_20140818_1100.fasta

RUN Rscript -e 'download.file(url="https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/3847/103/GCF_000004515.5_Glycine_max_v2.1/GCF_000004515.5_Glycine_max_v2.1_protein.faa.gz", destfile="blastdb/GCF_000004515.5_Glycine_max_v2.1_protein.faa.gz")' \
  && gzip -d blastdb/GCF_000004515.5_Glycine_max_v2.1_protein.faa.gz \
  && makeblastdb -parse_seqids -dbtype prot -taxid 3847 -in blastdb/GCF_000004515.5_Glycine_max_v2.1_protein.faa

COPY . .

ENTRYPOINT ["Rscript", "-e", "library(shiny); runApp(host='0.0.0.0', port=3838)"]

EXPOSE 3838
