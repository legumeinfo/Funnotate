<!-- --------------------------------------------------------------------------------- -->

# Funnotate

R/Shiny application for functional annotation of user-supplied FASTA sequences (nucleotide or protein).

Uses third-party tools including

* [ESTScan](http://estscan.sourceforge.net) to translate nucleotide sequences
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) alignment to reference gene databases (soybean, _Medicago_, _Arabidopsis_)
* [InterProScan](https://www.ebi.ac.uk/interpro/search/sequence/) (methods: TIGRFAM, ProDom, SMART, SUPERFAMILY, PANTHER, Gene3D, PIRSF, Pfam, Coils)
* [HMMer](http://hmmer.org) for assigning gene family
* [AHRD](https://github.com/groupschoof/AHRD/blob/master/README.textile) - Automated Human Readable Descriptions

View the application [here](http://dev.lis.ncgr.org:50003/en/shiny/Funnotate).

To review an old upload, append ?upload=<_upload_index_>

To review an old job, append ?job=<_job_id_>

<!-- --------------------------------------------------------------------------------- -->

