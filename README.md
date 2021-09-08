<!-- --------------------------------------------------------------------------------- -->

# Funnotate

R/Shiny application for functional annotation of user-supplied FASTA sequences (nucleotide or protein).

Uses third-party tools including

* [ESTScan](http://estscan.sourceforge.net) to translate nucleotide sequences
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) alignment to reference gene databases (soybean, _Medicago_, _Arabidopsis_)
* [InterProScan](https://www.ebi.ac.uk/interpro/search/sequence/) (methods: TIGRFAM, SMART, SUPERFAMILY, Gene3D, PIRSF, Pfam, Coils)
* [HMMer](http://hmmer.org) for assigning gene family
* [AHRD](https://github.com/groupschoof/AHRD/blob/master/README.textile) - Automated Human Readable Descriptions

<img src="www/funnotate-process.png" width="400px" height="272px">

Clicking a <img src="www/tools-512.png" width="16px" height="16px"> icon in the Gene Family column of the summary table calls [Lorax](https://github.com/LegumeFederation/lorax) to compute that family&rsquo;s phylogenetic tree. The resulting Phylogram page visualizations use these JavaScript libraries:

* [Chroma.js](https://github.com/gka/chroma.js/)
* [NVD3](https://nvd3.org/)
* [MSAViewer](https://github.com/wilzbach/msa/)
* [TnT Tree](https://tntvis.github.io/tnt.tree/index.html)

View the application [here](http://dev.lis.ncgr.org:50003/en/shiny/Funnotate/).

To review an old upload, append ?upload=<_upload_index_>

To review an old job, append ?job=<_job_id_>

<!-- --------------------------------------------------------------------------------- -->

