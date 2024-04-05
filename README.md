<!-- --------------------------------------------------------------------------------- -->

# Funnotate

R/Shiny application for functional annotation of user-supplied FASTA sequences (nucleotide or protein).

Uses third-party tools including

* [ESTScan](http://estscan.sourceforge.net) to translate nucleotide sequences (if needed)
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) alignment to reference gene databases (soybean, _Medicago_, _Arabidopsis_)
* [InterProScan](https://www.ebi.ac.uk/interpro/search/sequence/) (methods: TIGRFAM, SMART, SUPERFAMILY, Gene3D, PIRSF, Pfam, Coils)
* [HMMer](http://hmmer.org) for assigning gene family
* [AHRD](https://github.com/groupschoof/AHRD/blob/master/README.textile) (Automated Human Readable Descriptions) for best-hit extraction

<img src="static/funnotate-process.png" width="400px" height="272px">

Clicking a <img src="static/tools-512.png" width="16px" height="16px"> icon in the Gene Family column of the summary table calls [Lorax](https://github.com/LegumeFederation/lorax) to compute that family&rsquo;s phylogenetic tree. The resulting Phylogram page visualizations use these JavaScript libraries:

* [biojs-io-newick](https://github.com/daviddao/biojs-io-newick)
* [Chroma.js](https://github.com/gka/chroma.js/)
* [jQuery UI](https://jqueryui.com)
* [NVD3](https://nvd3.org/)
* [MSAViewer](https://github.com/wilzbach/msa/)
* [TnT Tree](https://tntvis.github.io/tnt.tree/index.html)

View the application [here](https://funnotate.legumeinfo.org).

### URL query strings

To review a previously uploaded sequence set (for example, to launch a new annotation job), append
<br>&nbsp;&nbsp; _?upload=<_upload_index_>_
<br>to the base URL.

To review an existing annotation job, append
<br>&nbsp;&nbsp; _?job=<_job_id_>_

To view the phylogram for a gene family with sequences inserted by an existing job, append
<br>&nbsp;&nbsp; _?family=<_family_>&job=<_job_id_>_

To view the phylogram for a gene family alone (with no user-supplied sequences), append
<br>&nbsp;&nbsp; _?family=<_family_>_
<br>To additionally highlight the proteins associated with specified genes from that family, append
<br>&nbsp;&nbsp; _?family=<_family_>&gene_name=<_gene1_>,<_gene2_>,<...>,<_geneN_>_

To automatically determine the gene family associated with one or more specified genes, append
<br>&nbsp;&nbsp; _?gene_name=<_gene1_>,<_gene2_>,<...>,<_geneN_>_
<br>If this finds more than one unique gene family, it will display the results in a table. Choose one to display its phylogram, with any related proteins highlighted.

To search for gene families by functional keywords, append
<br>&nbsp;&nbsp; _?search_
<br>or
<br>&nbsp;&nbsp; _?search=<_keywords_>_

### Example FASTA sequences

To use example FASTA files as input selections on the Funnotate home page (for example, to avoid having to hunt for them during a presentation), place them in your `static/examples/` subdirectory.

<!-- --------------------------------------------------------------------------------- -->

