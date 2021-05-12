# --------------------------------------------------------------
library(shinyjs)
library(DT)
# --------------------------------------------------------------

ui <- fluidPage(
  useShinyjs(),

  # base HTML
  h2("Funnotate"),
  img(src = "rik.png", style = "horizontal-align: left; vertical-align: center"),
  HTML("<b>&ldquo;Now that&rsquo;s what I call <i>annotation</i>!&rdquo;</b>"),
  p(""),

  # Home page
  conditionalPanel("output.page == 'home'",
    HTML("<p>Upload a file of protein or nucleotide sequences and have them annotated using <a href='https://legumeinfo.org' target='_blank'>LIS</a> legume resources.</p>"),
    p("Because this service involves several computationally intensive searches (see pipeline description below), results can take from several minutes to several hours, depending on the size of your query. Thanks for your patience."),
    p("Upload your protein or nucleotide FASTA file (max. 100 kbp)"),
    fileInput("browse", label = NULL),
    radioButtons("seqType", label = "Type of sequence", choices = c("nucleotide", "protein")),
    actionButton("upload", "Upload File"),
    p(""),
    wellPanel(HTML(paste(
      "Annotation pipeline:",
      "<br>",
      "<br>1. <a href='http://estscan.sourceforge.net' target='_blank'>ESTScan</a> if needed (to translate nucleotide uploads)",
      "<br>2. <a href='https://blast.ncbi.nlm.nih.gov/Blast.cgi' target='_blank'>BLAST</a> alignment to reference gene databases (soy, <i>Medicago</i>, <i>Arabidopsis</i>)",
      "<br>3. (optional) <a href='https://www.ebi.ac.uk/interpro/search/sequence/' target='_blank'>InterProScan</a> (methods: TIGRFAM, ProDom, SMART, SUPERFAMILY, PANTHER, Gene3D, PIRSF, Pfam, Coils)",
      "<br>4. Assignment of <a href='https://www.legumeinfo.org/search/phylotree' target='_blank'>gene family</a> (using <a href='http://hmmer.org' target='_blank'>HMMer</a>)",
      "<br>5. Best-hit extraction by <a href='https://github.com/groupschoof/AHRD/blob/master/README.textile' target='_blank'>AHRD</a> (Automated Human Readable Descriptions)",
      "<br>",
      "<br>Results include: Phytozome gene family, AHRD descriptor, best BLAST hits, GO and InterPro ID (the latter two only if InterProScan is run)",
      "<br>Complete output files for each analysis are also provided.",
      "<br>",
      "<br>Note that the gene families which are searched are the families displayed on LIS, <i>i.e.</i> those having at least one member which is an LIS species."
    )))
  ),

  # Upload page
  # TODO: show warnings and error messages
  conditionalPanel("output.page == 'upload'",
    textOutput("uploadedFile"),
    textOutput("sequenceType"),
    textOutput("numSequences"),
    textOutput("totalSequenceLength"),
    p(""),
    htmlOutput("uploadMessages"),
    p(tags$b("Please verify that this is correct before starting your job", style = "color: red"),
      "- or", actionLink("uhome", "start over")),
    checkboxInput("useInterpro", label = "Run InterPro analysis (significantly longer run time, but provides GO and other annotations)", value = FALSE, width = "100%"),
    actionButton("job", label = "Begin Annotation")
  ),

  # Job page
  # TODO: show intermediate results and error messages
  conditionalPanel("output.page == 'job'",
    htmlOutput("job2"),
    textOutput("uploadFile"),
    p(""),
    tableOutput("simpleTable"),
    htmlOutput("estscan"),
    htmlOutput("blast"),
    htmlOutput("interpro"),
    htmlOutput("hmm"),
    htmlOutput("ahrd"),
    p(""),
    htmlOutput("summaryLabel"),
    DTOutput("summaryTable"),
    HTML(paste("<br><sup>1</sup>GF Score = full-sequence E-value from hmmscan",
      "<br><sup>2</sup>BLAST Score: E-value is displayed, however best hit is chosen using bit score",
      "<br><sup>3</sup>AHRD Quality Scores are as follows:"
    )),
    wellPanel(
      HTML(paste("<p>",
        "The AHRD quality code consists of a three-character string,",
        "where each character is either * if the respective criteria are met",
        "or _ otherwise. The meaning by position is as follows:",
        "<br>1 - Bit score of the BLAST result is over 50 and e-value is better than e-10",
        "<br>2 - Overlap of the BLAST result is at least 60%",
        "<br>3 - Top token score of assigned descriptor is at least 0.5",
      "</p>")),
      p(""),
      p(paste("Note that the Best BLAST Hit may differ from the AHRD BLAST Hit",
        "because AHRD factors in the information content of the sequence description.")),
      p(""),
      HTML(paste("For further explanation of these codes and the AHRD algorithm, see",
        "<a href='https://github.com/groupschoof/AHRD/blob/master/README.textile' target='_blank'>AHRD Documentation</a>."))
    ),
    textOutput("jobDuration"),
    htmlOutput("jobMessages"),
    actionLink("jhome", "Start Over"),
    hr()
  )
)

# --------------------------------------------------------------

