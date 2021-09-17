# --------------------------------------------------------------
library(shinyjs)
library(DT)
# --------------------------------------------------------------

ui <- fluidPage(
  useShinyjs(),
  # Chroma.js
  singleton(tags$head(tags$script(src='https://cdnjs.cloudflare.com/ajax/libs/chroma-js/2.1.2/chroma.min.js', type='text/javascript'))),
  # TnT Tree
  singleton(tags$head(tags$link(href='http://tntvis.github.io/tnt.tree/build/tnt.tree.css', rel='stylesheet', type='text/css'))),
  singleton(tags$head(tags$script(src='http://d3js.org/d3.v3.min.js', type='text/javascript'))),
  singleton(tags$head(tags$script(src='http://tntvis.github.io/tnt.tree/build/tnt.tree.min.js', type='text/javascript'))),
  # NVD3 (for Taxa chart)
  singleton(tags$head(tags$link(href='https://cdn.rawgit.com/novus/nvd3/v1.8.1/build/nv.d3.css', rel='stylesheet', type='text/css'))),
  singleton(tags$head(tags$script(src='https://cdn.rawgit.com/novus/nvd3/v1.8.1/build/nv.d3.min.js', type='text/javascript'))),
  # MSA Viewer
  singleton(tags$head(tags$link(href='https://s3-eu-west-1.amazonaws.com/biojs/msa/latest/msa.min.gz.css', rel='stylesheet', type='text/css'))),
  singleton(tags$head(tags$script(src='https://s3-eu-west-1.amazonaws.com/biojs/msa/latest/msa.min.gz.js', type='text/javascript'))),
  # Ours (for phylotree and distance scale)
  singleton(tags$head(tags$link(href='css/phylogram.css', rel='stylesheet', type='text/css'))),
  # singleton(tags$head(tags$script(src='js/phylogram.js', type='text/javascript'))),
  extendShinyjs(script = "www/js/phylogram.js",
    functions = c("setPhylotree", "setPhylotreeLayout", "resetTaxa", "setMSA")
  ),

  # base HTML
  h2("Funnotate"),
  img(src = "rik.png", style = "horizontal-align: left; vertical-align: center"),
  HTML("<b>&ldquo;Now that&rsquo;s what I call <i>annotation</i>!&rdquo;</b>"),
  p(""),

  # Home page
  conditionalPanel("output.page == 'home'",
    HTML("<p>Annotate protein or nucleotide sequences using <a href='https://legumeinfo.org' target='_blank'>LIS</a> legume resources.</p>"),
    p("Because this service involves several computationally intensive searches (see pipeline description below), results can take from several minutes to several hours, depending on the size of your query. Thanks for your patience."),
    p("Upload your protein or nucleotide FASTA sequence(s) (max. 100 kbp)"),
    radioButtons("seqSource", label = NULL, choices = c("From text", "From file")),
    conditionalPanel("input.seqSource == 'From text'", textAreaInput("seqText", label = "Paste FASTA sequence(s) here:", width = "800px", height = "400px")),
    conditionalPanel("input.seqSource == 'From file'", fileInput("seqFile", label = NULL)),
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
  conditionalPanel("output.page == 'job'",
    htmlOutput("job2"),
    textOutput("uploadFile"),
    p(""),
    conditionalPanel("output.jobStatus == 'success'",
      tableOutput("simpleTable")
    ),
    htmlOutput("estscan"),
    htmlOutput("blast"),
    htmlOutput("ahrd"),
    htmlOutput("interpro"),
    htmlOutput("hmm"),
    htmlOutput("summary"),
    p(""),
    conditionalPanel("output.jobStatus == 'success'",
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
      )
    ),
    textOutput("jobDuration"),
    conditionalPanel("output.jobStatus != 'running'",
      actionLink("jhome", "Start Over")
    )
  ),

  # Phylogram page
  conditionalPanel("output.page == 'phylogram'",
    h3("Phylogram"),
    htmlOutput("phylogramFamilyInfo"),
    p(""),
    textOutput("phylogramStatus"),
    conditionalPanel("output.hasPhylotree == 'true'",
      HTML("<a href=https://legumeinfo.org/search/phylotree/userinfo target=_blank>Gene Family Help</a>"),
      checkboxGroupInput("phylogramToggleDisplay", label = NULL, inline = TRUE,
        choices = c("Taxa and Legend", "MSA Visualization"), selected = "Taxa and Legend"),
      # Taxa chart
      conditionalPanel("output.displayTaxaAndLegend == 'true'",
        hr(),
        actionButton("resetTaxa", label = "Reset Taxa Selection"),
        HTML("<svg id='taxa' height='300px'></svg>"),
        HTML(paste("<p style='font-size:9px; text-align: right;'>",
          "<a href='https://nvd3.org/' target='_blank'>NVD3</a>",
          " &bull; <a href='https://github.com/gka/chroma.js/' target='_blank'>Chroma.js</a>",
          "</p>"
        ))
      ),
      # MSA view
      conditionalPanel("output.displayMSA == 'true'",
        hr(),
        tags$div(id = "msa"),
        HTML(paste("<p style='font-size:9px; text-align: right;'>",
          "<a href='https://github.com/wilzbach/msa/' target='_blank'>MSA Viewer</a>",
          "</p>"
        ))
      ),
      # Phylotree
      hr(),
      radioButtons("phylotreeLayout", label = NULL, choices = c("Vertical layout", "Radial layout"), inline = TRUE),
      htmlOutput("phylotreeHilited"),
      HTML("<svg id='phylotreeDistanceScale'></svg>"),
      tags$div(id = "phylotree"),
      HTML(paste("<p style='font-size:9px; text-align: right;'>",
        "<a href='https://tntvis.github.io/tnt.tree/' target='_blank'>TnT Tree</a>",
        " &bull; <a href='https://github.com/gka/chroma.js/' target='_blank'>Chroma.js</a>",
        "</p>"
      )),
      hr()
    )
  )
)

# --------------------------------------------------------------

