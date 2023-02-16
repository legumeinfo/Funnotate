# --------------------------------------------------------------
library(shinyjs)
library(Biostrings)
library(cicerone)
library(DT)
library(jsonlite)
library(stringi)
# --------------------------------------------------------------

addResourcePath("static", "./static")

proxyclick <- paste(
  "$(function() {",
    "var $els = $('[data-proxy-click]');",
    "$.each($els, function(idx, el) {",
      "var $el = $(el);",
      "var $proxy = $('#' + $el.data('proxyClick'));",
      "$el.keydown(function(e) {",
        "if (e.keyCode == 13) $proxy.click();",
      "});",
    "});",
  "})"
)

ui <- function(req) {
  # The `req` object is a Rook environment
  # See https://github.com/jeffreyhorner/Rook#the-environment
  if (identical(req$REQUEST_METHOD, "GET")) {
    fluidPage(
      useShinyjs(),
      use_cicerone(),
      # Font Awesome (for spinning progress icons)
      singleton(tags$head(tags$link(href='https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css', rel='stylesheet', type='text/css'))),
      # Chroma.js
      singleton(tags$head(tags$script(src='https://cdnjs.cloudflare.com/ajax/libs/chroma-js/2.1.2/chroma.min.js', type='text/javascript'))),
      # TnT Tree
      singleton(tags$head(tags$link(href='https://tntvis.github.io/tnt.tree/build/tnt.tree.css', rel='stylesheet', type='text/css'))),
      singleton(tags$head(tags$script(src='https://d3js.org/d3.v3.min.js', type='text/javascript'))),
      singleton(tags$head(tags$script(src='https://tntvis.github.io/tnt.tree/build/tnt.tree.min.js', type='text/javascript'))),
      # NVD3 (for Taxa chart)
      singleton(tags$head(tags$link(href='https://cdn.jsdelivr.net/gh/novus/nvd3@v1.8.1/build/nv.d3.min.css', rel='stylesheet', type='text/css'))),
      singleton(tags$head(tags$script(src='https://cdn.jsdelivr.net/gh/novus/nvd3@v1.8.1/build/nv.d3.min.js', type='text/javascript'))),
      # MSA Viewer
      singleton(tags$head(tags$link(href='https://s3-eu-west-1.amazonaws.com/biojs/msa/latest/msa.min.gz.css', rel='stylesheet', type='text/css'))),
      singleton(tags$head(tags$script(src='https://s3-eu-west-1.amazonaws.com/biojs/msa/latest/msa.min.gz.js', type='text/javascript'))),
      # jQuery UI (for dialogs)
      singleton(tags$head(tags$link(href='https://code.jquery.com/ui/1.12.1/themes/base/jquery-ui.css', rel='stylesheet', type='text/css'))),
      singleton(tags$head(tags$script(src='https://code.jquery.com/ui/1.12.1/jquery-ui.min.js', type='text/javascript'))),
      # biojs-io-newick (for exporting tree in Newick format)
      singleton(tags$head(tags$script(src='static/js/biojs-io-newick.min.js', type='text/javascript'))),
      # Ours (for phylotree and distance scale)
      singleton(tags$head(tags$link(href='static/css/phylogram.css', rel='stylesheet', type='text/css'))),
      singleton(tags$head(tags$script(src='static/js/phylogram.js', type='text/javascript'))),
      extendShinyjs(script = "static/js/phylogram.js",
        functions = c("setPhylotree", "setPhylotreeLayout", "showSingletonNodes", "clearSubtreeFocus", "resetTaxa", "setMSA")
      ),
      tags$head(tags$script(HTML(proxyclick))),

      # base HTML
      singleton(tags$head(tags$title("Funnotate"))),
      h2(HTML("<a href='https://legumeinfo.org' target='_blank'><img src='static/lis-6044923.png' width='64px' height='64px'></a> Funnotate")),
      tags$div(id = "loading", HTML("<p>Loading, please wait. <i class='fa fa-spinner fa-spin' style='font-size: 32px;'></i></p>")),

      # Home page
      conditionalPanel("output.page == 'home'",
        HTML("<p>Annotate protein or nucleotide sequences using <a href='https://legumeinfo.org' target='_blank'>LIS</a> legume resources, and identify homologous gene families.</p>"),
        HTML("<p>(Or search for gene families using <a href='?search'>this functional keyword search</a>.)</p>"),
        p("Because this service involves several computationally intensive searches (see pipeline description below), results can take from several minutes to several hours, depending on the size of your query. Thanks for your patience."),
        p("Upload your protein or nucleotide FASTA sequence(s) (max. 100 kbp)"),
        radioButtons("seqSource", label = NULL, choices = c("From text", "From file")),
        conditionalPanel("input.seqSource == 'From text'", textAreaInput("seqText", label = "Paste FASTA sequence(s) here:", width = "800px", height = "400px")),
        conditionalPanel("input.seqSource == 'From file'", fileInput("seqFile", label = NULL)),
        radioButtons("seqType", label = "Type of sequence", choices = c("nucleotide", "protein")),
        actionButton("upload", "Upload Sequence(s)"),
        p(""),
        wellPanel(HTML(paste(
          "Annotation pipeline:",
          "<br>",
          "<br>1. <a href='http://estscan.sourceforge.net' target='_blank'>ESTScan</a> if needed (to translate nucleotide uploads)",
          "<br>2. <a href='https://blast.ncbi.nlm.nih.gov/Blast.cgi' target='_blank'>BLAST</a> alignment to reference gene databases (soy, <i>Medicago</i>, <i>Arabidopsis</i>)",
          "<br>3. (optional) <a href='https://www.ebi.ac.uk/interpro/search/sequence/' target='_blank'>InterProScan</a> (methods: TIGRFAM, SMART, SUPERFAMILY, Gene3D, PIRSF, Pfam, Coils)",
          "<br>4. Assignment of <a href='https://legacy.legumeinfo.org/search/phylotree' target='_blank'>gene family</a> (using <a href='http://hmmer.org' target='_blank'>HMMer</a>)",
          "<br>5. Best-hit extraction by <a href='https://github.com/groupschoof/AHRD/blob/master/README.textile' target='_blank'>AHRD</a> (Automated Human Readable Descriptions)",
          "<br>",
          "<br>Results include: Phytozome gene family, AHRD descriptor, best BLAST hits, GO and InterPro ID (the latter two only if InterProScan is run)",
          "<br>Complete output files for each analysis are also provided.",
          "<br>",
          "<br>Note that the gene families which are searched are the families displayed on LIS, <i>i.e.</i> those having at least one member which is an LIS species."
        )))
      ),

      # Gene family functional keyword search page
      conditionalPanel("output.page == 'geneFamilySearch'",
        tagAppendAttributes(
          textInput("familyKeywords", label = "Gene Family Search:", width = "256px", placeholder = "e.g. cysteine"),
	  `data-proxy-click` = "familySearch"
	),
        HTML(paste(
          "<p>",
          "Enter keyword(s) to search for functional descriptions in gene families.",
          "<br>Use quotes to search for phrases (\"cysteine methyl\"),",
          "<br>* for partial matches (cyst*),",
          "<br>OR (AND) to search for either (both) of two terms (cysteine OR methyl; cysteine AND methyl),",
          "<br>NOT to exclude a term (cysteine AND NOT methyl).",
          "</p>"
        )),
        HTML("<p>(Or go to the <a href='.'>Funnotate home page</a> to annotate your sequence(s) and identify homologous gene families.)</p>"),
        actionButton("familySearch", label = "Search"),
        p(),
        textOutput("familySearchMessage"),
        DTOutput("familyTable")
      ),

      # Gene family selection page
      conditionalPanel("output.page == 'geneFamilySelection'",
        h4("Gene Family Selection"),
        p("Click a gene family to display its phylogram."),
        DTOutput("familySelectTable")
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
        h3("Phylogram", id = "tour-intro"),
        htmlOutput("phylogramFamilyInfo"),
        p(""),
        htmlOutput("phylogramStatus"),
        conditionalPanel("output.hasPhylotree == 'true'",
          actionLink("tour", "Quick Interactive Tour"),
          HTML("<span style='font-size: 9px'>(uses <a href='https://cicerone.john-coene.com' target='_blank'>Cicerone</a>)</span>"),
          HTML("&bull; <a href=https://legacy.legumeinfo.org/search/phylotree/userinfo target=_blank>Gene Family Help</a>"),
          checkboxGroupInput("phylogramToggleDisplay", label = NULL, inline = TRUE,
            choices = c("Taxa and Legend", "MSA Visualization", "Phylotree info"),
            selected = c("Taxa and Legend", "Phylotree info")),
          # Taxa chart
          conditionalPanel("output.displayTaxaAndLegend == 'true'",
            hr(),
            tags$div(id = "tour-taxaAndLegend",
              actionButton("resetTaxa", label = "Reset Taxa Selection"),
              HTML("<svg id='taxa' height='300px'></svg>"),
              HTML(paste("<p style='font-size:9px; text-align: right;'>",
                "<a href='https://nvd3.org/' target='_blank'>NVD3</a>",
                " &bull; <a href='https://github.com/gka/chroma.js/' target='_blank'>Chroma.js</a>",
                "</p>"
              ))
            )
          ),
          # MSA view
          conditionalPanel("output.displayMSA == 'true'",
            hr(),
            tags$div(id = "tour-msa",
              tags$div(id = "msa", style = "resize: vertical; overflow: auto;"),
              HTML(paste("<p style='font-size:9px; text-align: right;'>",
                "<a href='https://github.com/wilzbach/msa/' target='_blank'>MSA Viewer</a>",
                "</p>"
              ))
            )
          ),
          # Phylotree
          hr(),
          conditionalPanel("output.displayPhylotreeInfo == 'true'",
            radioButtons("phylotreeLayout", label = NULL, choices = c("Vertical layout", "Radial layout"), inline = TRUE),
            conditionalPanel("output.focusOnSubtree == 'true'",
              HTML("You have focused on a subtree."),
              actionButton("resetSubtreeFocus", label = "Reset to full tree"),
              div(style = "display: inline-block; vertical-align: middle;", checkboxInput("showSingletonNodes", label = "Show singleton nodes", value = TRUE))
            ),
          ),
          tags$div(id = "tour-phylotree",
            conditionalPanel("output.displayPhylotreeInfo == 'true'", htmlOutput("phylotreeHilited")),
            HTML("<svg id='phylotreeDistanceScale'></svg>"),
            tags$div(id = "phylotree"),
            HTML(paste("<p style='font-size:9px; text-align: right;'>",
              "<a href='https://tntvis.github.io/tnt.tree/' target='_blank'>TnT Tree</a>",
              " &bull; <a href='https://github.com/daviddao/biojs-io-newick'>biojs-io-newick</a>",
              " &bull; <a href='https://github.com/gka/chroma.js/' target='_blank'>Chroma.js</a>",
              " &bull; <a href='https://jqueryui.com' target='_blank'>jQuery UI</a>",
              "</p>"
            ))
          ),
          hr()
        )
      )
    )

  } else if (identical(req$REQUEST_METHOD, "POST")) {
    # Handle the POST
    postBytes <- req$rook.input$read(-1)
    postString <- rawToChar(postBytes)

    postData <- list(seqSource = "POSTed sequence(s)")
    postData$rawFasta <- URLdecode(stri_match_first(postString, regex="fasta=([^&]+)&?")[, 2])
    postData$type <- stri_match_first(postString, regex="type=([^&]+)&?")[, 2]
    geneFamily <- stri_match_first(postString, regex="geneFamily=([^&]+)&?")[, 2]
    postData$geneFamily <- stri_match_first(geneFamily, regex = "(L_[A-Z0-9]{6})")[, 2]
    seqFile <- tempfile()
    write(postData$rawFasta, seqFile)
    if (postData$type == "n") {
      fasta <- readDNAStringSet(seqFile)
    } else {
      fasta <- readAAStringSet(seqFile)
    }
    postData$seqNames <- names(fasta)
    # trim the type and gene_family from the sequence names
    ii <- unlist(gregexpr("[ |\\+]type=", postData$seqNames))
    postData$seqNames <- mapply(function(seqName, i) {
      if (i < 0) return(seqName)
      substr(seqName, 1, i - 1)
    }, postData$seqNames, ii, USE.NAMES = FALSE)
    postData$seqs <- as.character(fasta)
    unlink(seqFile)
    postJson <- toJSON(postData, auto_unbox = TRUE)

    fluidPage(
      useShinyjs(),
      # base HTML
      h2("Funnotate"),
      # TODO: add new logo here
      p(paste(length(postData$seqs), "sequence(s),", sum(nchar(postData$seqs)), "total characters")),
      p(paste("Sequence Type:", ifelse(postData$type == "n", "nucleotide", "protein"))),
      p(paste("Gene Family:", postData$geneFamily)),
      p(""),
      actionButton("phylogramFromPost", "Compute Phylogram", onclick = sprintf("Shiny.onInputChange('postData', %s);", postJson))
    )
  }
}

# Enable POST requests and return the ui function
attr(ui, "http_methods_supported") <- c("GET", "POST")
ui

# --------------------------------------------------------------

