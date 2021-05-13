# --------------------------------------------------------------
source("backend.R")
# --------------------------------------------------------------

server <- function(input, output, session) {
  rv <- reactiveValues()

  # Handle state from URL
  # output$page is a hidden output used for conditionalPanel conditions
  output$page <- reactive({
    whichPage <- "home"

    urlFields <- parseQueryString(session$clientData$url_search)
    numUrlFields <- length(urlFields)
    if (numUrlFields == 1) {
      if (!is.null(urlFields$upload)) {
        uploadIndex <- as.integer(urlFields$upload)
        if (!is.na(uploadIndex)) {
          existingUpload <- readUpload(uploadIndex)
          if (!is.null(existingUpload)) {
            whichPage <- "upload"
            rv$upload <- existingUpload
            displayUpload()
          }
        }
      } else if (!is.null(urlFields$job)) {
        existingJob <- readJob(urlFields$job)
        if (!is.null(existingJob)) {
          whichPage <- "job"
          rv$job <- existingJob
          displayJob()
        }
      }
    }
    if (whichPage == "home") clearQueryString()
    whichPage
  })
  outputOptions(output, "page", suspendWhenHidden = FALSE)

  clearQueryString <- function() {
    # TODO: replace with baseURL from settings
    updateQueryString("http://dev.lis.ncgr.org:50003/en/shiny/Funnotate/")
  }

  # TODO: use updateQueryString(q, mode = "push") ?
  observeEvent(c(input$uhome, input$jhome), {
    updateTextAreaInput(session, "seqText", value = "")
    reset("seqFile") # fileInput clears text, but input$seqFile data still exist
    output$page <- renderText("home")
    clearQueryString()
  }, ignoreInit = TRUE)

  displayUpload <- function() {
    output$uploadedFile <- renderText(sprintf("Uploaded file: %s (%d bytes)", rv$upload$inputFileName, rv$upload$inputFileSize))
    output$sequenceType <- renderText(paste("Sequence type:", rv$upload$sequenceType))
    output$numSequences <- renderText(paste("Number of sequences:", rv$upload$numSequences))
    output$totalSequenceLength <- renderText(paste("Total sequence length:", rv$upload$totalSequenceLength))
    output$uploadMessages <- renderUI(HTML(paste(rv$upload$messages, collapse = "<br>")))
  }

  observeEvent(input$upload, {
    seqText <- trimws(input$seqText)
    seqType <- substr(input$seqType, 1, 1)
    if (input$seqSource == "From text") {
      # read sequence(s) from text area
      seqSize <- nchar(seqText)
      if (seqSize == 0) return() # TODO: post error message
      # save to temporary file
      tempSeqFile <- "temp/pasted-text.fasta"
      write(seqText, tempSeqFile)
      seq <- list(name = "Pasted text", size = seqSize, datapath = tempSeqFile)
      rv$upload <- createNewUpload(seq, seqType)
      system(paste("rm", tempSeqFile))
    } else {
      # read sequence(s) from file
      if (is.null(input$seqFile)) return() # TODO: post error message
      rv$upload <- createNewUpload(input$seqFile, seqType)
    }
    displayUpload()

    output$page <- renderText("upload")
    updateQueryString(sprintf("?upload=%d", rv$upload$index))
  })

  displayJob <- function() {
if (isActive(rv$job)) {
    funnotize <- function(filename) gsub("www/", "", filename)

    output$job2 <- renderUI(h4(paste("Annotation Job", rv$job$id)))
    seqUnits <- ifelse(rv$job$sequenceType == "nucleotide", "base pairs", "amino acids")
    output$uploadFile <- renderText(sprintf("Uploaded file: %s (%s, %d sequences, %d total %s)",
      rv$job$inputFileName, rv$job$sequenceType, rv$job$numSequences, rv$job$totalSequenceLength, seqUnits))
    # output$simpleTable populated below
    if (rv$job$sequenceType == "nucleotide") {
      output$estscan <- renderUI(HTML(sprintf("<a href='%s.trans' target='_blank'>ESTScan output</a>", rv$job$uploadFile)))
    } else {
      output$estscan <- renderUI("")
    }
    output$blast <- renderUI(HTML(paste("BLAST output files:",
      "<ul>",
        sprintf("<li><a href='%s' target='_blank'>TAIR10_pep_20101209</a>", funnotize(rv$job$blastFiles[1])),
        sprintf("<li><a href='%s' target='_blank'>Glyma.refseq_protein.fasta</a>", funnotize(rv$job$blastFiles[2])),
        sprintf("<li><a href='%s' target='_blank'>Mt4.0v1_GenesProteinSeq_20130731_1800.fasta</a>", funnotize(rv$job$blastFiles[3])),
      "</ul>"
    )))
    if (rv$job$useInterpro) {
      output$interpro <- renderUI(HTML(sprintf("<a href='%s' target='_blank'>InterPro output</a> (.txt file, RAW format)", funnotize(rv$job$iprFile))))
    } else {
      output$interpro <- renderUI("")
    }
    output$hmm <- renderUI(HTML(sprintf("<a href='%s' target='_blank'>HMM output</a> (.tbl file, tabular)", funnotize(rv$job$hmmFile))))
    output$ahrd <- renderUI(HTML(sprintf("<a href='%s' target='_blank'>AHRD output</a> (.txt file, tabular)", funnotize(rv$job$ahrdFile))))

    # Summary table
    summary <- createSummaryTable(rv$job) # list of simpleTable, columnNames, summaryTable, summaryTableOut
    output$simpleTable <- renderTable(summary$simpleTable, colnames = FALSE)
    output$summaryLabel <- renderUI(HTML(sprintf("<b>Summary Table (<a href='%s' target='_blank'>download</a>)", funnotize(rv$job$summaryFile))))
    output$summaryTable <- renderDT(summary$summaryTable, rownames = FALSE,
      colnames = summary$columnNames, escape = FALSE, options = list(pageLength = 10))
    outputOptions(output, "summaryTable", suspendWhenHidden = FALSE)
    if (!file.exists(rv$job$summaryFile)) {
      write.table(summary$summaryTableOut, rv$job$summaryFile, sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = gsub("<sup>\\d+</sup>", "", summary$columnNames))
    }
}

    # Format job duration
    elapsedTime <- as.POSIXct(rv$job$endTime) - as.POSIXct(rv$job$startTime)
    ss <- as.integer(round(as.double(elapsedTime, units = "secs")))
    hh <- ss %/% 3600
    ss <- ss - hh*3600
    mm <- ss %/% 60
    ss <- ss - mm*60
    if (hh > 0) {
      jd <- sprintf("%d:%02d:%02d", hh, mm, ss)
    } else {
      jd <- sprintf("%d:%02d", mm, ss)
    }
    output$jobDuration <- renderText(sprintf("Job duration: %s", jd))
    output$jobMessages <- renderUI(HTML(paste(rv$job$messages, collapse = "<br>")))
  }

  observeEvent(input$job, {
    output$jobMessages <- renderUI(HTML(""))

    # Prepare and run the job
    rv$job <- createNewJob(rv$upload)
    rv$job$useInterpro <- input$useInterpro
    rv$job <- runJob(rv$job)
    displayJob()

    output$page <- renderText("job")
    updateQueryString(sprintf("?job=%s", rv$job$id))
  })
}

# --------------------------------------------------------------

