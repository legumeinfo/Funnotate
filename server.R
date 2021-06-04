# --------------------------------------------------------------
library(promises)
library(future)
plan(multisession) # or multicore?

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
            displayUpload(existingUpload)
          }
        }
      } else if (!is.null(urlFields$job)) {
        existingJob <- readJob(urlFields$job)
        if (!is.null(existingJob)) {
          whichPage <- "job"
          displayJob(existingJob)
        }
      }
    }
    if (whichPage == "home") clearQueryString()
    whichPage
  })
  outputOptions(output, "page", suspendWhenHidden = FALSE)

  # output$jobStatus is a hidden output used for conditionalPanel conditions
  output$jobStatus <- reactive("new")
  outputOptions(output, "jobStatus", suspendWhenHidden = FALSE)

  clearQueryString <- function() {
    updateQueryString("")
  }

  # TODO: use updateQueryString(q, mode = "push") ?
  observeEvent(c(input$uhome, input$jhome), {
    updateTextAreaInput(session, "seqText", value = "")
    reset("seqFile") # fileInput clears text, but input$seqFile data still exist
    output$page <- renderText("home")
    clearQueryString()
  }, ignoreInit = TRUE)

  displayUpload <- function(upload) {
    output$uploadedFile <- renderText(sprintf("Uploaded file: %s (%d bytes)", upload$inputFileName, upload$inputFileSize))
    output$sequenceType <- renderText(paste("Sequence type:", upload$sequenceType))
    output$numSequences <- renderText(paste("Number of sequences:", upload$numSequences))
    output$totalSequenceLength <- renderText(paste("Total sequence length:", upload$totalSequenceLength))
    output$uploadMessages <- renderUI(HTML(paste(upload$messages, collapse = "<br>")))
  }

  observeEvent(input$upload, {
    seqText <- trimws(input$seqText)
    seqType <- substr(input$seqType, 1, 1)
    if (input$seqSource == "From text") {
      # read sequence(s) from text area
      seqSize <- nchar(seqText)
      if (seqSize == 0) return() # TODO: post error message
      # save to temporary file
      tempSeqFile <- tempfile(pattern="pasted-text", fileext="fasta")
      write(seqText, tempSeqFile)
      seq <- list(name = "Pasted text", size = seqSize, datapath = tempSeqFile)
      upload <- createNewUpload(seq, seqType)
      unlink(tempSeqFile)
    } else {
      # read sequence(s) from file
      if (is.null(input$seqFile)) return() # TODO: post error message
      upload <- createNewUpload(input$seqFile, seqType)
    }
    rv$upload <- upload
    displayUpload(upload)

    output$page <- renderText("upload")
    updateQueryString(sprintf("?upload=%d", rv$upload$index))
  })

  clearJobPage <- function() {
    output$job2 <- renderUI("")
    output$uploadFile <- renderText("")
    output$estscan <- renderUI("")
    output$blast <- renderUI("")
    output$interpro <- renderUI("")
    output$hmm <- renderUI("")
    output$ahrd <- renderUI("")
    output$jobDuration <- renderText("")
    output$jobMessages <- renderUI("")
  }

  getJobMessages <- function(job) {
    if (is.null(job)) {
      msgs <- c("The annotation job is queued and will begin shortly.",
        "The page will automatically refresh.")
      return(paste(msgs, collapse = "<br>"))
    } else if (!isDone(job)) {
      msgs <- c("Job status:", job$estStatus, job$blastStatus,
        job$ahrdStatus, job$iprStatus, job$hmmStatus, job$summaryStatus)
      msgs <- msgs[nchar(msgs) > 0]
      return(paste(msgs, collapse = "<br>"))
    }
    paste(job$messages, collapse = "<br>")
  }

  displayJob <- function(job) {
    output$jobStatus <- renderText(ifelse(isRunning(job), "running", "done"))
    clearJobPage()

    output$job2 <- renderUI(h4(paste("Annotation Job", job$id)))
    seqUnits <- ifelse(job$sequenceType == "nucleotide", "base pairs", "amino acids")
    output$uploadFile <- renderText(sprintf("Uploaded file: %s (%s, %d sequences, %d total %s)",
      job$inputFileName, job$sequenceType, job$numSequences, job$totalSequenceLength, seqUnits))

    # output$simpleTable populated below
    if (job$sequenceType == "nucleotide") {
      output$estscan <- renderUI(HTML(
        sprintf("<a href='%s' target='_blank'>ESTScan output</a>", funnotize(job$inputFile))
      ))
    }
    output$blast <- renderUI(HTML(
      paste("BLAST output files:",
        "<ul>",
        paste(sprintf("<li><a href='%s' target='_blank'>%s</a>",
          funnotize(job$blastFiles), basename(settings$blast$dbs)), collapse = " "),
        "</ul>"
      )
    ))
    if (job$useInterpro) {
      output$interpro <- renderUI(HTML(sprintf("<a href='%s' target='_blank'>InterPro output</a> (.txt file, RAW format)", funnotize(job$iprFile))))
    }
    output$hmm <- renderUI(HTML(sprintf("<a href='%s' target='_blank'>HMM output</a> (.tbl file, tabular)", funnotize(job$hmmFile))))
    output$ahrd <- renderUI(HTML(sprintf("<a href='%s' target='_blank'>AHRD output</a> (.txt file, tabular)", funnotize(job$ahrdFile))))

    # Summary table
    if (job$status == "success") {
      summary <- createSummaryTable(job) # list of simpleTable, columnNames, summaryTable, summaryTableOut
      output$simpleTable <- renderTable(summary$simpleTable, colnames = FALSE)
      output$summaryLabel <- renderUI(HTML(sprintf("<b>Summary Table (<a href='%s' target='_blank'>download</a>)", funnotize(job$summaryFile))))
      output$summaryTable <- renderDT(summary$summaryTable, rownames = FALSE,
        colnames = summary$columnNames, escape = FALSE, options = list(pageLength = 10))
      if (!file.exists(job$summaryFile)) {
        write.table(summary$summaryTableOut, job$summaryFile, sep = "\t", quote = FALSE,
          row.names = FALSE, col.names = gsub("<sup>\\d+</sup>", "", summary$columnNames))
      }
    }

    # Job status (if running) or messages (if success/failure)
    output$jobMessages <- renderUI({
      jobNow <- readJob(job$id)
      if (isRunning(jobNow)) invalidateLater(4000) # put delay time in settings?
      HTML(getJobMessages(jobNow))
    })

    # Format job duration
    output$jobDuration <- renderText({
      if (isDone(job)) {
        latest <- job$endTime
      } else {
        invalidateLater(4000) # put delay time in settings?
        latest <- Sys.time()
      }
      sprintf("Job duration: %s", elapsedTimeString(job$startTime, latest))
    })
  }

  observeEvent(input$job, {
    # Prepare and run the job
    job <- createNewJob(rv$upload, isolate(input$useInterpro))
    jobId <- job$id
    # For now (promises 1.0.1), use future() but eventually (promises 1.2.0) upgrade to future_promise()
    result <- future({
      runJob(job)
    }) %...>% {
      job <- readJob(jobId)
      displayJob(job)
    }

    # job <- readJob(jobId)
    displayJob(job)
    output$page <- renderText("job")
    updateQueryString(sprintf("?job=%s", job$id))

    # Return something other than the promise, so that Shiny remains responsive.
    # NULL
  })
}

# --------------------------------------------------------------

