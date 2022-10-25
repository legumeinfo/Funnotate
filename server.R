# --------------------------------------------------------------
library(promises)
library(future)
plan(multisession) # or multicore?

source("backend.R")
source("tour.R")
# --------------------------------------------------------------

server <- function(input, output, session) {
  rv <- reactiveValues()

  # Handle state from URL
  # output$page is a hidden output used for conditionalPanel conditions
  # TODO: refactor so as to test only once for existingJob
  output$page <- reactive({
    whichPage <- "home"

    urlFields <- parseQueryString(session$clientData$url_search)
    family <- urlFields$family
    if (!is.null(family)) {
      family <- URLdecode(family)
      # strip full-yuck prefix (if any) to extract the L_xxxxxx part of the family name
      family <- stri_match_first(urlFields$family, regex=".*(L_.+$)")[, 2]
    }
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
      } else if (!is.null(family)) {
        nullJob <- NULL
        existingPhylogram <- buildUserPhylogram(nullJob, family)
        if (!is.null(existingPhylogram)) {
          whichPage <- "phylogram"
          displayPhylogram(nullJob, existingPhylogram)
        }
      } else if (!is.null(urlFields$gene_name)) {
        df_gene_families <- genesToGeneFamiliesQuery(urlFields$gene_name)
        numRows <- nrow(df_gene_families)
        if (df_gene_families$description[numRows] == "Not found" || numRows > 1) {
          # post gene family selection page
          whichPage <- "geneFamilySelection"
          output$familySelectTable <- renderDT({
            dt <- datatable(df_gene_families[, c("link", "description", "genes")], escape = FALSE,
              options = list(pageLength = 10), rownames = FALSE, colnames = c("Gene Family", "Description", "Genes")
            )
            dt
          })
        } else {
          # go directly to phylogram page for this family, as for (family & gene_name) query below
          family <- df_gene_families$geneFamily[1]
          # strip full-yuck prefix (if any) to extract the L_xxxxxx part of the family name
          family <- stri_match_first(family, regex=".*(L_.+$)")[, 2]
          nullJob <- NULL
          existingPhylogram <- buildUserPhylogram(nullJob, family)
          if (!is.null(existingPhylogram)) {
            whichPage <- "phylogram"
            gene_name <- df_gene_families$gene_name[1]
            proteins <- genesToProteinsQuery(family, gene_name)
            existingPhylogram$seqNames <- proteins
            displayPhylogram(nullJob, existingPhylogram, proteins)
            # reset URL
            updateQueryString(sprintf("?family=%s&gene_name=%s", family, gene_name))
          }
	}
      } else if (!is.null(urlFields$search)) {
        # gene family functional keyword search
        whichPage <- "geneFamilySearch"
        keywords <- trimws(urlFields$search)
        updateTextInput(session, "familyKeywords", value = keywords)
      }
    } else if (numUrlFields == 2) {
      if (!is.null(urlFields$job) && !is.null(family)) {
        # Create Lorax tree for sequences that match the selected gene family
        existingJob <- readJob(urlFields$job)
        if (!is.null(existingJob)) {
          whichPage <- "phylogram"
          loraxResults <- buildUserPhylogram(existingJob, family)
          displayPhylogram(existingJob, loraxResults)
        }
      } else if (!is.null(family) && !is.null(urlFields$gene_name)) {
        nullJob <- NULL
        existingPhylogram <- buildUserPhylogram(nullJob, family)
        if (!is.null(existingPhylogram)) {
          whichPage <- "phylogram"
          proteins <- genesToProteinsQuery(family, urlFields$gene_name)
          existingPhylogram$seqNames <- proteins
          displayPhylogram(nullJob, existingPhylogram, proteins)
        }
      }
    }

    # Clean up and return page
    runjs("$('#loading').html('');") # clear page loading message
    if (whichPage == "home") clearQueryString()
    whichPage
  })
  outputOptions(output, "page", suspendWhenHidden = FALSE)

  # output$jobStatus is a hidden output used for conditionalPanel conditions
  output$jobStatus <- reactive("new")
  outputOptions(output, "jobStatus", suspendWhenHidden = FALSE)

  clearQueryString <- function() {
    colon_port <- ifelse(is.null(session$clientData$url_port), "", paste0(":", session$clientData$url_port))
    updateQueryString(sprintf("%s//%s%s%s", session$clientData$url_protocol, session$clientData$url_hostname, colon_port, session$clientData$url_pathname))
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
      seqFile <- tempfile()
      write(seqText, seqFile)
      seq <- list(name = "Pasted sequence(s)", size = seqSize, datapath = seqFile)
      upload <- createNewUpload(seq, seqType)
      unlink(seqFile)
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
    output$ahrd <- renderUI("")
    output$interpro <- renderUI("")
    output$hmm <- renderUI("")
    output$summary <- renderUI("")
    output$jobDuration <- renderText("")
  }

  statusDone <- function(taskStatus) {
    grepl("Done", taskStatus)
  }
  statusFailed <- function(taskStatus) {
    grepl("Failed", taskStatus)
  }
  # Color to display failed task, depends on whether job is salvageable
  statusColor <- function(jobStatus) {
    ifelse(jobStatus == "failure", "red", "orange")
  }

  displayJob <- function(job) {
    if (job$status == "success") {
      output$jobStatus <- renderText("success")
    } else if (isRunning(job)) {
      output$jobStatus <- renderText("running")
    } else {
      output$jobStatus <- renderText("failed")
    }
    clearJobPage()

    if (!is.null(job$geneFamily)) {
      # go straight to phylogram without computing gene family
      output$page <- reactive({
        loraxResults <- buildUserPhylogram(job, job$geneFamily)
        displayPhylogram(job, loraxResults)
        "phylogram"
      })
      updateQueryString(sprintf("?job=%s&family=%s", job$id, job$geneFamily))
      return
    }

    output$job2 <- renderUI(h4(paste("Annotation Job", job$id)))
    seqUnits <- ifelse(job$sequenceType == "nucleotide", "base pairs", "amino acids")
    output$uploadFile <- renderText(sprintf("Uploaded file: %s (%s, %d sequences, %d total %s)",
      job$inputFileName, job$sequenceType, job$numSequences, job$totalSequenceLength, seqUnits))

    # output$simpleTable populated below
    if (job$sequenceType == "nucleotide") {
      output$estscan <- renderUI({
        jobNow <- readJob(job$id)
        if (isRunning(jobNow)) invalidateLater(4000) # put delay time in settings?
        if (statusDone(jobNow$estscanStatus)) {
          estscanResult <- sprintf("<a href='%s' target='_blank'>ESTScan output</a>", jobNow$inputFile)
        } else if (statusFailed(jobNow$estscanStatus)) {
          estscanResult <- sprintf("<span style='color: %s'>%s</span>", statusColor(jobNow$status), jobNow$estscanStatus)
        } else {
          estscanResult <- jobNow$estscanStatus
        }
        HTML(estscanResult)
      })
    }

    output$blast <- renderUI({
      jobNow <- readJob(job$id)
      if (isRunning(jobNow)) invalidateLater(4000) # put delay time in settings?
      blastResults <- sapply(1:length(settings$blast$dbs), function(i) {
        if (statusDone(jobNow$blastStatus[i])) {
          result <- sprintf("<li><a href='%s' target='_blank'>%s</a>", jobNow$blastFiles[i], basename(settings$blast$dbs[i]))
        } else if (statusFailed(jobNow$blastStatus[i])) {
          result <- sprintf("<li><span style='color: %s'>%s</span>", statusColor(jobNow$status), jobNow$blastStatus[i])
        } else {
          result <- sprintf("<li>%s", jobNow$blastStatus[i])
        }
        result
      })
      HTML(paste("BLAST output",
        "<ul>",
        paste(blastResults, collapse = " "),
        "</ul>"
      ))
    })

    output$ahrd <- renderUI({
      jobNow <- readJob(job$id)
      if (isRunning(jobNow)) invalidateLater(4000) # put delay time in settings?
      if (statusDone(jobNow$ahrdStatus)) {
        ahrdResult <- sprintf("<a href='%s' target='_blank'>AHRD output</a> (.txt file, tabular)", jobNow$ahrdFile)
      } else if (statusFailed(jobNow$ahrdStatus)) {
        ahrdResult <- sprintf("<span style='color: %s'>%s</span>", statusColor(jobNow$status), jobNow$ahrdStatus)
      } else {
        ahrdResult <- jobNow$ahrdStatus
      }
      HTML(ahrdResult)
    })

    if (job$useInterpro) {
      output$interpro <- renderUI({
        jobNow <- readJob(job$id)
        if (isRunning(jobNow)) invalidateLater(4000) # put delay time in settings?
        if (statusDone(jobNow$iprStatus)) {
          iprResult <- ifelse(file.info(jobNow$iprFile)$size == 0, "No InterPro output",
            sprintf("<a href='%s' target='_blank'>InterPro output</a> (.txt file, RAW format)", jobNow$iprFile))
        } else if (statusFailed(jobNow$iprStatus)) {
          iprResult <- sprintf("<span style='color: %s'>%s</span>", statusColor(jobNow$status), jobNow$iprStatus)
        } else {
          iprResult <- jobNow$iprStatus
        }
        HTML(iprResult)
      })
    }

    output$hmm <- renderUI({
      jobNow <- readJob(job$id)
      if (isRunning(jobNow)) invalidateLater(4000) # put delay time in settings?
      if (statusDone(jobNow$hmmStatus)) {
        hmmResult <- sprintf("<a href='%s' target='_blank'>HMM output</a> (.tbl file, tabular)", jobNow$hmmFile)
      } else if (statusFailed(jobNow$hmmStatus)) {
        hmmResult <- sprintf("<span style='color: %s'>%s</span>", statusColor(jobNow$status), jobNow$hmmStatus)
      } else {
        hmmResult <- jobNow$hmmStatus
      }
      HTML(hmmResult)
    })

    output$summary <- renderUI({
      jobNow <- readJob(job$id)
      if (jobNow$status == "success") {
        # Summary table
        summary <- createSummaryTable(jobNow) # list of simpleTable, columnNames, summaryTable, summaryTableOut
        output$simpleTable <- renderTable(summary$simpleTable, colnames = FALSE)
        output$summaryLabel <- renderUI(HTML(sprintf("<b>Summary Table (<a href='%s' target='_blank'>download</a>)", jobNow$summaryFile)))
        output$summaryTable <- renderDT(summary$summaryTable, rownames = FALSE,
          colnames = summary$columnNames, escape = FALSE, options = list(pageLength = 10))
        if (!file.exists(jobNow$summaryFile)) {
          write.table(summary$summaryTableOut, jobNow$summaryFile, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = gsub("<sup>\\d+</sup>", "", summary$columnNames))
        }
      } else if (!statusDone(jobNow$summaryStatus)) {
        invalidateLater(4000) # put delay time in settings?
        HTML(jobNow$summaryStatus)
      }
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
    job <- createNewJob(rv$upload, input$useInterpro)
    jobId <- job$id
    # For now (promises 1.0.1), use future() but eventually (promises 1.2.0) upgrade to future_promise()
    future({
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

  observeEvent(input$phylogramFromPost, {
    # Directly create the phylogram from the input sequences and gene family

    # The logic behind putting the (nonstandard) upload here is to prevent InterMine users
    # from creating lots of upload files. Alternatively, we could use a cron job to clean them up.
    seq <- list(name = input$postData$seqSource, size = nchar(input$postData$rawFasta), # datapath = NULL,
      seqNames = unlist(input$postData$seqNames), sequences = unlist(input$postData$seqs))
    rv$upload <- createNewUpload(seq, input$postData$type)

    job <- createNewJobWithGeneFamily(rv$upload, input$postData$geneFamily)
    jobId <- job$id
    runJobWithGeneFamily(job)
    job <- readJob(jobId)

    # Create a minimal summary table (just Query and Gene.Family columns)
    df.summary <- data.frame(Query = sort(seq$seqNames), Gene.Family = input$postData$geneFamily, stringsAsFactors = FALSE)
    if (!file.exists(job$summaryFile)) {
      write.table(df.summary, job$summaryFile, sep = "\t", quote = FALSE, row.names = FALSE)
    }

    output$page <- reactive({
      loraxResults <- buildUserPhylogram(job, input$postData$geneFamily)
      displayPhylogram(job, loraxResults)
      "phylogram"
    })
    updateQueryString(sprintf("?job=%s&family=%s", job$id, input$postData$geneFamily))
    runjs("window.location.reload();")
  })

  displayPhylogram <- function(job, phylogramInfo, proteins = NULL) {
    if (is.null(phylogramInfo)) return()
    # Otherwise, display phylotree for selected gene family, with user sequence(s) inserted
    if (!phylogramInfo$done) {
      invalidateLater(3000) # put delay time in settings?
    } else {
      if (!is.null(phylogramInfo$tree)) output$hasPhylotree <- renderText("true")
    }

    # Display status/error message, if any
    output$phylogramStatus <- renderText(phylogramInfo$message)

    # Add links to gene family/IPR/GO descriptors
    gfds <- phylogramInfo$descriptor
    ii <- na.omit(unlist(stri_match_all(gfds, regex = "IPR\\d+")))
    for (ipr in ii) {
      link.ipr <- sprintf("<a href='http://www.ebi.ac.uk/interpro/entry/%s' target='_blank'>%s</a>", ipr, ipr)
      gfds <- gsub(ipr, link.ipr, gfds)
    }
    gg <- na.omit(unlist(stri_match_all(gfds, regex = "GO:\\d+")))
    for (go in gg) {
      link.go <- sprintf("<a href='http://amigo.geneontology.org/amigo/term/%s' target='_blank'>%s</a>", go, go)
      gfds <- gsub(go, link.go, gfds)
    }
    output$phylogramFamilyInfo <- renderUI(HTML(
      sprintf("<b>%s%s</b>: %s", legfed_prefix, phylogramInfo$family, paste(gfds, collapse = "; "))
    ))

    if (!is.null(phylogramInfo$tree)) {
      if (!is.null(phylogramInfo$seqNames)) {
        output$phylotreeHilited <- renderUI(HTML(paste(
          "Jump to highlighted feature:",
          paste(sprintf("<br>&bull; <a onclick='onScrollToHilite(\"%s\");'>%s</a>", phylogramInfo$seqNames, phylogramInfo$seqNames), collapse = "")
        )))
      }
      rv$tree <- phylogramInfo$tree
      js$setPhylotree(phylogramInfo$tree, "phylotree", paste(proteins, collapse = " ")) # sets tree data for both phylotree and taxa chart
    }

    if (!is.null(phylogramInfo$msa)) {
      js$setMSA(phylogramInfo$msa, "msa")
    }
  }

  observeEvent(input$resetTaxa, {
    js$resetTaxa("taxa")
  })

  observeEvent(input$phylotreeLayout, {
    layout <- stri_match_first(input$phylotreeLayout, regex = "^(.+) layout$")[, 2]
    js$setPhylotreeLayout(layout, "phylotree");
  })

  output$hasPhylotree <- reactive("false")
  outputOptions(output, "hasPhylotree", suspendWhenHidden = FALSE)

  output$focusOnSubtree <- reactive("false")
  outputOptions(output, "focusOnSubtree", suspendWhenHidden = FALSE)
  observeEvent(input$toggleSubtree, {
    output$focusOnSubtree <- renderText(ifelse(input$toggleSubtree, "true", "false"))
  })
  observeEvent(input$showSingletonNodes, {
    js$showSingletonNodes(input$showSingletonNodes, "phylotree")
  })
  observeEvent(input$resetSubtreeFocus, {
    js$clearSubtreeFocus("phylotree")
  })

  output$displayTaxaAndLegend <- reactive("false")
  output$displayMSA <- reactive("false")
  outputOptions(output, "displayTaxaAndLegend", suspendWhenHidden = FALSE)
  outputOptions(output, "displayMSA", suspendWhenHidden = FALSE)

  observeEvent(input$phylogramToggleDisplay, {
    output$displayTaxaAndLegend <- renderText(
      ifelse("Taxa and Legend" %in% input$phylogramToggleDisplay, "true", "false")
    )
    output$displayMSA <- renderText(
      ifelse("MSA Visualization" %in% input$phylogramToggleDisplay, "true", "false")
    )
  })

  observeEvent(input$familyKeywords, {
    # automatically click Search button if user entered keywords in the URL (which resets button value to 0)
    if (input$familySearch == 0 && nchar(input$familyKeywords) > 0) click("familySearch")
  })
  observeEvent(input$familySearch, {
    familyResults <- geneFamilySearchQuery(input$familyKeywords)
    nResults <- ifelse(is.null(familyResults), 0, nrow(familyResults))
    plural <- ifelse(nResults == 1, "", "s")
    what <- ifelse(nResults == 1, "its", "any")
    clickOne <- ifelse(nResults == 0, "", sprintf(" Click %s identifier to visualize its phylogram and MSA.", what))
    msg <- sprintf("%d result%s found for '%s'.%s", nResults, plural, input$familyKeywords, clickOne)
    output$familySearchMessage <- renderText(msg)
    output$familyTable <- renderDT({
      dt <- datatable(familyResults, escape = FALSE, options = list(pageLength = 100),
        rownames = FALSE, colnames = c("Identifier", "Description", "Score")
      )
      if (!is.null(familyResults)) dt <- dt %>% formatRound(3, 2) # otherwise an error occurs
      dt
    })
    # reset URL
    updateQueryString("?search")
  })

  observeEvent(input$tour, {
    updateCheckboxGroupInput(session, "phylogramToggleDisplay", selected = c("Taxa and Legend", "MSA Visualization"))
    tour$init()$start()
  })
}

# --------------------------------------------------------------

