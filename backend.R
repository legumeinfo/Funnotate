# --------------------------------------------------------------
library(Biostrings)
library(httr)
library(jsonlite)
library(stringi)
library(yaml)
# --------------------------------------------------------------

# Global settings
# setwd("/srv/shiny-server/Funnotate") # already there
settings <- read_yaml("settings.yml")

# Global data:
# Read table of InterPro data
df.interpro <- read.table("https://ftp.ebi.ac.uk/pub/databases/interpro/entry.list",
  header = TRUE, sep = "\t", colClasses = "character", stringsAsFactors = FALSE)
# Build table of GO terms
getGOTerms <- function() {
  ll.go <- readLines("https://ftp.ebi.ac.uk/pub/databases/interpro/interpro2go")
  ll.go <- ll.go[!startsWith(ll.go, "!")]
  mx.go <- t(sapply(ll.go, function(go) {
    stri_match_first(go, regex = "^.+GO:(.+) ; GO:(\\d+)$")
  }, USE.NAMES = FALSE))
  df.go <- as.data.frame(mx.go[, c(3, 2)])
  names(df.go) <- c("GO.term", "name")
  # deduplicate
  df.go <- unique(df.go)
  df.go
}
df.goterms <- getGOTerms()

# --------------------------------------------------------------

# For scrubbing upload files
scrubber <- list()
scrubber$n <- list(
  sequenceType = "nucleotide",
  goodChars = "acgtn",
  changeTo = "n"
)
scrubber$p <- list(
  sequenceType = "protein",
  goodChars = "abcdefghiklmnpqrstuvwxyz",
  changeTo = "x"
)

# Returns total number of bad characters in ss0
# compared to scrubbed characters in ss1
countBadChars <- function(ss0, ss1) {
  # assumes ss0 and ss1 are the same length,
  # as ss1 = ss0 with characters replaced
  cc0 <- unlist(strsplit(ss0, split = ""))
  cc1 <- unlist(strsplit(ss1, split = ""))
  sum(cc0 != cc1)
}

# isProbablyNucleotideSequence <- function(ss) {
#   ns <- nchar(ss)
# }

# File exists and is not empty
fileReallyExists <- function(filename) {
  fileExists <- file.exists(filename)
  result <- (fileExists && file.info(filename)$size > 0)
  if (fileExists && !result) system(paste("rm", filename))
  result
}
hmmFileReallyExists <- function(filename) {
  fileExists <- fileReallyExists(filename)
  if (!fileExists) return(FALSE)
  ll <- readLines(filename)
  ll <- ll[!startsWith(ll, "#")]
  result <- (length(ll) > 0)
  if (!result) system(paste("rm", filename))
  result
}

# Returns the file path relative to "www/"
funnotize <- function(filename) {
  gsub("www/", "", filename)
}

# --------------------------------------------------------------

# Current and elapsed time as string
currentTimeString <- function() {
  as.character(Sys.time())
}

elapsedTimeString <- function(t0, t1) {
  elapsedTime <- as.POSIXct(t1) - as.POSIXct(t0)
  ss <- as.integer(round(as.double(elapsedTime, units = "secs")))
  hh <- ss %/% 3600
  ss <- ss - hh*3600
  mm <- ss %/% 60
  ss <- ss - mm*60
  if (hh > 0) {
    ets <- sprintf("%d:%02d:%02d", hh, mm, ss)
  } else {
    ets <- sprintf("%d:%02d", mm, ss)
  }
  ets
}

# --------------------------------------------------------------

# fiData is the value of a Shiny fileInput control,
# seqType = n for nucleotide/DNA, p for protein sequence
createNewUpload <- function(fiData, seqType) {
  scrub <- scrubber[[seqType]]

  uploadDir <- "www/upload"
  if (!dir.exists(uploadDir)) {
    dir.create(uploadDir)
  }
  indexFile <- sprintf("%s/index", uploadDir)
  if (!file.exists(indexFile)) {
    write("1", indexFile)
  }
  nextIndex <- scan(indexFile, what = integer())
  write(as.character(nextIndex + 1), indexFile)

  # Create the upload - fields set to NA will be filled in below
  upload <- list(
    index = nextIndex,
    inputFileName = fiData$name,
    inputFileSize = fiData$size,
    seqType = seqType,
    sequenceType = scrub$sequenceType,
    numSequences = NA,
    totalSequenceLength = NA,
    totalBadChars = NA,
    uploadFile = sprintf("%s/upload_%d", uploadDir, nextIndex), # metadata for existing uploads
    inputFileScrubbed = sprintf("%s/input_%d", uploadDir, nextIndex),
    messages = c()
  )

  # Read the sequences (to a Biostrings::DNAStringSet or Biostrings::AAStringSet)
  if (seqType == "n") {
    fasta <- readDNAStringSet(fiData$datapath)
  } else {
    fasta <- readAAStringSet(fiData$datapath)
  }
  seqNames <- names(fasta)
  sequences <- as.character(fasta)
  upload$numSequences <- length(sequences) # or length(seqNames)
  upload$totalSequenceLength <- sum(width(sequences))
#  probablyNucleotides <- all(!is.na(stri_match_first(sequences, regex = "^[acgtn]*$")))
#  updateRadioButtons(session, "seqType", value = ifelse(probablyNucleotides, "nucleotide", "protein"))

  # Scrub the uploaded sequences
  badChars <- sprintf("(?i)[^%s]", scrub$goodChars)
  sequencesScrubbed <- gsub(badChars, scrub$changeTo, sequences)
  upload$totalBadChars <- countBadChars(sequences, sequencesScrubbed)
  if (seqType == "n") {
    fastaScrubbed <- DNAStringSet(sequencesScrubbed)
  } else {
    fastaScrubbed <- AAStringSet(sequencesScrubbed)
  }
  names(fastaScrubbed) <- seqNames
  writeXStringSet(fastaScrubbed, upload$inputFileScrubbed, width = 60) # TODO: match original width if possible

  # Error checking
  if (upload$totalSequenceLength > 100000) {
    upload$messages <- c(upload$messages, "Total sequence length exceeds 100 kbp.")
  }
  if (upload$totalBadChars > 0) {
    upload$messages <- c(upload$messages,
      sprintf("%d invalid characters (other than %s) were found. These will be changed to '%s'.",
        upload$totalBadChars, scrub$goodChars, scrub$changeTo)
    )
  }
  if (seqType == "n") {
    bb <- (width(sequences) %% 3 != 0)
    if (FALSE && any(bb)) {
      upload$messages <- c(upload$messages, "The following sequences do not look like nucleotide sequences:", sequences[bb])
    }
    # } else if (...) {
    #   scrubbingErrors <- append(scrubbingErrors, list(
    #     "The following sequences do not look like protein sequences:",
    #     sequences[bb]
    #   ))
  }

  write_yaml(upload, upload$uploadFile)
  upload
}

# --------------------------------------------------------------

createNewJob <- function(upload, useInterpro) {
  # Generate an as-yet-unused job id
  jobsDir <- "www/job"
  dd <- list.dirs(jobsDir)
  while (TRUE) {
    # Why is job id like "A1B2C"? Andrew says it is to guard against "obscenities".
    aa <- sample(LETTERS, 3, replace = TRUE)
    nn <- sample(1:9, 2, replace = TRUE) # why not allow 0s? Easily confused with Os?
    jobId <- paste(aa[1], nn[1], aa[2], nn[2], aa[3], sep = "")
    if (!(jobId %in% dd)) break
  }
  jobDir <- sprintf("%s/%s", jobsDir, jobId)

  # Create the job - fields set to NA will be filled in later
  job <- list(
    id = jobId,
    dir = jobDir,
    inputFileName = upload$inputFileName,
    inputFileSize = upload$inputFileSize,
    sequenceType = upload$sequenceType,
    numSequences = upload$numSequences,
    totalSequenceLength = upload$totalSequenceLength,
    totalBadChars = upload$totalBadChars,
    useInterpro = useInterpro,
    inputFile = upload$inputFileScrubbed,
    jobFile = sprintf("%s/job_%s", jobDir, jobId), # metadata for existing jobs
    blastStatus = sprintf("%s: Queued", basename(settings$blast$dbs)),
    ahrdStatus = "AHRD: Queued",
    hmmStatus = "HMMer: Queued",
    summaryStatus = "Postprocessing: Queued",
    blastFiles = sprintf("%s/blast_%s_%d", jobDir, jobId, 1:length(settings$blast$dbs)),
    ahrdFile = sprintf("%s/ahrd_%s.txt", jobDir, jobId),
    hmmFile = sprintf("%s/hmm_%s.tbl", jobDir, jobId),
    summaryFile = sprintf("%s/summary_%s.txt", jobDir, jobId),
    status = "new",
    startTime = currentTimeString(),
    endTime = ""
  )
  # Insert optional fields in the right position
  if (upload$sequenceType == "nucleotide") {
    j <- which(names(job) == "jobFile")
    job <- append(job, list(estscanStatus = "ESTScan: Queued"), j)
  }
  if (job$useInterpro) {
    j <- which(names(job) == "ahrdStatus")
    job <- append(job, list(iprStatus = "InterPro: Queued"), j)
    j <- which(names(job) == "ahrdFile")
    job <- append(job, list(iprFile = sprintf("%s/ipr_%s.txt", job$dir, job$id)), j)
  }
  job
}

isActive <- function(job) {
  job$status != "failure"
}
isDone <- function(job) {
  job$status %in% c("failure", "success")
}
isRunning <- function(job) {
  is.null(job) || !isDone(job)
}

failJob <- function(job) {
  job$status <- "failure"
  job$endTime <- currentTimeString()
  job
}

# --------------------------------------------------------------

runESTScan <- function(job) {
  inputFileTrans <- paste0(job$inputFile, ".trans")
  estscanCmd <- sprintf("%s -M %s -t %s %s", # -o /dev/null
    settings$estscan$exe, settings$estscan$matrix, inputFileTrans, job$inputFile)
  job$estscanStatus <- "ESTScan: Running"
  writeJob(job)
  system(estscanCmd)
  # Remove (first) semicolon after sequence name, if any
  system(sprintf("perl -pi -e 's/;//' %s", inputFileTrans))
  # TODO: Clean up the header (sequence name)?
  if (fileReallyExists(inputFileTrans)) {
    j <- which(names(job) == "inputFile") - 1
    job <- append(job, list(originalInputFile = job$inputFile), j)
    job$inputFile <- inputFileTrans
    job$estscanStatus <- "ESTScan: Done"
  } else {
    job$estscanStatus <- "ESTScan: Failed, no output file"
    job$status <- "failure"
    job$endTime <- currentTimeString()
  }
  writeJob(job)
  job
}

runBLAST <- function(job) {
  for (i in 1:length(settings$blast$dbs)) {
    blastDb.i <- basename(settings$blast$dbs[i])
    job$blastStatus[i] <- sprintf("%s: Running", blastDb.i)
    writeJob(job)
    #blastCmd.i <- sprintf("%s -db %s -query %s -out %s -outfmt 6 -num_threads %d",
    #  settings$blast$exe, settings$blast$dbs[i], inputFile, job$blastFiles[i], settings$num_threads)
    blastCmd.i <- sprintf("%s -p blastp -d %s -i %s -o %s -m 8",
    # blastCmd.i <- sprintf("%s -p blastp -d %s -i %s -o %s -e 0.0001 -v 200 -b 200 -m 0 -a 4",
      settings$blast$exe, settings$blast$dbs[i], job$inputFile, job$blastFiles[i])
    system(blastCmd.i)
    if (fileReallyExists(job$blastFiles[i])) {
      job$blastStatus[i] <- sprintf("%s: Done", blastDb.i)
    } else {
      job$blastStatus[i] <- sprintf("%s: Failed, no output file", blastDb.i)
      #job$status <- "failure"
      job$endTime <- currentTimeString()
    }
    writeJob(job)
  }
  if (sum(file.exists(job$blastFiles)) == 0) {
    job$status <- "failure"
    writeJob(job)
  }
  job
}

runAHRD <- function(job) {
  ahrdTmpYmlFile <- NA # for cleanup
  # Read YAML file and update certain parameters
  ahrdYml <- read_yaml(settings$ahrd$yml)
  ahrdYml$proteins_fasta <- job$inputFile
  ahrdYml$output <- job$ahrdFile
  job$ahrdStatus <- "AHRD: Running"
  writeJob(job)
  b <- 1
  for (bdb in ahrdYml$blast_dbs) {
    # AHRD expects the databases in FASTA (text) format
    db <- bdb$database
    # match db with index
    i <- which(basename(settings$blast$dbs) == basename(db))
    ahrdYml$blast_dbs[[b]]$file <- job$blastFiles[i]
    b <- b + 1
  }
  # remove missing blast_dbs
  bb <- sapply(ahrdYml$blast_dbs, function(bdb) file.exists(bdb$file))
  ahrdYml$blast_dbs <- ahrdYml$blast_dbs[bb]
  # write the modified YAML to a temporary file
  ahrdTmpYmlFile <- sprintf("temp/ahrd_%s.yml", job$id)
  write_yaml(ahrdYml, ahrdTmpYmlFile)
  ahrdCmd <- sprintf("%s -Xmx2g -jar %s %s", settings$ahrd$java, settings$ahrd$jar, ahrdTmpYmlFile)
  system(ahrdCmd)
  # clean up temporary file
  if (!is.na(ahrdTmpYmlFile) && file.exists(ahrdTmpYmlFile)) system(paste("rm", ahrdTmpYmlFile))
  # done
  if (fileReallyExists(job$ahrdFile)) {
    job$ahrdStatus <- "AHRD: Done"
  } else {
    job$ahrdStatus <- "AHRD: Failed, no output file"
    job$status <- "failure"
    job$endTime <- currentTimeString()
  }
  writeJob(job)
  job
}

runInterPro <- function(job) {
  iprXml <- sprintf("temp/ipr_%s.xml", job$id)
  iprCmdXml <- sprintf("%s -i %s -o %s -f XML %s", settings$interpro$exe, job$inputFile, iprXml, settings$interpro$params)
  job$iprStatus <- "InterPro: Running"
  writeJob(job)
  system(iprCmdXml)
  if (file.exists(iprXml)) {
    iprCmdRaw <- sprintf("%s -i %s -mode convert -f RAW -o %s", settings$interpro$exe, iprXml, job$iprFile)
    system(iprCmdRaw)
    # clean up temporary file
    if (!is.na(iprXml) && file.exists(iprXml)) system(paste("rm", iprXml))
    # done
    if (fileReallyExists(job$iprFile)) {
      job$iprStatus <- "InterPro: Done"
    } else {
      job$iprStatus <- "InterPro: Failed to convert XML to raw (txt) - no matches"
      #job$status <- "failure"
      job$endTime <- currentTimeString()
    }
  } else {
    job$iprStatus <- "InterPro: Failed, no XML output"
    #job$status <- "failure"
    job$endTime <- currentTimeString()
  }
  writeJob(job)
  job
}

runHMMer <- function(job) {
  hmmCmd <- sprintf("%s --cpu %d --tblout %s %s %s",
    settings$hmmer$exe, settings$num_threads, job$hmmFile, settings$hmmer$db, job$inputFile)
  job$hmmStatus <- "HMMer: Running"
  writeJob(job)
  system(hmmCmd)
  if (hmmFileReallyExists(job$hmmFile)) {
    job$hmmStatus <- "HMMer: Done"
  } else {
    job$hmmStatus <- "HMMer: Failed, no output"
    #job$status <- "failure"
    job$endTime <- currentTimeString()
  }
  writeJob(job)
  job
}

runJob <- function(job) {
  if (!dir.exists(job$dir)) dir.create(job$dir)

  # ESTScan
  if (job$sequenceType == "nucleotide") {
    if (isActive(job)) job <- runESTScan(job)
  }

  # BLAST
  if (isActive(job)) job <- runBLAST(job)

  # AHRD
  if (isActive(job)) job <- runAHRD(job)

  # InterPro
  if (job$useInterpro) {
    if (isActive(job)) job <- runInterPro(job)
  }

  # HMMer
  if (isActive(job)) job <- runHMMer(job)

  # Job completed!
  if (isActive(job)) {
    job$status <- "success"
    job$endTime <- currentTimeString()
  }
  writeJob(job)
  job
}

# --------------------------------------------------------------

# Read an existing upload
readUpload <- function(index) {
  uploadFile <- sprintf("www/upload/upload_%d", index)
  if (!file.exists(uploadFile)) return(NULL)
  upload <- read_yaml(uploadFile)
  upload
}

# Read an existing job
readJob <- function(jobId) {
  jobFile <- sprintf("www/job/%s/job_%s", jobId, jobId)
  if (!file.exists(jobFile)) return(NULL)
  job <- read_yaml(jobFile)
  job
}

# Write an existing job
writeJob <- function(job) {
  write_yaml(job, job$jobFile)
}

# --------------------------------------------------------------

# Read output files and create summary table
createSummaryTable <- function(job) {
  if (!endsWith(job$summaryStatus, "Done")) {
    job$summaryStatus <- "Postprocessing: Running"
    writeJob(job)
  }

  # AHRD
  df.ahrd <- read.table(job$ahrdFile, skip = 2, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df.summary <- df.summary.txt <- df.ahrd[, 1:4]
  colnames.summary <- c("Query", "AHRD BLAST Hit", "AHRD Quality<sup>3</sup>", "AHRD Descriptor")
  blankColumn <- rep("", nrow(df.summary))

  # InterPro
  if (job$useInterpro) {
    if (!file.exists(job$iprFile)) {
      df.i <- data.frame(iit1 = blankColumn, go1 = blankColumn, iit2 = blankColumn, go2 = blankColumn, stringsAsFactors = FALSE)
    } else {
      df.ipr <- read.table(job$iprFile, header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
      # Loop over all sequences, match with first column of df.ahrd
      df.i <- as.data.frame(do.call(rbind, lapply(df.ahrd[, 1], function(q) {
        df.ii <- df.ipr[df.ipr[, 1] == q, ]
        iit <- setdiff(unique(df.ii[, 12]), "NULL")
        go <- stri_match_all(df.ii[, 14], regex = "\\(GO:(\\d+)\\)")
        go <- setdiff(unique(unlist(sapply(go, function(g) g[, 2], USE.NAMES = FALSE))), NA)
        iit1 <- ifelse(length(iit) == 0, "", paste(sprintf("<a href='https://www.ebi.ac.uk/interpro/entry/%s' target='_blank'>%s</a>", iit, iit), collapse = ", "))
        go1 <- ifelse(length(go) == 0, "", paste(sprintf("<a href='http://amigo.geneontology.org/amigo/term/GO:%s' target='_blank'>%s</a>", go, go), collapse = ", "))
        iit2 <- ifelse(length(iit) == 0, "", paste(iit, collapse = ","))
        go2 <- ifelse(length(go) == 0, "", paste(go, collapse = ","))
        # InterPro id, GO terms (HTML and text format)
        c(iit1, go1, iit2, go2)
      })), stringsAsFactors = FALSE)
    }
    df.summary <- cbind(df.summary, data.frame(iit = df.i[, 1], go = df.i[, 2], stringsAsFactors = FALSE))
    df.summary.txt <- cbind(df.summary.txt, data.frame(iit = df.i[, 3], go = df.i[, 4], stringsAsFactors = FALSE))
    colnames.summary <- c(colnames.summary, "InterPro-ID", "GO Terms")
    nGo <- sum(df.i[, 2] != "")
  } else {
    nGo <- 0
  }
  df.simpleTable <- data.frame(
    label = c("Number of annotated sequences", "Number with GO assignment"),
    value = c(
      sprintf("%d (%2.1f%%)", nrow(df.summary), 100*nrow(df.summary)/job$numSequences),
      sprintf("%d (%2.1f%%)", nGo, 100*nGo/job$numSequences)
    ),
    stringsAsFactors = FALSE
  )

  # HMMer
  if (!file.exists(job$hmmFile)) {
    df.h <- data.frame(gf1 = blankColumn, gfs1 = blankColumn, gf2 = blankColumn, gfs2 = blankColumn, stringsAsFactors = FALSE)
  } else {
    ll.hmm <- readLines(job$hmmFile)
    ll.hmm <- ll.hmm[!startsWith(ll.hmm, "#")]
    df.hmm <- as.data.frame(do.call(rbind, strsplit(ll.hmm, split = "\\s+")), stringsAsFactors = FALSE)
    # Loop over all sequences, match with first column of df.ahrd
    df.h <- as.data.frame(do.call(rbind, lapply(df.ahrd[, 1], function(q) {
      df.hi <- df.hmm[df.hmm[, 3] == q, ]
      if (nrow(df.hi) == 0) {
        gf1 <- gfs1 <- gf2 <- gfs2 <- ""
      } else {
        family <- df.hi[1, 1]
        gf1 <- paste(
          sprintf("<a href='https://legumeinfo.org/chado_phylotree/legfed_v1_0.%s' title='View the phylotree for this family' target='_blank'>%s</a>", family, family),
          sprintf("<a href='?job=%s&family=%s' title='Rebuild family phylotree including your sequence' target='%s.%s'><img src='tools-512.png' width='16px' height='16px' style='vertical-align: top'></a>", job$id, family, family, job$id)
        )
        gf2 <- df.hi[1, 1]
        gfs1 <- gfs2 <- df.hi[1, 5]
      }
      # gene families, gene family score (HTML and text format)
      c(gf1, gfs1, gf2, gfs2)
    })), stringsAsFactors = FALSE)
  }
  df.summary <- cbind(df.summary, data.frame(gf = df.h[, 1], gfs = df.h[, 2], stringsAsFactors = FALSE))
  df.summary.txt <- cbind(df.summary.txt, data.frame(gf = df.h[, 3], gfs = df.h[, 4], stringsAsFactors = FALSE))
  colnames.summary <- c(colnames.summary, "Gene Family", "GF Score<sup>1</sup>")

  # BLAST
  df.blast <- list()
  ii <- which(file.exists(job$blastFiles))
  for (i in ii) {
    df.blast[[i]] <- read.table(job$blastFiles[i], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  }
  # Loop over all sequences, match with first column of df.ahrd
  df.b <- as.data.frame(do.call(rbind, lapply(df.ahrd[, 1], function(q) {
    # Loop over the BLAST results and choose the one with the highest score (column 12)
    bbh <- bs <- score <- ""
    for (i in ii) {
      if (!(q %in% df.blast[[i]][, 1])) next
      df.bi <- df.blast[[i]][df.blast[[i]][, 1] == q, ]
      m <- which(df.bi[, 12] == max(df.bi[, 12]))[1]
      if (score == "" || df.bi[m, 12] > score) {
        score <- df.bi[m, 12]
        bs <- df.bi[m, 11]
        bbh <- df.bi[m, 2]
      }
    }
    # best BLAST hit, BLAST score
    c(bbh, bs)
  })), stringsAsFactors = FALSE)
  df.bl <- data.frame(bbh = df.b[, 1], bs = df.b[, 2], stringsAsFactors = FALSE)
  df.summary <- cbind(df.summary, df.bl)
  df.summary.txt <- cbind(df.summary.txt, df.bl)
  colnames.summary <- c(colnames.summary, "Best BLAST Hit", "BLAST Score<sup>2</sup>")

  # Sort by sequence name
  df.summary <- df.summary[order(df.summary[, 1]), ]
  df.summary.txt <- df.summary.txt[order(df.summary.txt[, 1]), ]

  if (!endsWith(job$summaryStatus, "Done")) {
    job$summaryStatus <- "Postprocessing: Done"
    writeJob(job)
  }

  list(simpleTable = df.simpleTable, columnNames = colnames.summary,
    summaryTable = df.summary, summaryTableOut = df.summary.txt)
}

# --------------------------------------------------------------

# Write sequences and their names for a given (job, gene family) to a temporary file for input to Lorax
# TODO: better description, refactoring
buildUserPhylogram <- function(job, family) {
  # TODO: error checking...

  # output: append user phylogram information or status/error messages, as appropriate
  userPhylogramInfo <- list(family = family, done = FALSE)
  df.summary.txt <- read.table(job$summaryFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  ff.matches <- (df.summary.txt$Gene.Family == family)
  userPhylogramInfo$descriptor <- df.summary.txt$AHRD.Descriptor[ff.matches]
  userPhylogramInfo$ipr <- df.summary.txt$InterPro.ID[ff.matches]
  userPhylogramInfo$go <- df.summary.txt$GO.Terms[ff.matches]

  # Check status of phylotree computation
  family_job <- paste(family, job$id, sep = ".")

  statusUrl <- sprintf("%s/trees/%s/FastTree/status", settings$lorax$url, family_job)
  statusResponse <- GET(statusUrl)
  statusResult <- 42
  tryCatch({
    statusResult <- fromJSON(rawToChar(statusResponse$content))
  }, error = function(e) {
    #statusResult <- 42
  })
  statusCode <- statusResponse$status_code

  if (statusResult == -1) {
    userPhylogramInfo$message <- "Computing phylogenetic tree, please be patient."
    return(userPhylogramInfo)

  } else if (statusResult == 0) {
    # Done computing phylogenetic tree, now parse it
    treeUrl <- sprintf("%s/trees/%s/FastTree/tree.nwk", settings$lorax$url, family_job)
    treeResponse <- GET(treeUrl)
    newickTree <- rawToChar(treeResponse$content)

    userSeqNames <- paste0("USR.", stri_match_all(newickTree, regex = sprintf("%s\\.([^\\:]+)\\:", job$id))[[1]][, 2])
    newickTree <- gsub(job$id, "USR", newickTree)

    # Return the multiple sequence alignment as well
    msaUrl <- sprintf("%s/trees/%s/alignment", settings$lorax$url, family_job)
    msaResponse <- GET(msaUrl)
    msa <- rawToChar(msaResponse$content)
    msa <- gsub(job$id, "USR", msa)

    # success
    userPhylogramInfo$seqNames <- userSeqNames
    userPhylogramInfo$tree <- trimws(newickTree)
    userPhylogramInfo$msa <- trimws(msa)
    userPhylogramInfo$done <- TRUE
    return(userPhylogramInfo)

  } else if (statusCode == 404) {
    # Phylogenetic tree computation not yet started

    # Find matching sequences for family, and write them to a temporary file to upload to Lorax
    seqNames <- df.summary.txt$Query[ff.matches]
    numMatchingSequences <- length(seqNames)
    if (numMatchingSequences == 0) {
      userPhylogramInfo$message <- paste("No matching sequences for", family)
      userPhylogramInfo$done <- TRUE
      return(userPhylogramInfo)
    }
    seqNames <- gsub("[|:]", ".", seqNames) # to prevent FASTA parser from treating them as delimiters
    # ---
    fasta <- readAAStringSet(job$inputFile) # original or translated protein sequence
    allSeqNames <- names(fasta)
    allSequences <- as.character(fasta)
    # ---
    seqFilename <- sprintf("temp/%s.fasta", family_job)
    if (file.exists(seqFilename)) {
      # TODO: do we need this?
      system(paste("rm", seqFilename))
      system(paste("touch", seqFilename))
    }
    for (sn in seqNames) {
      i <- which(grepl(sn, allSeqNames))
      write(paste0(">", allSeqNames[i]), seqFilename, append = TRUE)
      # split FASTA sequences into lines of at most 60 characters
      ss <- unlist(stri_match_all(allSequences[i], regex = ".{1,60}"))
      for (s in ss) write(s, seqFilename, append = TRUE)
    }

    # Upload (POST) matching sequences to Lorax
    sequencesUrl <- sprintf("%s/trees/%s/sequences", settings$lorax$url, family_job)
    sequencesResponse <- POST(sequencesUrl, body = list(peptide = upload_file(seqFilename)), verbose())
    sequencesCode <- sequencesResponse$status_code

    # Clean up
    if (file.exists(seqFilename)) system(paste("rm", seqFilename))

    if (sequencesCode == 500) {
      # Internal Server Error
      userPhylogramInfo$message <- sprintf("Error %d: Unable to compute tree for %s (no sequences).", sequencesCode, family)
      userPhylogramInfo$done <- TRUE
      return(userPhylogramInfo)
    } else if (sequencesCode != 200) {
      userPhylogramInfo$message <- sprintf("Error %d: Sequence upload for tree computation was not successful.", sequencesCode)
      userPhylogramInfo$done <- TRUE
      return(userPhylogramInfo)
    } else {
      # Tree computation using hmmalign
      hmmalignUrl <- sprintf("%s/trees/%s/hmmalign_FastTree", settings$lorax$url, family_job)
      hmmalignResponse <- GET(hmmalignUrl)
      hmmalignCode <- hmmalignResponse$status_code
      if (hmmalignCode != 200) {
        userPhylogramInfo$message <- sprintf("Error %d: Launch of tree computation was not successful.", hmmalignCode)
        userPhylogramInfo$done <- TRUE
      } else {
        userPhylogramInfo$message <- "Tree computation launched successfully."
      }
      return(userPhylogramInfo)
    }

  } else {
    userPhylogramInfo$message <- sprintf("Error %d: Unable to compute tree.", statusCode)
    userPhylogramInfo$done <- TRUE
    return(userPhylogramInfo)
  }
}

# --------------------------------------------------------------

