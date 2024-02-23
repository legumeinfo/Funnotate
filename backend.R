# --------------------------------------------------------------
library(Biostrings)
library(httr)
library(InterMineR)
library(jsonlite)
library(stringi)
library(stringr)
library(yaml)
# --------------------------------------------------------------

# Global settings
# setwd("/srv/shiny-server/Funnotate") # already there
settings <- read_yaml("settings.yml")

# Global data:
# Read table of gene families
getGeneFamilies <- function() {
  gfFile <- tempfile()
  download.file(settings$gene_families_data, gfFile)
  df.gf <- read.table(gzfile(gfFile, "rt"), header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  unlink(gfFile)
  names(df.gf) <- c("name", "descriptor")
  df.gf$name <- gsub("-consensus", "", df.gf$name)
  df.gf
}
df.geneFamilies <- getGeneFamilies()

legfed_prefix <- "legfed_v1_0."

# LegumeMine service
legumeMine <- initInterMine(mine = listMines()["LegumeMine"])

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
  # assumes ss0 and ss1 have the same number of sequences (ns),
  # and each sequence has the same number of characters,
  # as ss1 = ss0 with characters replaced
  ns <- length(ss0) # == length(ss1)
  cc0 <- str_split(ss0, "")
  cc1 <- str_split(ss1, "")
  sum(sapply(1:ns, function(i) sum(cc0[[i]] != cc1[[i]]) ))
}

isProbablyNucleotideSequence <- function(lines) {
  ss <- lines[!startsWith(lines, ">")]
  ss <- tolower(str_split_1(str_flatten(ss), ""))
  nn <- str_split_1(scrubber$n$goodChars, "")
  all(ss %in% nn)
}

# Create directory if it does not exist
requireDirectory <- function(dirpath) {
  if (!dir.exists(dirpath)) dir.create(dirpath)
}

# File exists and is not empty
fileReallyExists <- function(filename) {
  fileExists <- file.exists(filename)
  result <- (fileExists && file.info(filename)$size > 0)
  if (fileExists && !result) file.remove(filename)
  result
}
hmmFileReallyExists <- function(filename) {
  fileExists <- fileReallyExists(filename)
  if (!fileExists) return(FALSE)
  ll <- readLines(filename)
  ll <- ll[!startsWith(ll, "#")]
  result <- (length(ll) > 0)
  if (!result) file.remove(filename)
  result
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

logError <- function(errMsg) {
  errFile <- "static/error.log"
  write(currentTimeString(), errFile, append = TRUE)
  write(errMsg, errFile, append = TRUE)
  write("----------------------------------------", errFile, append = TRUE)
}

# --------------------------------------------------------------

# fiData =
#   1. For uploaded FASTA files, the value of a Shiny fileInput control (name, size, datapath).
#   2. For FASTA sequences pasted into the text input, (name = source, size, datapath = temporary file).
#   3. For FASTA sequences posted from InterMine, (name = source, size, seqNames, sequences)
#      to avoid creating the upload files until necessary.
# seqType = n for nucleotide/DNA, p for protein sequence
createNewUpload <- function(fiData, seqType) {
  scrub <- scrubber[[seqType]]

  uploadDir <- "static/upload"
  requireDirectory(uploadDir)
  indexFile <- sprintf("%s/index", uploadDir)
  if (!file.exists(indexFile)) {
    write("1", indexFile)
  }
  nextIndex <- scan(indexFile, what = integer())

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

  if (!is.null(fiData$datapath)) {
    # Read the sequences (to a Biostrings::DNAStringSet or Biostrings::AAStringSet)
    if (seqType == "n") {
      fasta <- readDNAStringSet(fiData$datapath)
    } else {
      fasta <- readAAStringSet(fiData$datapath)
    }
    seqNames <- names(fasta)
    sequences <- as.character(fasta)
  } else {
    # We already extracted the sequences from the POST
    seqNames <- fiData$seqNames
    sequences <- fiData$sequences
  }
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

  write(as.character(nextIndex + 1), indexFile)
  write_yaml(upload, upload$uploadFile)
  upload
}

# --------------------------------------------------------------

createNewJob <- function(upload, useInterpro) {
  # Generate an as-yet-unused job id
  jobsDir <- "static/job"
  requireDirectory(jobsDir)
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

createNewJobWithGeneFamily <- function(upload, geneFamily) {
  # Generate an as-yet-unused job id
  jobsDir <- "static/job"
  requireDirectory(jobsDir)
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
    inputFile = upload$inputFileScrubbed,
    geneFamily = geneFamily,
    jobFile = sprintf("%s/job_%s", jobDir, jobId), # metadata for existing jobs
    summaryFile = sprintf("%s/summary_%s.txt", jobDir, jobId),
    status = "new",
    startTime = currentTimeString(),
    endTime = ""
  )
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
  ahrdTmpYmlFile <- tempfile()
  write_yaml(ahrdYml, ahrdTmpYmlFile)
  ahrdCmd <- sprintf("%s -Xmx2g -jar %s %s", settings$ahrd$java, settings$ahrd$jar, ahrdTmpYmlFile)
  system(ahrdCmd)
  # clean up temporary file
  if (!is.na(ahrdTmpYmlFile) && file.exists(ahrdTmpYmlFile)) unlink(ahrdTmpYmlFile)
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
  iprXml <- tempfile()
  iprCmdXml <- sprintf("%s -i %s -o %s -f XML %s", settings$interpro$exe, job$inputFile, iprXml, settings$interpro$params)
  job$iprStatus <- "InterPro: Running"
  writeJob(job)
  system(iprCmdXml)
  if (file.exists(iprXml)) {
    iprCmdRaw <- sprintf("%s -i %s -mode convert -f RAW -o %s", settings$interpro$exe, iprXml, job$iprFile)
    system(iprCmdRaw)
    # clean up temporary file
    if (!is.na(iprXml) && file.exists(iprXml)) unlink(iprXml)
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
  requireDirectory(job$dir)

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

runJobWithGeneFamily <- function(job) {
  requireDirectory(job$dir)

  # ESTScan
  if (job$sequenceType == "nucleotide") {
    if (isActive(job)) job <- runESTScan(job)
  }

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
  uploadFile <- sprintf("static/upload/upload_%d", index)
  if (!file.exists(uploadFile)) return(NULL)
  upload <- read_yaml(uploadFile)
  upload
}

# Read an existing job
readJob <- function(jobId) {
  jobFile <- sprintf("static/job/%s/job_%s", jobId, jobId)
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
  df.ahrd <- read.table(job$ahrdFile, skip = 2, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
  df.summary <- df.summary.txt <- df.ahrd[, 1:4]
  colnames.summary <- c("Query", "AHRD BLAST Hit", "AHRD Quality<sup>3</sup>", "AHRD Descriptor")
  blankColumn <- rep("", nrow(df.summary))

  # InterPro
  if (job$useInterpro) {
    if (!file.exists(job$iprFile)) {
      df.i <- data.frame(iit1 = blankColumn, go1 = blankColumn, iit2 = blankColumn, go2 = blankColumn, stringsAsFactors = FALSE)
    } else {
      df.ipr <- read.table(job$iprFile, header = FALSE, sep = "\t", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
      hasGOTermColumn <- (ncol(df.ipr) >= 14)
      # Loop over all sequences, match with first column of df.ahrd
      df.i <- as.data.frame(do.call(rbind, lapply(df.ahrd[, 1], function(q) {
        df.ii <- df.ipr[df.ipr[, 1] == q, ]
        iit <- setdiff(unique(df.ii[, 12]), "NULL")
        iit1 <- ifelse(length(iit) == 0, "", paste(sprintf("<a href='https://www.ebi.ac.uk/interpro/entry/%s' target='_blank'>%s</a>", iit, iit), collapse = ", "))
        iit2 <- ifelse(length(iit) == 0, "", paste(iit, collapse = ","))
        if (hasGOTermColumn) {
          go <- stri_match_all(df.ii[, 14], regex = "\\(GO:(\\d+)\\)")
          go <- setdiff(unique(unlist(sapply(go, function(g) g[, 2], USE.NAMES = FALSE))), NA)
          go1 <- ifelse(length(go) == 0, "", paste(sprintf("<a href='http://amigo.geneontology.org/amigo/term/GO:%s' target='_blank'>%s</a>", go, go), collapse = ", "))
          go2 <- ifelse(length(go) == 0, "", paste(go, collapse = ","))
        } else {
          go1 <- go2 <- ""
        }
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
    df.hmm <- as.data.frame(do.call(rbind, str_split(ll.hmm, "\\s+")), stringsAsFactors = FALSE)
    # Loop over all sequences, match with first column of df.ahrd
    df.h <- as.data.frame(do.call(rbind, lapply(df.ahrd[, 1], function(q) {
      df.hi <- df.hmm[df.hmm[, 3] == q, ]
      if (nrow(df.hi) == 0) {
        gf1 <- gfs1 <- gf2 <- gfs2 <- ""
      } else {
        family <- df.hi[1, 1]
        gf1 <- paste(
          sprintf("<a href='?family=%s' title='View the phylotree for this family' target='_blank'>%s</a>", family, family),
          sprintf("<a href='?job=%s&family=%s' title='Rebuild family phylotree including your sequence' target='%s.%s'><img src='static/tools-512.png' width='16px' height='16px' style='vertical-align: top'></a>", job$id, family, family, job$id)
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
    df.blast[[i]] <- read.table(job$blastFiles[i], header = FALSE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
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

# Given an MSA and a Newick tree with the same sequences,
# return the MSA rearranged into the same order as in the tree
msaOrderedLikeTree <- function(msa_in, tree) {
  # sequence names from the tree (top to bottom)
  tree_names <- stri_match_all(tree, regex = "[\\(\\,]([^\\(\\,\\:]+)\\:")[[1]][, 2]

  # msa_in is a single string with rows separated by newlines,
  # so convert it to individual rows (like a FASTA file)
  msa1 <- str_split_1(msa_in, "\n")
  nr <- length(msa1) # total number of rows (lines)
  ss <- which(startsWith(msa1, ">")) # indices of rows containing sequence name (= first row of each sequence)
  ee <- c(ss[-1] - 1, nr) # indices of last row of each sequence

  # put the sequences in the correct order (that of the tree)
  oo <- order(match(substring(msa1[ss], 2), tree_names))
  msa2 <- character(nr)
  l <- 0
  for (o in oo) {
    nl <- ee[o] - ss[o] + 1 # number of lines in the oth sequence
    msa2[l + 1:nl] <- msa1[ss[o]:ee[o]]
    l <- l + nl
  }

  # combine all rows into a single string (separated by newlines)
  msa_out <- paste(msa2, collapse = "\n")
  msa_out
}

# Write sequences and their names for a given (job, gene family) to a temporary file for input to Lorax
# TODO: better description, refactoring
buildUserPhylogram <- function(job, family) {
  # TODO: error checking...

  # output: append user phylogram information or status/error messages, as appropriate
  userPhylogramInfo <- list(family = family, done = FALSE)
  i.match <- which(grepl(family, df.geneFamilies$name))
  userPhylogramInfo$descriptor <- ifelse(length(i.match) == 0, "unknown", df.geneFamilies$descriptor[i.match])

  if (is.null(job)) {
    # Display precomputed phylotree and MSA for the given family
    treeUrl <- sprintf("%s/trees/%s/FastTree/tree.nwk", settings$lorax$url, family)
    treeResponse <- GET(treeUrl)
    msaUrl <- sprintf("%s/trees/%s/alignment", settings$lorax$url, family)
    msaResponse <- GET(msaUrl)
    # failure
    if (treeResponse$status_code != 200 || msaResponse$status_code != 200) return(NULL)
    # success
    newickTree <- rawToChar(treeResponse$content)
    msa <- rawToChar(msaResponse$content)
    userPhylogramInfo$tree <- trimws(newickTree)
    userPhylogramInfo$msa <- msaOrderedLikeTree(trimws(msa), userPhylogramInfo$tree)
    userPhylogramInfo$done <- TRUE
    return(userPhylogramInfo)
  }

  # Check status of phylotree computation
  # (wrap in tryCatch() to detect errors arising in Lorax)
  tryCatch({
    family_job <- paste(family, job$id, sep = ".")
    statusUrl <- sprintf("%s/trees/%s/FastTree/status", settings$lorax$url, family_job)
    statusResponse <- GET(statusUrl)
    statusResult <- 42
    tryCatch({
      statusResult <- fromJSON(rawToChar(statusResponse$content))
    }, error = function(e2) {
      # ...
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
      userPhylogramInfo$msa <- msaOrderedLikeTree(trimws(msa), userPhylogramInfo$tree)
      userPhylogramInfo$done <- TRUE
      return(userPhylogramInfo)

    } else if (statusCode == 404) {
      # Phylogenetic tree computation not yet started

      # Find matching sequences for family
      df.summary.txt <- read.table(job$summaryFile, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
      ff.matches <- (df.summary.txt$Gene.Family == family)
      seqNames <- df.summary.txt$Query[ff.matches]
      numMatchingSequences <- length(seqNames)
      # If there are no matches, let the user know (and return)
      if (numMatchingSequences == 0) {
        userPhylogramInfo$message <- paste("No matching sequences for", family)
        userPhylogramInfo$done <- TRUE
        logError(sprintf("%s (job %s)", userPhylogramInfo$message, job$id))
        return(userPhylogramInfo)
      }

      # Read the user's original (or translated) protein sequences
      fasta <- readAAStringSet(job$inputFile)
      allSeqNames <- names(fasta)
      # To correctly match user sequences by name,
      # remove any characters after and including the first whitespace
      allSeqNames <- sapply(allSeqNames, function(sn) str_split_i(sn, " ", 1))
      names(fasta) <- allSeqNames
      allSequences <- as.character(fasta)

      # Write matching sequences to a temporary file to upload to Lorax
      seqFile <- tempfile()
      file.create(seqFile)
      for (sn in seqNames) {
        sn_out <- gsub("[|:]", ".", sn) # to prevent downstream Newick tree parser from treating them as delimiters
        write(paste0(">", sn_out), seqFile, append = TRUE)
        # split FASTA sequences into lines of at most 60 characters
        ss <- unlist(stri_match_all(allSequences[sn], regex = ".{1,60}"))
        for (s in ss) write(s, seqFile, append = TRUE)
      }

      # Upload (POST) matching sequences to Lorax
      sequencesUrl <- sprintf("%s/trees/%s/sequences", settings$lorax$url, family_job)
      sequencesResponse <- POST(sequencesUrl, body = list(peptide = upload_file(seqFile)), verbose())
      sequencesCode <- sequencesResponse$status_code

      # Clean up
      if (file.exists(seqFile)) unlink(seqFile)

      if (sequencesCode == 500) {
        # Internal Server Error
        userPhylogramInfo$message <- sprintf("Error %d: Unable to compute tree for %s (no sequences).", sequencesCode, family)
        userPhylogramInfo$done <- TRUE
        logError(sprintf("%s (job %s)", userPhylogramInfo$message, job$id))
        return(userPhylogramInfo)
      } else if (sequencesCode != 200) {
        userPhylogramInfo$message <- sprintf("Error %d: Sequence upload for tree computation was not successful.", sequencesCode)
        userPhylogramInfo$done <- TRUE
        logError(sprintf("%s (family %s, job %s)", userPhylogramInfo$message, family, job$id))
        return(userPhylogramInfo)
      } else {
        # Tree computation using hmmalign
        hmmalignUrl <- sprintf("%s/trees/%s/hmmalign_FastTree", settings$lorax$url, family_job)
        hmmalignResponse <- GET(hmmalignUrl)
        hmmalignCode <- hmmalignResponse$status_code
        if (hmmalignCode != 200) {
          userPhylogramInfo$message <- sprintf("Error %d: Launch of tree computation was not successful.", hmmalignCode)
          userPhylogramInfo$done <- TRUE
          logError(sprintf("%s (family %s, job %s)", userPhylogramInfo$message, family, job$id))
        } else {
          userPhylogramInfo$message <- "Tree computation launched successfully."
        }
        return(userPhylogramInfo)
      }

    } else {
      userPhylogramInfo$message <- sprintf("Error %d: Unable to compute tree.", statusCode)
      userPhylogramInfo$done <- TRUE
      logError(sprintf("%s (family %s, job %s)", userPhylogramInfo$message, family, job$id))
      return(userPhylogramInfo)
    }
  }, warning = function(w) {
    # Report Lorax warnings (note that w alone returns a less specific message)
    userPhylogramInfo$message <- paste(w, "<br>The Funnotate sysadmin has been notified.")
    userPhylogramInfo$done <- TRUE
    logError(sprintf("%s (family %s, job %s)", trimws(as.character(w)), family, job$id))
    return(userPhylogramInfo)
  }, error = function(e) {
    # Report Lorax errors (note that e alone returns a less specific message)
    userPhylogramInfo$message <- paste(e, "<br>The Funnotate sysadmin has been notified. Please try again later.")
    userPhylogramInfo$done <- TRUE
    logError(sprintf("%s (family %s, job %s)", trimws(as.character(e)), family, job$id))
  })
}

# --------------------------------------------------------------

geneFamilySearchQuery <- function(keywords) {
  keywords <- trimws(keywords)
  if (nchar(keywords) == 0) return(NULL)

  familyRequest <- sprintf("%s/service/search?q=%s&searchBag=&facet_Category=GeneFamily", legumeMine@mine, URLencode(keywords))
  familyResponse <- GET(familyRequest)
  json <- fromJSON(rawToChar(familyResponse$content))
  results_per_page <- 100 # set somewhere in LegumeMine
  num_pages <- 1 + (json$totalHits - 1) %/% results_per_page
  familyResults <- json$results$fields
  familyResults$relevance <- json$results$relevance
  if (num_pages > 1) {
    for (p in 2:num_pages) {
      s <- (p - 1)*results_per_page # start field in URL must have zero offset
      familyResponse_p <- GET(paste0(familyRequest, "&start=", s))
      json_p <- fromJSON(rawToChar(familyResponse_p$content))
      familyResults_p <- json_p$results$fields
      familyResults_p$relevance <- json_p$results$relevance
      familyResults <- rbind(familyResults, familyResults_p)
    }
  }

  if (!is.null(familyResults)) {
    families <- stri_match_first(familyResults$primaryIdentifier, regex = "^.+\\.(.+)$")[, 2]
    familyResults$primaryIdentifier <- sprintf("<a href='?family=%s' target='_blank'>%s</a>", families, familyResults$primaryIdentifier)
  }
  familyResults
}

# --------------------------------------------------------------

genesToProteinsQuery <- function(family, genes) {
  # add prefix to restore full-yuck gene family name for the query
  family <- paste0(legfed_prefix, family)

  # convert genes to character vector
  genes <- str_split_1(URLdecode(genes), ",")

  # find all genes in family
  family_constraints = setConstraints(
    paths = "GeneFamily.primaryIdentifier",
    operators = "=",
    values = list(family)
  )
  genes_query = setQuery(
    select = c("GeneFamily.genes.primaryIdentifier"),
    where = family_constraints
  )
  all_genes <- runQuery(legumeMine, genes_query)

  # match against user-supplied genes
  matched_genes <- base::intersect(all_genes$GeneFamily.genes.primaryIdentifier, genes)
  if (length(matched_genes) == 0) return(NULL)

  # match against proteins
  gene_constraints = setConstraints(
    paths = "Gene.primaryIdentifier",
    operators = "=",
    values = list(matched_genes)
  )
  protein_query = setQuery(
    select = c("Gene.proteins.primaryIdentifier"),
    where = gene_constraints
  )
  proteins <- runQuery(legumeMine, protein_query)
  if (length(proteins) == 0) return(NULL)
  proteins$Gene.proteins.primaryIdentifier
}

# --------------------------------------------------------------

genesToGeneFamiliesQuery <- function(genes) {
  # convert genes to character vector
  genes <- str_split_1(URLdecode(genes), ",")

  # find gene families associated with genes
  gene_constraints = setConstraints(
    paths = "Gene.primaryIdentifier",
    operators = "=",
    values = list(genes)
  )
  gene_families_query = setQuery(
    select = c(
      "Gene.primaryIdentifier",
      "Gene.geneFamilyAssignments.geneFamily.primaryIdentifier",
      "Gene.geneFamilyAssignments.geneFamily.description"
    ),
    where = gene_constraints
  )
  query_results <- runQuery(legumeMine, gene_families_query) # data frame if successful, list of length 0 if not
  if (length(query_results) == 0) {
    return(data.frame(geneFamily = "", description = "Not found",
      genes = paste0(genes, collapse = "<br>"), gene_name = "", link = ""))
  }

  # rearrange query_results to be by gene family instead of by gene
  names(query_results) <- c("gene", "geneFamily", "description")
  query_results <- query_results[startsWith(query_results$geneFamily, legfed_prefix), ]
  if (nrow(query_results) == 0) {
    return(data.frame(geneFamily = "", description = "Not found",
      genes = paste0(genes, collapse = "<br>"), gene_name = "", link = ""))
  }
  ff <- sort(unique(query_results$geneFamily))
  dd <- sapply(ff, function(f) query_results$description[min(which(query_results$geneFamily == f))], USE.NAMES = FALSE)
  gg_comma <- sapply(ff, function(f) paste(query_results$gene[which(query_results$geneFamily == f)], collapse = ","), USE.NAMES = FALSE)
  gg_br <- gsub(",", "<br>", gg_comma)
  df_gene_families <- data.frame(geneFamily = ff, description = dd, genes = gg_br, gene_name = gg_comma,
    link = sprintf("<a href='?family=%s&gene_name=%s'>%s</a>", ff, gg_comma, ff))

  # add non-family row for any genes not found
  gg_found <- str_split_1(gg_comma, ",")
  gg_not_found <- setdiff(genes, gg_found)
  if (length(gg_not_found) > 0) {
    df_gene_families <- rbind(df_gene_families, data.frame(geneFamily = "", description = "Not found",
      genes = paste0(gg_not_found, collapse = "<br>"), gene_name = "", link = ""))
  }

  df_gene_families
}

# --------------------------------------------------------------

