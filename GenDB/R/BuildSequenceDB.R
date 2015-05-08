#BuildSequenceDB

BuildSequenceDB <- function(dataDir, ext="*.fna", range=c(1,Inf), reg.exp=FALSE, ...)
{
#Check arguments
	if (length(range)==1) range <- c(1, range)

#Search for genome files
	cat("Searching for FASTA files ...\n")
	if (reg.exp) fasta.ext <- glob2rx(ext)
	files <- dir(dataDir, ext, recursive=TRUE, full.names=FALSE)
	if (length(files)==0) stop("No files found")
	first <- max(1, range[1])
	last <- min(c(length(files), range[2]))
	if (last<first)
    stop("the specified range has resulted in no genomes, so no database has been generated")
  files <- files[first:last]
 
#Initialise database
  nEntries <- last-first+1
  db <- data.frame(
    Name=character(nEntries),
    Description=character(nEntries),
    Organism=character(nEntries),
    Accession=character(nEntries),
    Version=character(nEntries),
    File=files,
    FileType="fna",
    Length=integer(nEntries),
    AmbiguousBases=integer(nEntries),
    IsGenome=logical(nEntries),
    IsSequence=logical(nEntries),
    IsChromosome=logical(nEntries),
    IsPlasmid=logical(nEntries),
    IsComplete=logical(nEntries),
    IsDraft=logical(nEntries),
    Registered=logical(nEntries),
    stringsAsFactors=FALSE
  )
  rownames(db) <- as.character(first:last)

#Process genomes
  basesProcessed <- 0
  createTime <- Sys.time()
  startTime <- proc.time()
  for (i in 1:nEntries)
  {
    if (i>1)
      cat(paste("Processing genome", i, "of", nEntries, "...", elapsedTime, "s elapsed\n"))
    else
      cat(paste("Processing genome", i, "of", nEntries, "...\n"))
    file <- file.path(dataDir, files[i])
    genome <- NULL
    res <- try(genome <- read.fasta(file))
    if (class(res)=="try-erorr")
    {
      print(res)
      stop("try error bail out")
    } #if
    if (is.null(genome))
    {
      cat("Invalid FASTA:", file, "\n")
      next #invalid genome load, skip it
    } #if
    if (is.list(genome) && length(genome)>0) genome <- genome[[1]]
    if(is.null(attr(genome, "Annot")))
    {
      cat("Invalid Annot:", file, "\n")
      next #invalid genome load, skip it
    } #if
    info <- ParseFastaHeader(attr(genome, "Annot")) #extract FASTA header info and parse it
    db$Name[i] <- info $Name
    db$Description[i] <- info$Description
    db$Organism[i] <- info$Organism
    db$Accession[i] <- info$Accession
    db$Version[i] <- info$Version
	  db$Length[i] <- length(genome)
	  db$AmbiguousBases[i] <- db$Length[i]-sum(genome %in% c("a", "c", "g", "t"))
    db$IsGenome[i] <- info$IsGenome
    db$IsSequence[i] <- info$IsSequence
    db$IsChromosome[i] <- info$IsChromosome
    db$IsPlasmid[i] <- info$IsPlasmid
    db$IsComplete[i] <- info$IsComplete
    db$IsDraft[i] <- info$IsDraft
    db$Registered <- TRUE
		elapsedTime = round((proc.time()-startTime)[3], 3)
	} #for i

#Set class and other attributes of database
	attr(db, "DataDir") <- dataDir
	attr(db, "StartCreationDate") <- createTime
	attr(db, "EndCreationDate") <- Sys.time()
	class(db) <- c("SeqDB", "GenomeDB", "data.frame")

	cat("Done.\n", sum(db$Registered), " of ", nEntries, " sequences successfully registered.\n")
	invisible(db) #return database
} #function


ParseFastaHeader <- function(header)
{
	if (is.null(header)) return(list())
	tokens <- regmatches(header, regexec("^.*\\|.*\\|(.*)\\|(.*)\\|(.*)$", header))[[1]]
	definition <- tokens[4]
	info <- list(Source=tokens[1], Definition=definition, Accession=sub("\\.\\d+$", "", tokens[3]), Version=tokens[3])
	info <- ParseGenBankDefinition(definition, info)
	info 
} #function
