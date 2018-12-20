#BuildSequenceDB

BuildSequenceDB <- function(dataDir, ext="*.fna.gz", range=c(1,Inf), 
reg.exp=FALSE, save=TRUE, ...)
{
#Check arguments
	if (length(range)==1) range <- c(1, range)

#Search for genome files
	cat("Searching for FASTA files ...\n")
	if (reg.exp) ext <- glob2rx(ext)
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
    #Registered=logical(nEntries),
    stringsAsFactors=FALSE
  )
  rownames(db) <- as.character(first:last)

#Process genomes
  basesProcessed <- 0
  createTime <- Sys.time()
  startTime <- proc.time()
db <- data.frame()
  for (i in 1:nEntries)
  {
    if (i>1)
      cat(paste("Processing genome", i, "of", nEntries, "...", elapsedTime, "s elapsed\n"))
    else
      cat(paste("Processing genome", i, "of", nEntries, "...\n"))
    file <- file.path(dataDir, files[i])
    seqs <- NULL
    res <- try(seqs <- read.fasta(file))
    if (class(res)=="try-erorr")
    {
      print(res)
      stop("try error bail out")
    } #if
    if (is.null(seqs) || length(seqs)==0)
    {
      cat("Invalid FASTA:", file, "\n")
      next #invalid genome load, skip it
    } #if
    for (j in seq_along(seqs)) {
    	genome <- seqs[[j]]
    	if(is.null(attr(genome, "Annot")))
    	{
      	cat("Invalid Annot:", file, "\n")
      	next #invalid genome load, skip it
    	} #if
    	info <- ParseFastaHeader(attr(genome, "Annot")) #extract FASTA header info and parse it
    	seqInfo <- list(
    		Name=info$Name,
    		Description=info$Description,
    		Organism=info$BinomialName,
    		Accession=info$Accession,
    		Version=info$Version,
    		File=files[i],
    		FileType="fasta",
    		FileOffset=j,
	  		Length=length(genome),
	  		AmbiguousBases=length(genome)-sum(!is.na(genome) & genome %in% c("a", "c", "g", "t")),
    		IsGenome=info$IsGenome,
    	IsSequence=info$IsSequence,
    	IsChromosome=info$IsChromosome,
    		IsPlasmid=info$IsPlasmid,
    		IsComplete=info$IsComplete,
    		IsDraft=info$IsDraft
    	)
    	db <- rbind(db, seqInfo, stringsAsFactors=FALSE)
    } #for j
		elapsedTime = round((proc.time()-startTime)[3], 3)
	} #for i

#Set class and other attributes of database
	attr(db, "DataDir") <- dataDir
	attr(db, "StartCreationDate") <- createTime
	attr(db, "EndCreationDate") <- Sys.time()
	class(db) <- c("SeqDB", "GenomeDB", "data.frame")

	if (save) save(db, file=file.path(dataDir, "db.rda"))
	cat("Done.\n", nrow(db), "sequences in", nEntries, "files successfully registered.\n")
	invisible(db) #return database
} #function


ParseFastaHeader <- function(header)
{
	if (is.null(header)) return(list())
#	tokens <- regmatches(header, regexec("^.*\\|.*\\|(.*)\\|(.*)\\|(.*)$", header))[[1]]
	tokens <- strsplit(header, "\\|")[[1]]
	if (length(tokens)==0) return(list())
  if (length(tokens)==4) {
	  definition <- tokens[4]
	    source <- sub("^>", "", tokens[1])
	    version <- tokens[3]
	} else {
		tokens <- tokens[1]
	  source <- ""
	  def <- regmatches(tokens, regexec("^>?([^ ]+) +(.*)$", tokens))[[1]]
	  if (length(def)<2) return(list())
	  version <- def[2]
	  if (length(def)==3) {
	    definition <- def[3]
	  } else {
	    definition <- ""
	  }
	}
	info <- list(Source=source, Definition=definition, Accession=sub("\\.\\d+$", "", version), Version=version)
	info <- ParseGenBankDefinition(definition, info)
	info 
} #function
