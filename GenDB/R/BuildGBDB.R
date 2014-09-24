BuildGBDB <- function(dataDir, ext="*.gbk", range=c(1,Inf), reg.exp=FALSE, ...)
{
#Check arguments
	if (length(range)==1) range <- c(1, range)

#Search for GenBank files
	cat("Searching for FASTA files ...\n")
	if (reg.exp) fasta.ext <- glob2rx(fasta.ext)
	files <- dir(dataDir, fasta.ext, recursive=TRUE, full.names=TRUE)
	first <- max(1, range[1])
	last <- min(c(length(files), range[2]))
	if (last<first)
    stop("the specified range has resulted in no genomes, so no database has been generated")
  files <- files[first:last]
 
#Get total bases in genomes
	totalBases <- 0
	
#Initialise database
  nEntries <- last-first+1
  db <- data.frame(
    Name=character(nEntries),
    Description=character(nEntries),
    Accession=character(nEntries),
    Version=character(nEntries),
    File=files[1:nEntries],
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
    genomeFile <- files[i]
    genome <- read.fasta(genomeFile)
    if (is.list(genome) && length(genome)>0) genome <- genome[[1]]
    if(is.null(attr(genome, "Annot"))) next #invalid genome load, skip it
    genomeInfo <- ParseFastaHeader(attr(genome, "Annot")) #extract FASTA header info and parse it
    db$Name[i] <- genomeInfo $Name
    db$Description[i] <- genomeInfo$Description
    db$Accession[i] <- genomeInfo$Accession
    db$Version[i] <- genomeInfo$Version
	  db$Length[i] <- length(genome)
	  db$UnknownBases[i] <- db$Length[i]-sum(genome %in% c("a", "c", "g", "t"))
    db$IsGenome[i] <- genomeInfo$IsGenome
    db$IsSequence[i] <- genomeInfo$IsSequence
    db$IsChromosome[i] <- genomeInfo$IsChromosome
    db$IsPlasmid[i] <- genomeInfo$IsPlasmid
    db$IsComplete[i] <- genomeInfo$IsComplete
    db$IsDraft[i] <- genomeInfo$IsDraft
    db$Registered <- TRUE
		elapsedTime = round((proc.time()-startTime)[3], 3)
	} #for i

	db <- structure(db, 	DataDir=dataDir, StartCreationDate=createTime, 
	  EndCreationDate=Sys.time(), class=c("SeqDB", "data.frame"))

	cat("Done.\n", sum(db$Registered), " of ", nEntries, " sequences successfully registered.\n")
	invisible(db) #return database
} #function
