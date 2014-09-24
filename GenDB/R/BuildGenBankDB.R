BuildGenBankDB <- function(dataDir, ext="*.gbk", range=c(1,Inf), reg.exp=FALSE, ...)
{
#Check arguments
	if (length(range)==1) range <- c(1, range)

#Search for GenBank files
	cat("Searching for GenBank files ...\n")
	if (reg.exp) ext <- glob2rx(ext)
	files <- dir(dataDir, ext, recursive=TRUE, full.names=TRUE)
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
    GI=character(nEntries),
    Content=character(nEntries),
    PhysicalStructure=character(nEntries),
    GB=character(nEntries),
    LastUpdate=character(nEntries),
    File=files,
    FileType="gbk",
    Length=integer(nEntries),
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
    file <- files[i]
    info <- read.gbk(file, entries=FALSE, sequence=FALSE)
    db$Name[i] <- info$Name
    db$Description[i] <- info$Description
    db$Organism[i] <- info$Organism
    db$Accession[i] <- info$Accession
    db$Version[i] <- info$Version
    db$GI[i] <- info$GI
    db$Content[i] <- info$Content
    db$PhysicalStructure[i] <- info$PhysicalStructure
    db$GB[i] <- info$GB
    db$LastUpdate[i] <- info$LastUpdate
	  db$Length[i] <- info$Length
    db$IsGenome[i] <- info$IsGenome
    db$IsSequence[i] <- info$IsSequence
    db$IsChromosome[i] <- info$IsChromosome
    db$IsPlasmid[i] <- info$IsPlasmid
    db$IsComplete[i] <- info$IsComplete
    db$IsDraft[i] <- info$IsDraft
    db$Registered[i] <- TRUE
		elapsedTime = round((proc.time()-startTime)[3], 3)
	} #for i

#Set class and other attributes of database
	attr(db, "DataDir") <- dataDir
	attr(db, "StartCreationDate") <- createTime
	attr(db, "EndCreationDate") <- Sys.time()
	class(db) <- c("SeqDB", "GenomeDB", "data.frame")

#	db <- structure(db, 	DataDir=dataDir, StartCreationDate=createTime, 
#	  EndCreationDate=Sys.time(), class=c("GenBankDB", "GenomeDB", "data.frame"))

	cat("Done.\n", sum(db$Registered), " of ", nEntries, " sequences successfully registered.\n")
	invisible(db) #return database
} #function
