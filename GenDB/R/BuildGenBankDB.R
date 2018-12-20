#BuildGenBankDB.R

BuildGenBankDB <- function(dataDir, ext="*.gbff.gz", range=c(1,Inf), reg.exp=FALSE, ...)
{
#Check arguments
	if (length(range)==1) range <- c(1, range)

#Search for GenBank files
	cat("Searching for GenBank files ...\n")
	if (reg.exp) ext <- glob2rx(ext)
	files <- dir(dataDir, ext, recursive=TRUE, full.names=FALSE)
	if (length(files)==0) stop("No files found")
	first <- max(1, range[1])
	last <- min(c(length(files), range[2]))
	if (last<first)
    stop("the specified range has resulted in no genomes, so no database has been generated")
  files <- files[first:last]
 
#Initialise database
  nFiles <- last-first+1
  db <- NULL #for holding database of GenBank record info

#Process genomes
  basesProcessed <- 0
  createTime <- Sys.time()
  startTime <- proc.time()
  elapsedTime <- 0
  for (fn in 1L:nFiles)
  {
    if (fn>1)
      cat(paste("Processing genome", fn, "of", nFiles, "...", elapsedTime, "s elapsed\n"))
    else
      cat(paste("Processing genome", fn, "of", nFiles, "...\n"))
    file <- file.path(dataDir, files[fn])
    info <- try(read.gbff(file, entries=FALSE, sequence=FALSE))
    if (class(info)=="try-erorr")
    {
      print(info)
      stop("read.gbff error")
    } #if
    if (is.null(info) || class(info)!="list" || any(sapply(info, class)!="gbk.info"))
    {
      cat("Invalid gbk.info:", file, "\n")
      next #invalid gbk info loaded, skip it
    } #if
    nRecords <- length(info)
    if (nRecords==0) continue
    seqInfo <- data.frame(
      Name=character(nRecords),
      Description="",
      Organism="",
      Domain="",
      Phylum="",
      Class="",
      Order="",
      Family="",
      Genus="",
      Species="",
      Strain="",
      Substrain="",
      taxid="",
      Accession="",
      Version="",
      GI="",
      Content="",
      PhysicalStructure="",
      GB="",
      LastUpdate="",
      File=files[fn],
      RecordNo=1L:nRecords,
      FileType="gbff",
      Length=0L,
      IsGenome=FALSE,
      IsSequence=FALSE,
      IsChromosome=FALSE,
      IsPlasmid=FALSE,
      IsComplete=FALSE,
      IsDraft=FALSE,
      Registered=FALSE,
      stringsAsFactors=FALSE
    )
    for (i in 1L:nRecords)
    {
      seqInfo$Name[i] <- info[[i]]$Name
      seqInfo$Description[i] <- info[[i]]$Description
      seqInfo$Organism[i] <- info[[i]]$BinomialName
      seqInfo$Domain[i] <- info[[i]]$Domain
      seqInfo$Phylum[i] <- info[[i]]$Phylum
      seqInfo$Class[i] <- info[[i]]$Class
      seqInfo$Order[i] <- info[[i]]$Order
      seqInfo$Family[i] <- info[[i]]$Family
      seqInfo$Genus[i] <- info[[i]]$Genus
      seqInfo$Species[i] <- info[[i]]$Species
      seqInfo$Strain[i] <- info[[i]]$Strain
      seqInfo$Substrain[i] <- info[[i]]$Substrain
      seqInfo$taxid[i] <- info[[i]]$taxid
      seqInfo$Accession[i] <- info[[i]]$Accession
      seqInfo$Version[i] <- info[[i]]$Version
      seqInfo$GI[i] <- info[[i]]$GI
      seqInfo$Content[i] <- info[[i]]$Content
      seqInfo$PhysicalStructure[i] <- info[[i]]$PhysicalStructure
      seqInfo$GB[i] <- info[[i]]$GB
      seqInfo$LastUpdate[i] <- info[[i]]$LastUpdate
	    seqInfo$Length[i] <- info[[i]]$Length
      seqInfo$IsGenome[i] <- info[[i]]$IsGenome
      seqInfo$IsSequence[i] <- info[[i]]$IsSequence
      seqInfo$IsChromosome[i] <- info[[i]]$IsChromosome
      seqInfo$IsPlasmid[i] <- info[[i]]$IsPlasmid
      seqInfo$IsComplete[i] <- info[[i]]$IsComplete
      seqInfo$IsDraft[i] <- info[[i]]$IsDraft
      seqInfo$Registered[i] <- TRUE
    } #for i
    db <- rbind(db, seqInfo)
		elapsedTime <- round((proc.time()-startTime)[3], 3)
	} #for fn

#Set class and other attributes of database
	attr(db, "DataDir") <- dataDir
	attr(db, "StartCreationDate") <- createTime
	attr(db, "EndCreationDate") <- Sys.time()
	class(db) <- c("GenBankDB", "GenomeDB", "data.frame")

#	db <- structure(db, 	DataDir=dataDir, StartCreationDate=createTime, 
#	  EndCreationDate=Sys.time(), class=c("GenBankDB", "GenomeDB", "data.frame"))

	cat("Done.\n", sum(db$Registered), "sequences out of", nrow(db), "from", nFiles, "assemblies successfully registered.\n")
	invisible(db) #return database
} #function
