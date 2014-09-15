library(seqinr)

#BuildSequenceDB
#Compile a database of genomes from a directory structure containing
#FASTA files.
#
#USAGE
#BuildSequenceDB(dataDir, limit=c(1,Inf), filter="*.fna")
#
#INPUTS
#dataDir:  The source directory at which to start a recursive search for 
#    FASTA files.
#
#OPTIONAL INPUTS
#limit:  Range of genomes to process.  If not specified or left empty, 
#    it defaults to c(1, Inf), which means process all genomes.  This option 
#    is only meant for speed testing and will not generate a proper database.
#filter:  File extension used for identifying files to process.  The default 
#    value is ".fna".
#
#FUNCTION DEPENDENCIES
#GetFileList, readfastafile, dnaread

#*function
BuildSequenceDB <- function(dataDir, limit=c(1,Inf), filter="*.fna")
{
#Check arguments
	if (length(limit)==1) limit <- c(1, limit)

#Search for genome files
	cat("Searching for FASTA files ...\n")
	files = dir(dataDir, glob2rx(filter), recursive=TRUE, full.names=TRUE)
	first <- max(1, limit[1])
	last <- min(c(length(files), limit[2]))
	if (last<=0)
  stop("No genomes, so no database generated.")

#Get total bases in genomes
	totalBases <- 0
	#for (i in first:last)
	#{
	#	totalBases <- totalBases + files{i}.filesize;
	#} #for i

#Check for destination directory
#	if (!dir.exists(dbDir))
#	{
#		cat("Database directory does not exist, so creating it.\n")
#		dir.create(dbDir)
#	} #if

#Initialise database
db <- data.frame(
	Name=character(last),
	Description=character(last),
	Accession=character(last),
	File=files[1:last],
	Length=integer(last),
	UnknownBases=integer(last),
	IsGenome=logical(last),
	IsSequence=logical(last),
	IsChromosome=logical(last),
	IsPlasmid=logical(last),
	IsComplete=logical(last),
	IsDraft=logical(last),
	stringsAsFactors=FALSE
)

#Process genomes
	basesProcessed <- 0
	createTime <- Sys.time()
	startTime <- proc.time()
	for (i in first:last)
	{
		if (i>first)
			cat(paste("Processing genome", i, "of", last, "...", elapsedTime, "s elapsed\n"))
  	else
			cat(paste("Processing genome", i, "of", last, "...\n"))
  genomeFile <- files[i]
  genome <- read.fasta(genomeFile)
  if (is.list(genome) && length(genome)>=1) genome <- genome[[1]]
  genomeInfo <- GetFastaHeaderInfo(genome)
  db$Name[i] <- genomeInfo $Name
  db$Description[i] <- genomeInfo$Description
  db$Accession[i] <- genomeInfo$Accession
	db$Length[i] <- length(genome)
	db$UnknownBases[i] <- db$Length[i]-sum(genome %in% c("a", "c", "g", "t"))
  db$IsGenome[i] <- genomeInfo$IsGenome
  db$IsSequence[i] <- genomeInfo$IsSequence
  db$IsChromosome[i] <- genomeInfo$IsChromosome
  db$IsPlasmid[i] <- genomeInfo$IsPlasmid
  db$IsComplete[i] <- genomeInfo$IsComplete
  db$IsDraft[i] <- genomeInfo$IsDraft
		elapsedTime = round((proc.time()-startTime)[3], 3)
	} #for i

	db <- structure(db, 	DataDir=dataDir, StartCreationDate=createTime, EndCreationDate=Sys.time())
	class(db) <- c("SeqDB", "data.frame")

	cat("Done.\n")
	return(db)
} #function

#*function
dir.exists <- function(...)
{
	file.exists(...) & file.info(...)$isdir
} #function

#*function
GetFastaHeaderInfo.originalcode <- function(fasta.data)
{
	header <- attr(fasta.data, "Annot")
	tokens <- regmatches(header, regexec("^.*\\|.*\\|(.*)\\|(.*)\\|\\s*(.*),\\s*(.*)$", header))[[1]]
	if (length(tokens)==0) tokens <- regmatches(header, regexec("^.*\\|.*\\|(.*)\\|(.*)\\|\\s*(.*)$", header))[[1]]
	tokens <- tokens[-1]
	genome <- list(Source=tokens[1], Accession=tokens[2], Name=tokens[3])
	if (length(tokens)>=4)
	{
		genome$Description <- tokens[4]
		genome$IsGenome <- grepl("genome", genome$Description, ignore.case=TRUE)
		genome$IsSequence <- grepl("sequence", genome$Description, ignore.case=TRUE)
		genome$IsChromosome <- grepl("chromosome", genome$Name, ignore.case=TRUE)
		genome$IsPlasmid <- grepl("plasmid", genome$Name, ignore.case=TRUE)
		genome$IsComplete <- grepl("complete", genome$Description, ignore.case=TRUE)
	} else {
		genome$Description <- ""
		genome$IsGenome <- FALSE
		genome$IsSequence <- FALSE
		genome$IsChromosome <- FALSE
		genome$IsPlasmid <- FALSE
		genome$IsComplete <- FALSE
	} #if
	genome$IsDraft <- any(grepl("draft", c(genome$Name, genome$Description), ignore.case=TRUE))
	return(genome)
} #function

#*function
GetFastaHeaderInfo <- function(fasta.data)
{
	header <- attr(fasta.data, "Annot")
	if (is.null(header)) return(list())
	tokens <- regmatches(header, regexec("^.*\\|.*\\|(.*)\\|(.*)\\|(.*)$", header))[[1]]
	definition <- tokens[4]
	info <- list(Source=tokens[1], Definition=definition, Accession=sub("\\.\\d+$", "", tokens[3]), Version=tokens[3])
	info <- ParseGenBankDefinition(definition, info)
	info 
} #function

#*function
ParseGenBankDefinition <- function(definition, info=list())
{
	if (missing(definition) || !is.character(definition))
		stop("definition must be of type character")
	if (!is.list(info)) stop("info must of a list")
	tokens <- regmatches(definition, regexec("^\\s*(.*),\\s*(.*)$", definition))[[1]]
	if (length(tokens)==0) tokens <- regmatches(definition, regexec("^\\s*(.*)\\s*$", definition))[[1]]
	tokens <- tokens[-1]
	info$Name <- tokens[1]
	if (length(tokens)>=2)
	{
		info$Description <- tokens[2]
		info$IsGenome <- grepl("genome", info$Description, ignore.case=TRUE)
		info$IsSequence <- grepl("sequence", info$Description, ignore.case=TRUE)
		info$IsChromosome <- grepl("chromosome", info$Name, ignore.case=TRUE)
		info$IsPlasmid <- grepl("plasmid", info$Name, ignore.case=TRUE)
		info$IsComplete <- grepl("complete", info$Description, ignore.case=TRUE)
	} else {
		info$Description <- ""
		info$IsGenome <- FALSE
		info$IsSequence <- FALSE
		info$IsChromosome <- FALSE
		info$IsPlasmid <- FALSE
		info$IsComplete <- FALSE
	} #if
	info$IsDraft <- grepl("draft", info$Definition, ignore.case=TRUE)
	info
} #function

#*function
PRINT.SeqDB <- function(db)
{
	cat("Sequence Database\nData directory: ", attr(db, "DataDir"))
	cat("\nCreation started: ", format(attr(db, "StartCreationDate")))
	cat("\nCreation finished: ", format(attr(db, "EndCreationDate")))
	cat("\nNumber of sequences: ", nrow(db))
	cat("\nTotal number of bases:  ", sum(as.numeric(db$Length)))
	cat("\n")
} #function

#*function
GetIndexFromName <- function(db, name, first.only=FALSE, ...)
{
	if (missing(db)) stop("Database not specified.")
	if (missing(name) || length(name)!=1 | !is.character(name) || nchar(name)==0)
		stop("An invalid name has been specified.")
idx <- which(grepl(name, db$Name, ignore.case=TRUE, ...))
if (length(idx)==0) warning(paste("No sequences found that match", name))
if (first.only) idx[1]
else idx
} #function

#*function
GetIndexFromAccession <- function(db, accession, exact=FALSE, first.only=FALSE, ...)
{
	if (missing(db)) stop("Database not specified.")
	if (missing(accession) || length(accession)!=1 | !is.character(accession))
		stop("An invalid accession number has been specified.")
	if (!exact)
	{
		m <- regexec("(.*)\\..*", accession)
		if (m[[1]][1]!=-1) 
			accession <- regmatches(accession, m)[[1]][2]
	} #if
idx <- which(grepl(accession, db$Accession, ignore.case=TRUE, ...))
if (length(idx)==0) warning(paste("No sequences found that match", accession))
if (first.only) idx[1]
else idx
} #function

#*function
LoadSequence <- function(db, key)
{
	if (length(key)==0 || length(key)>1) stop("Invalid key (name or sequence index).")
	if (is.character(key))
		idx <- GetIndex(db, key)
	else
	{
		if (!is.numeric(key)) stop("An invalid key has been specified.")
		idx <- key
	} #if
	if (length(idx)==0 || idx<1 || idx>nrow(db))
		stop("No sequence has been found or an invalid sequence index has been specified.")
	genome <- read.fasta(db$File[idx])
	if (length(genome)==1) genome[[1]]
	else genome
} #function

#*function
subset.SeqDB <- function(x, subset, select, drop=FALSE, ...)
{
	if (missing(x)) stop("No database has been specified.")
    if (missing(subset)) 
        r <- TRUE
    else {
        e <- substitute(subset)
        r <- eval(e, x, parent.frame())
        if (!is.logical(r)) 
            stop("'subset' must be logical")
        r <- r & !is.na(r)
    }
    if (missing(select)) 
        vars <- TRUE
    else {
        nl <- as.list(seq_along(x))
        names(nl) <- names(x)
        vars <- eval(substitute(select), nl, parent.frame())
    }
    structure(x[r, vars, drop = drop], 
	DataDir=attr(x, "DataDir"),
	StartCreationDate=Sys.time(), EndCreationDate=Sys.time())
} #function

#*function
GetGenomeIndices <- function(db)
{
	if (missing(db)) stop("No database has been specified.")
	if (!all.equal(class(db), c("SeqDB", "data.frame"))) stop("The specified database is not valid.")
	which(db$IsGenome)
} #function

#*function
GetGenomes <- function(db)
{
	if (missing(db)) stop("No database has been specified.")
	if (!all.equal(class(db), c("SeqDB", "data.frame"))) stop("The specified database is not valid.")
	subset(db, db$IsGenome)
} #function
