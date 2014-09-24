BuildSequenceDB <- function(dataDir, fasta.ext="*.fna", range=c(1,Inf), reg.exp=FALSE, ...)
{
#Check arguments
	if (length(range)==1) range <- c(1, range)

#Search for genome files
	cat("Searching for FASTA files ...\n")
	if (reg.exp) fasta.ext <- glob2rx(fasta.ext)
	files <- dir(dataDir, fasta.ext, recursive=TRUE, full.names=TRUE)
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
	  db$AmbiguousBases[i] <- db$Length[i]-sum(genome %in% c("a", "c", "g", "t"))
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

ParseFastaHeader <- function(header)
{
	if (is.null(header)) return(list())
	tokens <- regmatches(header, regexec("^.*\\|.*\\|(.*)\\|(.*)\\|(.*)$", header))[[1]]
	definition <- tokens[4]
	info <- list(Source=tokens[1], Definition=definition, Accession=sub("\\.\\d+$", "", tokens[3]), Version=tokens[3])
	info <- ParseGenBankDefinition(definition, info)
	info 
} #function

ParseGenBankDefinition <- function(definition, info=list())
{
	if (missing(definition) || !is.character(definition))
		stop("definition must be of type character")
	if (!is.list(info)) stop("info must be a list")
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

summary.SeqDB <- function(db, ...)
{
	x <- list()
	x$dataDir <- attr(db, "DataDir")
	x$startTime <- attr(db, "StartCreationDate")
	x$endTime <- attr(db, "EndCreationDate")
	x$nEntries <- nrow(db)
	x$registered <- sum(db$Registered)
	x$totalBases <- sum(as.numeric(db$Length))
	x$ambiguous <- sum(as.numeric(db$AmbiguousBases))
  bools <- subset(db, select=c("IsGenome", "IsSequence", "IsChromosome", "IsPlasmid", "IsComplete", "IsDraft"))
  tab <- sapply(bools, function(l) c("are"=sum(l), "are not"=sum(!l)))
  colnames(tab) <- substr(colnames(tab), 3, nchar(colnames(tab)))
  x$characteristics <- tab
	class(x) <- "summary.SeqDB"
	x
} #function

print.summary.SeqDB <- function(x, ...)
{
	cat("Sequence Database\nData directory: ", x$dataDir)
	cat("\nCreation started: ", format(x$startTime))
	cat("\nCreation finished: ", format(x$endTime))
	cat("\nNumber of database entries: ", x$nEntries)
	cat("\nNumber of successfully recorded sequences: ", x$registered)
	cat("\nTotal number of bases:  ", x$totalBases)
	cat("\nTotal number of ambiguous bases: ", x$ambiguous)
	cat("\n")
  if (!is.null(x$characteristics))
  print(x$characteristics)
  invisible(x)
} #function

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

LoadSequence <- function(db, key)
{
	if (length(key)==0 || length(key)>1) stop("Invalid key (name or sequence index).")
	if (is.character(key))
		idx <- GetIndexFromAccession(db, key)
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

GetGenomeIndices <- function(db, complete=TRUE)
{
	if (missing(db)) stop("no database has been specified")
	if (!all.equal(class(db), c("SeqDB", "data.frame"))) stop("the specified database is not valid")
	if (complete) which(db$IsGenome & db$IsComplete)
	else which(db$IsGenome)
} #function

GetGenomes <- function(db, complete=TRUE)
{
	if (missing(db)) stop("no database has been specified")
	if (!all.equal(class(db), c("SeqDB", "data.frame"))) stop("the specified database is not valid")
	if (complete) subset(db, db$IsGenome & db$IsComplete)
	else subset(db, db$IsGenome)
} #function
