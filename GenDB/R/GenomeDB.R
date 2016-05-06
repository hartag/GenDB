#GenomeDB.R

summary.GenomeDB <- function(db, ...)
{
	x <- list()
	x$dataDir <- attr(db, "DataDir")
	x$startTime <- attr(db, "StartCreationDate")
	x$endTime <- attr(db, "EndCreationDate")
	x$nEntries <- nrow(db)
	x$registered <- sum(db$Registered)
#	x$totalBases <- sum(as.numeric(db$Length))
#	x$ambiguous <- sum(as.numeric(db$AmbiguousBases))
  bools <- subset(db, select=c("IsGenome", "IsSequence", "IsChromosome", "IsPlasmid", "IsComplete", "IsDraft"))
  tab <- sapply(bools, function(l) c("are"=sum(l), "are not"=sum(!l)))
  colnames(tab) <- substr(colnames(tab), 3, nchar(colnames(tab)))
  x$characteristics <- tab
	class(x) <- "summary.GenomeDB"
	x
} #function

print.summary.GenomeDB <- function(x, ...)
{
	cat("Flat File Genome Database\nData directory: ", x$dataDir)
	cat("\nCreation started: ", format(x$startTime))
	cat("\nCreation finished: ", format(x$endTime))
	cat("\nNumber of database entries: ", x$nEntries)
	cat("\nNumber of successfully recorded sequences: ", x$registered)
#	cat("\nTotal number of bases:  ", x$totalBases)
#	cat("\nTotal number of ambiguous bases: ", x$ambiguous)
	cat("\n")
  if (!is.null(x$characteristics))
  print(x$characteristics)
  invisible(x)
} #function

GetIndexFromName <- function(db, name, exact=FALSE, first.only=FALSE, ...)
{
	if (missing(db)) stop("Database not specified.")
	if (!("GenomeDB" %in% class(db))) stop("the specified database is not valid")
	if (!is.character(name))
		stop("an invalid value for name has been given.")
	idx <- sapply(name, function(n) {
	  ind <- grep(n, db$Name, fixed=exact, ignore.case=!exact, ...)
	  if (length(ind)==0)
	  {
	  	ind <- as.character(NA)
	  	names(ind) <- n
	  }
	  else
	  	names(ind) <- db$Name[ind]
	  ind
	}, simplify=FALSE, USE.NAMES=FALSE)
	if (first.only) idx <- lapply(idx, "[", 1)
	unlist(idx)
} #function

GetIndexFromAccession <- function(db, accession, exact=FALSE, first.only=FALSE, ...)
{
	if (missing(db)) stop("Database not specified.")
	if (!("GenomeDB" %in% class(db))) stop("the specified database is not valid")
	if (!is.character(accession))
		stop("an invalid value for accession has been given.")
		accessionField <- if (class(db)[1]=="GenBankDB") "Version" else "Accession"
	idx <- sapply(accession, function(acc) {
	  ind <- grep(acc, db[[accessionField]], fixed=exact, ignore.case=!exact, ...)
	  if (length(ind)==0)
	  {
	  	ind <- as.character(NA)
	  	names(ind) <- acc
	  }
	  else
	  	names(ind) <- db$Accession[ind]
	  ind
	}, simplify=FALSE, USE.NAMES=FALSE)
	if (first.only) idx <- lapply(idx, "[", 1)
	unlist(idx)
} #function

subset.GenomeDB <- function(x, subset, select, drop=FALSE, ...)
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
	if (!("GenomeDB" %in% class(db))) stop("the specified database is not valid")
	if (complete) which(db$IsGenome & db$IsComplete)
	else which(db$IsGenome)
} #function

GetGenomes <- function(db, complete=TRUE)
{
	if (missing(db)) stop("no database has been specified")
	if (!("GenomeDB" %in% class(db))) stop("the specified database is not valid")
	if (complete) subset(db, db$IsGenome & db$IsComplete)
	else subset(db, db$IsGenome)
} #function


LoadData <- function(x, key, ...)
{
	if (length(key)==0 || length(key)>1) stop("Invalid key (name or sequence index).")
	if (is.character(key))
		idx <- GetIndexFromAccession(x, key)
	else
	{
		if (!is.numeric(key)) stop("An invalid key has been specified.")
		idx <- key
	} #if
	if (length(idx)==0 || idx<1 || idx>nrow(x))
		stop("No sequence has been found or an invalid sequence index has been specified.")
	inputFile <- file.path(attr(x, "DataDir"), x$File[idx])
  if (class(x)[1]=="GenBankDB")
	  return(read.gbk(inputFile, ...))
  if (class(x)[1]=="SeqDB")
  {
	  data <- read.fasta(inputFile, ...)
	  if (length(data)==1) return(data[[1]])
	  else return(data)
	} #if
	stop("invalid database type")
} #function

Location <- function(db)
{
	if (missing(db)) stop("no database has been specified")
	if (!("GenomeDB" %in% class(db))) stop("the specified database is not valid")
  attr(db, "DataDir")
} #function

`Location<-` <- function(db, value)
{
	if (missing(db)) stop("no database has been specified")
	if (!("GenomeDB" %in% class(db))) stop("the specified database is not valid")
	if (missing(value)) stop("a new location has not been specified")
  attr(db, "DataDir") <- value
  invisible(db) #return database with new location
} #function

