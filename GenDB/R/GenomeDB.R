#GenomeDB.R

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
	r <- regexec("^\\s*(\\S+\\s\\S+)\\s*", definition)
  m <- regmatches(definition, r)[[1]]
  info$Organism <- if (length(m)<2) definition else m[2]
	info
} #function

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
  if (class(x)[1]=="GenBankDB")
	  return(read.gbk(x$File[idx], ...))
  if (class(x)[1]=="SeqDB")
  {
	  data <- read.gbk(x$File[idx], ...)
	  if (length(data)==1) return(data[[1]])
	  else return(data)
	} #if
	stop("invalid database type")
} #function

