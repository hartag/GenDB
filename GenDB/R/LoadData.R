#LoadData.R


#Alternative generic function implementation
#LoadData <- function(x, ...)
#  UseMethod("LoadData")


LoadData.GenBankDB <- function(x, key, ...)
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
	InputFile <- file.path(attr(x, "DataDir"), x$File[idx])
	read.gbff(inputFile, ...)
} #function


LoadData.SeqDB <- function(x, key, ...)
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
	InputFile <- file.path(attr(x, "DataDir"), x$File[idx])
	data <- read.fasta(inputFile, ...)
	if (length(data)==1) data[[1]]
	else data
} #function
