#read.gbff.R

read.gbff <- function(file, entries=FALSE, regions=c("gene", "CDS", "misc_feature"), sequence=FALSE)
{
#Check arguments
  regions <- match.arg(regions)
  
#Load GenBank file
  if (entries)
    chunk <- readLines(file, warn=FALSE) #read file
  else
    chunk <- readLines(file, n=200, warn=FALSE) #read file
  lineCount <- length(chunk)

#Get genome information
##Get and parse Definition line
  header <- regexec("DEFINITION\\s+(.+)$", chunk[2])
  m <- regmatches(chunk[2], header)[[1]]
  definition <- m[[2]]
  if (grepl("^\\s+\\S.*$", chunk[3])[[1]])
  {
    header <- regexec("^\\s+(\\S.*)$", chunk[3])
    m <- regmatches(chunk[3], header)[[1]]
    definition <- paste(definition, m[[2]])
  } #if
  info <- list(Definition=definition)
  info <- ParseGenBankDefinition(definition, info)
  
##Extract information from first line
  header <- regexec("LOCUS\\s+([^ \t]+)\\s+(\\d+)\\s+bp\\s+([^ \t]+)\\s+([^ \t]+)\\s+(...)\\s+(.+)$", chunk[1])
  m <- regmatches(chunk[1], header)[[1]]
  info$Accession <- m[[2]]
  info$Length <- as.integer(m[[3]])
  info$Content <- m[[4]]
  info$PhysicalStructure <- m[[5]]
  info$GB <- m[[6]]
  info$LastUpdate <- m[[7]]
  info$File <- file
  verLine <- which(grepl("^VERSION", chunk[1:8]))[1]
  if (!is.null(verLine))
  {
    header <- regexec("VERSION\\s+(.+)\\s+(\\S+).*$", chunk[verLine])
    m <- regmatches(chunk[verLine], header)[[1]]
    info$Version <- m[[2]]
    info$GI <- m[[3]]
  } #if

##Get taxanomic classification information
  line <- which(grepl("^SOURCE", chunk[1:20]))[1]
  taxInfo <- NULL
  taxRanks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  if (!is.null(line) && is.numeric(line)) {
    if (grepl("ORGANISM", chunk[line+1L])) {
      taxLines <- paste0(chunk[line+2L], chunk[line+3L])
      taxLines <- gsub("^[^A-Za-z]+|[^A-Za-z]+$", "", taxLines)
      taxInfo <- strsplit(taxLines, ";\\s*")[[1]]
      names(taxInfo) <- taxRanks[seq_along(taxInfo)]
    } #if
  } #if
  info$Domain <- NA
  info$Phylum <- NA
  info$Class <- NA
  info$Order <- NA
  info$Family<- NA
  info$Genus <- NA
  if (!is.null(taxInfo)) {
    info$Domain <- taxInfo[1]
    info$Phylum <- taxInfo[2]
    info$Class <- taxInfo[3]
    info$Order <- taxInfo[4]
    info$Family<- taxInfo[5]
    info$Genus <- taxInfo[6]
  } #if
##Get species
#  info$Species <- strsplit(info$BinomialName, "\\s+")[[1]]
#  if (length(info$Species)>=2) info$Species <- info$Species[2]
  info$Species <- info$BinomialName
##Get strain and substrain
  info$Strain <- NA
  line <- which(grepl("/strain=\"", chunk))[1]
  if (!is.null(line) && is.numeric(line)) {
    info$Strain <- gsub("^\\s*/strain=\"|\"[^\"]*$", "", chunk[line])
  } #if
  info$Substrain <- NA
  line <- which(grepl("/sub_strain=\"", chunk))[1]
  if (!is.null(line) && is.numeric(line)) {
    info$Substrain <- gsub("^\\s*/sub_strain=\"|\"[^\"]*$", "", chunk[line])
  } #if
##Get taxanomic ID
  info$taxid <- NA
  line <- which(grepl("/db_xref=\"taxon:", chunk))[1]
  if (!is.null(line) && is.numeric(line)) {
    info$taxid <- gsub("^\\s*/db_xref=\"taxon:|\"[^\"]*$", "", chunk[line])
  } #if
    
#Read in the genome sequence
  if (sequence)
  {
  seqStart <- grep("^ORIGIN", chunk)
  seqEnd <- grep("^//", chunk)
  if (length(seqStart)>1 || length(seqEnd)>1)
    warning("gbk file contains more than one ORIGIN.  Only the first one has been read.")
  if (length(seqStart)==0 || length(seqEnd)==0)
    stop("gbk file contains no ORIGIN.")
  seq <- unlist(strsplit(chunk[(seqStart[1]+1):(seqEnd[1]-1)], ""))
  seq <- seq[seq %in% letters]
  seq.rc <- rev(comp(seq))
  } #if

#Get extents of regions
  if (entries)
  {
  m <- regexec(paste0(regions, "\\s+(\\d+)\\.\\.(\\d+)"), chunk)
  cds <- regmatches(chunk, m)
  cds[sapply(cds,length)==0] <- NULL
  start <- as.numeric(sapply(cds, function(x) x[[2]]))
  stop <- as.numeric(sapply(cds, function(x) x[[3]]))
  genes1<- data.frame(start=start, stop=stop)
  m <- regexec(paste0(regions, "\\s+complement\\((\\d+)\\.\\.(\\d+)\\)"), chunk)
  cds <- regmatches(chunk, m)
  cds[sapply(cds,length)==0] <- NULL
  stop <- -as.integer(sapply(cds, function(x) x[[2]]))
  start <- -as.integer(sapply(cds, function(x) x[[3]]))
  genes2 <- data.frame(start=start, stop=stop)
genes2 <- info$Length+1L+genes2
genes2 <- genes2[order(genes2$start),]
if (sequence)
{
  genes1$start.codon <- sapply(genes1$start, function(i) paste(seq[i:(i+2)], collapse=""))
  genes1$stop.codon <- sapply(genes1$stop, function(i) paste(seq[(i-2):i], collapse=""))
  genes2$start.codon <- sapply(genes2$start, function(i) paste(seq.rc[i:(i+2)], collapse=""))
  genes2$stop.codon <- sapply(genes2$stop, function(i) paste(seq.rc[(i-2):i], collapse=""))
} #if
gene.list <- rbind(genes1, genes2)
gene.list <- data.frame(gene.list, strand=c(rep(1, dim(genes1)[1]), rep(2, dim(genes2)[1])), rf=((gene.list$start-1) %% 3)+1)
  } #if

#Finish preparing return
if (entries) info$regions<- gene.list
if (sequence)
{
  info$primary <- seq
  info$complementary <- seq.rc
} #if
#info$taxonomy <- taxInfo
class(info) <- "gbk.info"
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
	r <- regexec("^\\s*(\\S+\\s\\S+)\\s*", definition)
  m <- regmatches(definition, r)[[1]]
  info$BinomialName <- if (length(m)<2) definition else m[2]
	info
} #function
