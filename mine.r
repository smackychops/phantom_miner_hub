#!/usr/bin/env Rscript

library(RISmed)
library(ggplot2)
library(data.table)
library(stringr)
library(parallel)
library(stringr)
library(openssl)
library(rmarkdown)
library(parallel)

ggnice <- function(){
	theme(panel.grid.major = element_blank(), 
	      panel.grid.minor = element_blank(), 
	      panel.background = element_blank(), 
	      axis.line = element_line(colour = "black"), 
	      axis.text = element_text(color="black"))
}

ncore <- detectCores()

searchvector <- commandArgs(trailingOnly = TRUE)
search <- paste(searchvector,  collapse=" ")

# ####################################################################
# Get Data
getPubMed <- function(term,   limit=NULL)
{
	# Get number of papers
	res <- EUtilsSummary(term,  type='esearch',  db='pubmed')
	qc <- QueryCount(res)
	message(paste("Found ", qc, " papers in PubMed", sep=""))

	# Download them in chunks of 1000
	retmax <- 1000
	if(!is.null(limit))
	{
		if(qc>limit){qc <- limit}
	}
	hits <- list()
	if(qc<retmax){dovec<-c(0)}else{dovec<-c(0, seq(retmax, qc, by=retmax))}
	for(i in dovec)
	{
		message(paste("Downloading,  retstart=", i, sep=""))
		Sys.sleep(3)
		res2 <- EUtilsSummary(term,  type='esearch',  db='pubmed', retstart=i, retmax=retmax)
		get <- EUtilsGet(res2, type="efetch", db="pubmed")
		hits <- c(hits, get)
	}
	hits
}
hash <- md5(search)
rd <- file.path("rd", paste0(hash, ".rd"))
if(!file.exists(rd))
{
	message("Doing search and will save results to ", rd)
	pm <- getPubMed(search)
	save(pm, compress=TRUE, file=rd)
} else
{
	message("Found we did this search,  pulling saved data from ", rd)
	load(rd)
}

# This is a results list of tables to pre-build here,  then feed to the reporting script
res <- list(search=search)

title <- paste(searchvector, collapse="_")

res$title <- title

# ####################################################################
# Year Distribution
years <- as.data.frame(table(do.call(c, lapply(pm,  FUN=function(x) YearPubmed(x)))))
colnames(years) <- c("year", "papers")
years$year <- as.numeric(as.character(years$year))
years <- years[order(years$year, decreasing=T), ]

res$years <- years

# ####################################################################
# Affiliations ## not working - JJL
# Seems to only be returning one per article,  but MEDLINE should have one for each author
# Could do some text mining to pull out most frequent words to get countries,  universities,  etc
#aff <- as.data.frame(table(do.call(c, lapply(pm,  FUN=function(x) str_replace(Affiliation(x), ", .*", "")))))
#aff <- as.data.frame(table(do.call(c, lapply(pm,  FUN=function(x) Affiliation(x)))))
#aff <- aff[order(aff$Freq, decreasing=T), ]
#colnames(aff) <- c("place", "papers")
#res$aff <- aff

# ####################################################################
# Journals+Places
journal <- as.data.frame(table(do.call(c, lapply(pm,  FUN=function(x) MedlineTA(x)))))
journal <- journal[order(journal$Freq, decreasing=T), ]
colnames(journal) <- c("journal", "papers")
#head(journal, 30)
res$journal <- journal

getFreq <- function(pm, fun)
{
	tab <- as.data.frame(table(do.call(c, lapply(pm,  FUN=fun))))
	colnames(tab) <- c("name", "papers")
	tab <- tab[order(tab$papers, decreasing=T), ]
	tab
}

res$country <- getFreq(pm, Country)
res$grants <- getFreq(pm, GrantID)
#head(getFreq(pm, CollectiveName), 50)
head(getFreq(pm, Agency), 50)
head(getFreq(pm, Language), 50)
head(getFreq(pm, PublicationStatus), 50)
#head(getFreq(pm, PublicationType), 50)

# Citations
cites <- do.call(c, lapply(pm, Cited))
pmids <- do.call(c, lapply(pm, PMID))
jor <- do.call(c, lapply(pm, MedlineTA))
art <- do.call(c, lapply(pm, ArticleTitle))
yr <- do.call(c, lapply(pm, YearPubmed))
#typ <- do.call(c, lapply(pm, function(x) sapply(PublicationType(x), function(y) y[1])))

citation <- data.table(pmid=pmids, pmc_cites=cites, journal=jor, year=yr, title=art)
citation <- citation[order(citation$pmc_cites, decreasing=T), ]
res$citation <- citation
#getFreq(pm, PublicationType)

# ####################################################################
# MESH Terms

# Single Terms
getMesh <- function(x)
{
	mesh <- Mesh(x)
	mesh <- mesh[!is.na(mesh)]
	mesh <- do.call(c, lapply(mesh, function(x) as.character(x[x$Type=="Descriptor", ]$Heading)))
}
mesh <- do.call(c, lapply(pm,  FUN=getMesh))
mesh <- as.data.frame(table(mesh))
mesh <- mesh[order(mesh$Freq, decreasing=T), ]
res$mesh <- mesh

# Pairs
getMeshPairs <- function(x, n=2)
{
	m <- Mesh(x)
	m <- m[!is.na(m)]
	m <- lapply(m, function(x) as.character(x[x$Type=="Descriptor", ]$Heading))
	m <- lapply(m, function(x) x[!(x %in% c("Animals", "Humans", "Male", "Female"))])
	combos <- function(y)
	{
		if(length(y) > n)
		{
			allpairs <- t(combn(y, m=n))
			mc <- t(apply(allpairs, 1, function(c) c[order(c, decreasing=F)]))
		}
	}
	do.call(rbind, lapply(m, combos))
}
m2 <- do.call(rbind, mclapply(pm,  function(x) getMeshPairs(x, 2), mc.cores=ncore))
#m3 <- do.call(rbind, lapply(pm,  function(x) getMeshPairs(x, 3)))
#m4 <- do.call(rbind, lapply(pm,  function(x) getMeshPairs(x, 4)))
#m5 <- do.call(rbind, mclapply(pm,  function(x) getMeshPairs(x, 5), mc.cores=ncore))
#m6 <- do.call(rbind, mclapply(pm,  function(x) getMeshPairs(x, 6), mc.cores=ncore))
#m7 <- do.call(rbind, mclapply(pm,  function(x) getMeshPairs(x, 7), mc.cores=ncore))

countCombos <- function(x)
{
	x <- data.table(x)
	cols <- paste("V", 1:ncol(x), sep="")
	setkeyv(x, cols)
	x$pair <- ""
	x <- x[, list(papers=length(pair)), by=cols]
	x <- x[order(x$papers, decreasing=T), ]
}
#head(countCombos(m2), 25)
#head(countCombos(m3), 25)
#cc3 <- countCombos(m3)
#write.csv(cc3[1:3000, ], file="cc3.csv", row.names=F)
#head(countCombos(m4), 25)
#head(countCombos(m5), 25)
#head(countCombos(m6), 25)
#head(countCombos(m7), 25)

res$mesh_pairs <- countCombos(m2)

# Probably want to filter out stupid stuff here like humans,  animals,  male,  female
# Can also explore bringing in the qualifer terms

# ####################################################################
# Authors

# Name Frequency
Names <- function(x)
{
	peeps <- Author(x)
	peeps <- do.call(c, lapply(peeps,  FUN=function(x) paste(x$LastName, x$ForeName,  sep=",  ")))
	peeps
}
res$authors <- getFreq(pm, Names)

# First Author Frequency
FirstAuthor <- function(x)
{
	peeps <- Author(x)
	getFirst <- function(y)
	{
		if(!is.na(y[1, ]$LastName))
		{
			y <- y[y$order==1, ]
			paste(y$LastName, y$ForeName,  sep=",  ")
		}
	}
	peeps <- do.call(c, lapply(peeps,  FUN=getFirst))
	peeps
}
res$authors_first <- getFreq(pm, FirstAuthor)

# Last Author Frequency
LastAuthor <- function(x)
{
	peeps <- Author(x)
	getFirst <- function(y)
	{
		if(!is.na(y[1, ]$LastName))
		{
			y <- y[y$order==max(y$order), ]
			paste(y$LastName, y$ForeName,  sep=",  ")
		}
	}
	peeps <- do.call(c, lapply(peeps,  FUN=getFirst))
	peeps
}
res$authors_last <- getFreq(pm, LastAuthor)

# Author Pairs
getAuthorPairs <- function(x, n=2)
{
	peeps <- Author(x)
	peeps <- peeps[!is.na(peeps)]
	peeps <- lapply(peeps,  FUN=function(x) paste(x$LastName, x$ForeName,  sep=",  "))
	combos <- function(y)
	{
		if(length(y) > n)
		{
			allpairs <- t(combn(y, m=n))
			mc <- t(apply(allpairs, 1, function(c) c[order(c, decreasing=F)]))
		}
	}
	do.call(rbind, mclapply(peeps, combos, mc.cores=ncore))
}
a2 <- do.call(rbind, mclapply(pm,  function(x) getAuthorPairs(x, 2), mc.cores=ncore))
#head(countCombos(a2), 25)
res$author_pairs <- countCombos(a2)

a3 <- do.call(rbind, lapply(pm,  function(x) getAuthorPairs(x, 3)))
#head(countCombos(a3), 25)
res$author_trios <- countCombos(a3)

#a4 <- do.call(rbind, mclapply(pm,  function(x) getAuthorPairs(x, 4), mc.cores=ncore))
#head(countCombos(a4), 25)

#a5 <- do.call(rbind, lapply(pm,  function(x) getAuthorPairs(x, 5)))
#head(countCombos(a5), 25)

# Could be useful to go back and query each author in the co-authors with author of interest and drop in a list of their top MESH terms (if different from those with the author)
# Find out which topics the person is co-authoring on with who - make pairs of author pairs and MESH pairs (hah!)

# Plot graph
# Worry about this later
# Need some way to get edge weights into Gephi, etc
#library(igraph)
#ig <- graph.data.frame(pairs,  directed=FALSE)
#library(ggplot2)
#library(network)
#library(GGally)
#library(gridExtra)
#ng <- network(pairs,  directed=FALSE)
#pdf(file="pubmed-net.pdf", width=10.5,  height=8,  paper="special")
#ggnet(ng,  mode="fruchtermanreingold",  size=3,  segment.alpha=0.9,  alpha=0.2,  label.nodes=T,  color="black") + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),  "cm"))
#ggnet(ng,  mode="fruchtermanreingold",  size=12,  segment.alpha=0.9,  alpha=0.5,  label.nodes=F,  color="black",  weight.method="degree",  top8.nodes=T) + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),  "cm"))
#dev.off()

# Abstract text word cloud?
library(tm)
library(SnowballC)

abs <- do.call(c, lapply(pm, AbstractText))

corp <- Corpus(VectorSource(abs))
#corp <- tm_map(corp,  PlainTextDocument)
corp <- tm_map(corp,  removePunctuation)
corp <- tm_map(corp,  removeWords,  c("the", "this", stopwords('english')))
corp <- tm_map(corp,  stemDocument)

tdm <- TermDocumentMatrix(corp)

rs <- rowSums(as.matrix(tdm))
rs <- sort(rs, decreasing=T)
res$abswc <- data.table(term=names(rs), n=rs)

outpath <- file.path("rd", paste0(hash, "_res.rd"))
message("Saving results into ", outpath)
save(res, compress=TRUE, file=outpath)

rmarkdown::render('report.Rmd', output_file = paste(Sys.Date(),  res$title,  ".html",  sep="_"))
print(warnings())
print(res$title)
