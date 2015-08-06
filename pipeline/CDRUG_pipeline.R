########################
## Benjamin Haibe-Kains
## All rights Reserved
## February 11, 2013
########################



## This script runs the following scripts:
##  1. 'CDRUG_foo.R' defining all functions required for the rest of the analyses
##  2. 'normalization_cgp.R' performing curation, annotation and normalization of CGP data
##  3. 'normalization_ccle.R' performing curation, annotation and normalization of CCLE data
##  4. 'CDRUG_format.R' performing additional curation to identify common cell lines, tissue types and drugs investigated in CGP and CCLE
##  5. 'CDRUG_pipeline.R' computing all the correlations and generating all the tables/figures for the paper


## remove all existing objects from the workspace
rm(list=ls(all=TRUE))

require(parallel) || stop("Library parallel is not available")

## define functions required for the analysis pipeline
source(file.path("CDRUG_foo.R"))

########################
## global parameters

## set random seed to ensuer reproducibility of the resuls
set.seed(54321)

## number of cpu cores available for the analysis pipeline
## set to 'NULL' if all the available cores should be used
nbcore <- 32
availcore <- parallel::detectCores()
if (is.null(nbcore) || nbcore > availcore) { ncore <- availcore }
options("mc.cores"=nbcore)

## list of characters to be removed from row and column names
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

## directory where all the analysis results will be stored
saveres <- "saveres"
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE, recursive=TRUE) }

## minimum number of samples to compute correlation
minsample <- 10

## max fdr, threshold used to identify genes with significant concordance index
myfdr <- 0.20

## method to estimate the association gene-drug, controlled for tissue type
# genedrugm <- c("lm", "cindex")
genedrugm <- "lm"

## GSEA parameters
## path to the GSEA java executable
gsea.exec <- file.path("gsea2-2.0.13.jar")
## number of random geneset permutations
gsea.nperm <- 1000
## file containing the geneset definitions
genesets.filen <- file.path("c5.all.v4.0.entrez.gmt")
## minimum size for a geneset to be analyzed
min.geneset.size <- 15
## maximum size for a geneset to be analyzed
max.geneset.size <- 250

## number of most variant genes to consider for correlation
topvar <- 1000

if(!file.exists("CDRUG_log.txt")) {
  steps <- c("CDRUG_normalization_cgp", "CDRUG_normalization_ccle", "CDRUG_normalization_gsk", "CDRUG_format", "CDRUG_correlation_data", "CDRUG_analysis", "CDRUG_analysisbis", "CDRUG_analysis_gsk")
  progress.log <- cbind(steps, rep("...", length(steps)))
  dimnames(progress.log) <- list(paste("step", 1:length(steps), sep="."), c("script", "progress"))
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
} else {
    progress.log <- read.table(file=file.path("CDRUG_log.txt"), sep="\t", header=TRUE, stringsAsFactor=FALSE)
}

########################
## curation, annotation and normalization of CGP data

message("\n-----------------------------\n| Normalization of CGP data |\n-----------------------------")
if (progress.log["step.1", "progress"] != "done") {
  progress.log["step.1", "progress"] <- "in progress"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
  source("CDRUG_normalization_cgp.R")
  progress.log["step.1", "progress"] <- "done"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
}
message("\t-> DONE")

########################
## curation, annotation and normalization of CCLE data

message("\n------------------------------\n| Normalization of CCLE data |\n------------------------------")
if (progress.log["step.2", "progress"] != "done") {
  progress.log["step.2", "progress"] <- "in progress"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
  source("CDRUG_normalization_ccle.R")
  progress.log["step.2", "progress"] <- "done"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
}
message("\t-> DONE")

message("\n-----------------------------\n| Normalization of GSK data |\n-----------------------------")
if (progress.log["step.3", "progress"] != "done") {
  progress.log["step.3", "progress"] <- "in progress"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
  source("CDRUG_normalization_gsk.R")
  progress.log["step.3", "progress"] <- "done"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
}
message("\t-> DONE")

########################
## common cell lines, tissue types and drugs investigated in CGP and CCLE

message("\n-------------------------------------\n| Intersection between GGP and CCLE |\n-------------------------------------")
if (progress.log["step.4", "progress"] != "done") {
  progress.log["step.4", "progress"] <- "in progress"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
  source("CDRUG_format.R")
  progress.log["step.4", "progress"] <- "done"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
}
message("\t-> DONE")

########################
## script performing teh correlation analyses at the level of gene expressions

message("\n---------------------------------------------------------\n| Correlations of gene expressions between GGP and CCLE |\n---------------------------------------------------------")
if (progress.log["step.5", "progress"] != "done") {
  progress.log["step.5", "progress"] <- "in progress"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
  source("CDRUG_correlation_data.R")
  progress.log["step.5", "progress"] <- "done"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
}
message("\t-> DONE")

########################
## script generating all the figures and tables reported in the manuscript

message("\n-------------------------------------\n| Correlations between GGP and CCLE |\n-------------------------------------")
if (progress.log["step.6", "progress"] != "done") {
  progress.log["step.6", "progress"] <- "in progress"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
  source("CDRUG_analysis.R")
  progress.log["step.6", "progress"] <- "done"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
}
message("\t-> DONE")

########################
## script generating the figures reporting the resuts on tesing the most likely source of discrepancies

message("\n-----------------------------------------------\n| Look for the most likely source of discrepancies |\n-----------------------------------------------")
if (progress.log["step.7", "progress"] != "done") {
  progress.log["step.7", "progress"] <- "in progress"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
  source("CDRUG_analysisbis.R")
  progress.log["step.7", "progress"] <- "done"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
}
message("\t-> DONE")

########################
## script generating the figures for the GSK vs CGP vs CCLE comparison

message("\n-------------------------------------\n| Correlations between GSK and GGP/CCLE |\n-------------------------------------")
if (progress.log["step.8", "progress"] != "done") {
  progress.log["step.8", "progress"] <- "in progress"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
  source("CDRUG_analysis_gsk.R")
  progress.log["step.8", "progress"] <- "done"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("CDRUG_log.txt"), quote=FALSE)
}
message("\t-> DONE")

## save session info
write(toLatex(sessionInfo(), locale = FALSE), file="sessionInfoR.tex", append=FALSE)

