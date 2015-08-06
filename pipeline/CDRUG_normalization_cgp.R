########################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## February 11, 2013
########################

# rm(list=ls())

datapath <- file.path("data", "CGP")
rawpath <- file.path(datapath, "raw_ge")
if(!file.exists(rawpath)) { dir.create(rawpath, showWarnings=FALSE, recursive=TRUE) }
  
require(gdata) || stop("Library gdata is not available!")
require(R.utils) || stop("Library R.utils is not available!")


########################
## download data
########################
ftpdir <- "ftp://ftp.ebi.ac.uk//pub/databases/microarray/data/experiment/MTAB/E-MTAB-783/"
myfn <- file.path(rawpath, "celfile_timestamp.RData")
if(!file.exists(myfn)) {
  message("Download genomic data")
  
  require(R.utils) || stop("Library R.utils is not available!")
  
  dir.create(file.path(rawpath, "dwl"), showWarnings=FALSE)
  
  ## download and compress CEL files
  celfile.timestamp <- celfn <- NULL
  i <- 1
  while(i <= 9) {
    ## assuming there are only 9 zip archives (need to check if the update version has more)
   dwl.status <- download.file(url=sprintf("%s/E-MTAB-783.raw.%i.zip", ftpdir, i), destfile=file.path(rawpath, "dwl", sprintf("E-MTAB-783.raw.%i.zip", i)))
   if(dwl.status != 0) {
     message("\t-> download failed, let's try again ...")
     file.remove(file.path(rawpath, "dwl", sprintf("E-MTAB-783.raw.%i.zip", i)))
     i <- i - 1
    } else {
       ## unzip archive
       fff <- unzip(zipfile=file.path(rawpath, "dwl", sprintf("E-MTAB-783.raw.%i.zip", i)), list=TRUE)
       celfile.timestamp <- c(celfile.timestamp, as.character(fff[ ,"Date"]))
       celfn <- c(celfn, as.character(fff[ ,"Name"]))
       res <- unzip(zipfile=file.path(rawpath, "dwl", sprintf("E-MTAB-783.raw.%i.zip", i)), exdir=rawpath)
       ## compress each CEL file individually using gzip
       library(R.utils)
       sapply(file.path(rawpath, as.character(fff[ ,"Name"])), R.utils::gzip, overwrite=TRUE)
       i <- i + 1
     }
  }
  celfile.timestamp <- t(sapply(strsplit(celfile.timestamp, split=" "), function(x) { return(x) }))
  dimnames(celfile.timestamp) <- list(celfn, c("file.day", "file.hour"))
   
  # unlink(file.path(rawpath, "dwl"), recursive=TRUE)
  write.csv(celfile.timestamp, file=file.path(rawpath, "celfile_timestamp.csv"))
  save(list=c("celfile.timestamp"), compress=TRUE, file=myfn)
}

## download sample information
message("Download sample information")
dwl.status <- download.file(url=sprintf("%s/E-MTAB-783.sdrf.txt", ftpdir), destfile=file.path(rawpath, "dwl", "E-MTAB-783.sdrf.txt"))
if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
file.copy(from=file.path(rawpath, "dwl", "E-MTAB-783.sdrf.txt"), to=file.path(rawpath, "E-MTAB-783.sdrf.txt"))
  
## download drug sensitivity (release 2)
message("Download drug sensitivity measurements")
dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/releases/release-2.0/gdsc_manova_input_w2.csv", destfile=file.path(rawpath, "dwl", "gdsc_manova_input_w2.csv"))
if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
file.copy(from=file.path(rawpath, "dwl", "gdsc_manova_input_w2.csv"), to=file.path(rawpath, "cgp_drug_sensitivity.csv"))

## download drug concentration (release 2)
message("Download screening drug concentrations")
dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/current_release/gdsc_compounds_conc_w2.csv", destfile=file.path(rawpath, "dwl", "gdsc_compounds_conc_w2.csv"))
if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
file.copy(from=file.path(rawpath, "dwl", "gdsc_compounds_conc_w2.csv"), to=file.path(rawpath, "cgp_drug_concentration.csv"))
  
## download cell line annotations and COSMIC IDs
## annotations from COSMIC cell line project
dda <- "v65_280513"
dwl.status <- download.file(url=sprintf("ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicCellLineProject_%s.tsv.gz", dda), destfile=file.path(rawpath, "dwl", sprintf("CosmicCellLineProject_%s.tsv.gz", dda)))
if(dwl.status != 0) { stop("Download failed, please rerun the pipeline! It may be that there is a new version of the file CosmicCellLineProject, please look at ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/ and update the script accordingly ...") }
## untar
res <- R.utils::gunzip(filename=file.path(rawpath, "dwl", sprintf("CosmicCellLineProject_%s.tsv.gz", dda)), overwrite=TRUE)
file.copy(from=file.path(rawpath, "dwl", sprintf("CosmicCellLineProject_%s.tsv", dda)), to=file.path(rawpath, "cosmic_celline_collection.csv"))
cosmic.celline <- read.csv(file=file.path(rawpath, "cosmic_celline_collection.csv"), sep="\t", stringsAsFactors=FALSE)
cosmic.celline[cosmic.celline == "" | cosmic.celline == " " | cosmic.celline == "  "] <- NA
cosmic.celline <- cosmic.celline[- c(grep("row selected", cosmic.celline[ ,1]), grep("rows selected", cosmic.celline[ ,1])), , drop=FALSE]
## remove cell line with no name
cosmic.celline <- cosmic.celline[!is.na(cosmic.celline[ , "Sample.name"]), , drop=FALSE]
## merge the gene targets
dupln <- unique(cosmic.celline[ , "Sample.name"][duplicated(cosmic.celline[ , "Sample.name"])])
tt <- cosmic.celline
iix.rm <- NULL
for(i in 1:length(dupln)) {
  duplix <- cosmic.celline[ ,"Sample.name"] == dupln[i]
  iix <- sort((which(duplix)), decreasing=FALSE)[1]
  iix.rm <- c(iix.rm, setdiff(which(duplix), iix))
  tt[iix, "Gene.name"] <- paste(cosmic.celline[duplix, "Gene.name"], collapse="///")
  tt[iix, "UniProt.ID"] <- paste(cosmic.celline[duplix, "UniProt.ID"], collapse="///")
  tt[iix, "Zygosity"] <- paste(cosmic.celline[duplix, "Zygosity"], collapse="///")
  tt[iix, "CDS_MUT_SYNTAX"] <- paste(cosmic.celline[duplix, "CDS_MUT_SYNTAX"], collapse="///")
  tt[iix, "AA_MUT_SYNTAX"] <- paste(cosmic.celline[duplix, "AA_MUT_SYNTAX"], collapse="///")
  tt[iix, "NCBI36.genome.position"] <- paste(cosmic.celline[duplix, "NCBI36.genome.position"], collapse="///")
  tt[iix, "GRCh37.genome.position"] <- paste(cosmic.celline[duplix, "GRCh37.genome.position"], collapse="///")
}
tt <- tt[-iix.rm, , drop=FALSE]
rownames(tt) <- tt[ , "Sample.name"]
cosmic.celline <- tt
## annotations from GDSC (Genomics of Drug Sensitivity in Cancer)
dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/current_release/gdsc_cell_lines_w2.csv", destfile=file.path(rawpath, "dwl", "gdsc_cell_lines_w2.csv"))
file.copy(from=file.path(rawpath, "dwl", "gdsc_cell_lines_w2.csv"), to=file.path(rawpath, "cgp_celline_collection.csv"))
if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
gdsc.celline <- read.csv(file=file.path(rawpath, "cgp_celline_collection.csv"), stringsAsFactors=FALSE)
gdsc.celline[gdsc.celline == "" | gdsc.celline == " " | gdsc.celline == "  "] <- NA
gdsc.celline <- gdsc.celline[!is.na(gdsc.celline[ , "CELL_LINE_NAME"]), , drop=FALSE]
dupln <- unique(gdsc.celline[ , "CELL_LINE_NAME"][duplicated(gdsc.celline[ , "CELL_LINE_NAME"])])
gdsc.celline <- gdsc.celline[!duplicated(gdsc.celline[ , "CELL_LINE_NAME"]), , drop=FALSE]
rownames(gdsc.celline) <- gdsc.celline[ , "CELL_LINE_NAME"]
## merge GDSC and COSMIC annotations through COSMIC_ID
iix <- which(!is.na(gdsc.celline[ , "COSMIC_ID"]) & !is.element(gdsc.celline[ , "COSMIC_ID"], cosmic.celline[ , "ID_sample"]))
tt <- data.frame(matrix(NA, nrow=nrow(cosmic.celline) + length(iix), ncol=ncol(cosmic.celline), dimnames=list(c(rownames(cosmic.celline), rownames(gdsc.celline)[iix]), colnames(cosmic.celline))))
tt[rownames(cosmic.celline), ] <- cosmic.celline
tt[rownames(gdsc.celline)[iix], "Sample.name"] <- gdsc.celline[iix, "CELL_LINE_NAME"]
tt[rownames(gdsc.celline)[iix], "ID_sample"] <- gdsc.celline[iix, "COSMIC_ID"]
celline.cgp <- tt

## download drug information
message("Download drug information")
dwl.status <- download.file(url="http://www.nature.com/nature/journal/v483/n7391/extref/nature11005-s2.zip", destfile=file.path(rawpath, "dwl", "nature11005-s2.zip"))
if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
ff <- as.character(unzip(zipfile=file.path(rawpath, "dwl", "nature11005-s2.zip"), list=TRUE)[1, 1])
unzip(zipfile=file.path(rawpath, "dwl", "nature11005-s2.zip"), exdir=file.path(rawpath, "dwl"))
file.copy(from=file.path(rawpath, "dwl", ff), to=file.path(rawpath, "nature_supplementary_information.xls"))

########################
## normalize and format data
########################
myfn <- file.path(saveres, "cgp_frma.RData")
if(!file.exists(myfn)) {

  require(affy) || stop("Library affy is not available!")
  require(Hmisc) || stop("Library Hmisc is not available!")
  require(genefu) || stop("Library genefu is not available!")
  require(frma) || stop("Library frma is not available!")
  require(hthgu133afrmavecs) || stop("Library hthgu133afrmavecs is not available!")
  data(hthgu133afrmavecs)
  require(hthgu133acdf) || stop("Library hthgu133acdf is not available!")
  data(hthgu133acdf)

  ## CEL file names
  celfn <- list.celfiles(rawpath, full.names=TRUE)
  celfns <- list.celfiles(rawpath, full.names=FALSE)
  ## experiments' names
  names(celfn) <- names(celfns) <- gsub(".CEL.gz", "", celfns)
  ## chip type and date
  chipt <- sapply(celfn, celfileChip)
  chipd <- t(sapply(celfn, celfileDateHour))
  ## reorder CEL files by hybridization time or timestamp
  myx <- NULL
  if(any(!complete.cases(chipd))) {
    ## all hybridization dates are not available
    load(file.path(rawpath, "celfile_timestamp.RData"))
    if(!all(is.element(celfns, paste(rownames(celfile.timestamp), "gz", sep=".")))) { stop("Timestamp is not available for all CEL files!") }
      celfile.timestamp <- celfile.timestamp[match(celfns, paste(rownames(celfile.timestamp), "gz", sep=".")), , drop=FALSE] 
      myx <- order(celfile.timestamp[ ,"file.day"], celfile.timestamp[ ,"file.hour"], decreasing=FALSE)
  } else {
      myx <- order(chipd[ ,"day"], chipd[ ,"hour"], decreasing=FALSE)
  }
  celfn <- celfn[myx]
  celfns <- celfns[myx]
  chipt <- chipt[myx]
  chipd <- chipd[myx, , drop=FALSE]
  celfile.timestamp <- celfile.timestamp[myx, , drop=FALSE]


  ## read info about drugs and experiments

  ## phenotype for the drugs
  message("Read drug sensitivity measurements")
  myfn2 <- file.path(saveres, "cgp_drug_sensitivity.RData")
  if(!file.exists(myfn2)) {
    drugpheno <- read.csv(file.path(rawpath, "cgp_drug_sensitivity.csv"), stringsAsFactors=FALSE)
    drugpheno[drugpheno == "" | drugpheno == " "] <- NA
    save(list="drugpheno", compress=TRUE, file=myfn2)
  } else { load(myfn2) }
  ## format column names
  coln2 <- unlist(drugpheno[1, ,drop=TRUE])
  coln2[coln2 == ""] <- NA
  drugpheno <- drugpheno[-1, ,drop=FALSE]
  coln <- colnames(drugpheno)
  coln2[is.na(coln2)] <- coln[is.na(coln2)]
  coln2 <- genefu::rename.duplicate(x=coln2, sep="_dupl")$new.x
  myx <- sapply(sapply(strsplit(coln2, "_"), function(x) { return(x[[1]]) }), Hmisc::all.is.numeric)
  coln2[myx] <- paste("drugid", gsub(pattern=badchars, replacement="_", x=toupper(coln2[myx])), sep="_")
  colnames(drugpheno) <- coln2
  ## drug identifiers and names
  dn <- toupper(gsub(badchars, "", sapply(strsplit(coln, "_"), function(x) { return(x[[1]]) })))
  ## manual curation for drug names starting with a figure
  dn[!is.na(dn) & dn == "X17AAG"] <- "17AAG"
  dn[!is.na(dn) & dn == "X681640"] <- "681640"
  did <- sapply(strsplit(coln2, "_"), function(x) { if(x[[1]] == "drugid") { return(x[[2]]) } else { return(NA) } })
  drugnid <- cbind("drug.name"=dn, "drug.id"=did)[!is.na(did) & !duplicated(did), ]
  rownames(drugnid) <- paste("drugid", drugnid[ , "drug.id"], sep="_")

  ## cell line identifiers
  dupln <- duplicated(drugpheno[ ,"Cell.Line"])
  if(sum(dupln) > 1) { warning("some cell lines are duplicated, only the first instance is kept") }
  drugpheno <- drugpheno[!dupln, , drop=FALSE]
  if(any(!is.element(drugpheno[ ,"Cell.Line"], celline.cgp[ , "Sample.name"]))) { stop("Some cell line names are not included in the COSMIC database") }
  celln <- drugpheno[ ,"Cell.Line"]
  drugpheno <- data.frame("cellid"=celln, drugpheno)
  rownames(drugpheno) <- celln
  
  ## protein coding variants
  ## Genetic mutation data for cancer genes. Includes MSI status (1 = unstable and 0 = stable) and gene-fusions. A binary code 'x::y' description is used for each gene where 'x' identifies a coding variant and 'y' indicates copy number information from SNP6.0 data. For gene fusions, cell lines are identified as fusion not-detected (0) or the identified fusion is given. The following abbreviations are used: not analysed (na), not detected or wild-type (wt), no copy number information (nci).
  ## we assume that AKT2 and WT1 are the first and last genes in the file
  rangeg <- which(colnames(drugpheno) == "AKT2"):which(colnames(drugpheno) == "WT1")
  mutation <- as.matrix(drugpheno[ , rangeg, drop=FALSE])
  mutation <- apply(X=mutation, MARGIN=c(1, 2), FUN=function(x) {
    x <- unlist(strsplit(x, split="::"))
    if(length(x) == 2) {
      if(!is.na(x[[1]]) && (x[[1]] == "na" || x[[1]] == "p.?" || x[[1]] == "p.0?")) {
        x <- NA
      } else {
        x <- x[[1]]
      }
    } else { x <- NA }
    return(x)
  })

  ## info about each experiment
  message("Read sample information")
  sampleinfo <- read.csv(file.path(rawpath, "E-MTAB-783.sdrf.txt"), sep="\t", stringsAsFactors=FALSE)
  sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
  ## curate cell line names
  sampleinfo[sampleinfo[ , "Source.Name"] == "MZ2-MEL.", "Source.Name"] <- "MZ2-MEL"
  iix <- which(!duplicated(sampleinfo[ , "Source.Name"]) & !is.element(sampleinfo[ , "Source.Name"], celline.cgp[ , "Sample.name"]))
  if(length(iix) > 0) {
    ## enrich the list of cell lines
    tt <- matrix(NA, nrow=length(iix), ncol=ncol(celline.cgp), dimnames=list(sampleinfo[iix, "Source.Name"], colnames(celline.cgp)))
    tt[ , "Sample.name"] <- sampleinfo[iix, "Source.Name"]
    celline.cgp <- rbind(celline.cgp, tt)
  }
  fn <- gsub(patter="[.]CEL", replacement="", x=sampleinfo[ ,"Array.Data.File"])
  if(any(!is.element(fn[!is.na(fn)], names(celfns)))) { stop("some CEL files are missing for the CGP project") }
  rownames(sampleinfo) <- fn
  sampleinfo <- sampleinfo[names(celfn), , drop=FALSE]
  sampleinfo <- data.frame("samplename"=names(celfns), "filename"=celfns, "chiptype"=chipt, "hybridization.date"=chipd[ ,"day"], "hybridization.hour"=chipd[ ,"hour"], "file.day"=celfile.timestamp[ ,"file.day"], "file.hour"=celfile.timestamp[ ,"file.hour"], "batch"=NA, "cellid"=sampleinfo[ , "Source.Name"], sampleinfo)
  sampleinfo2 <- sampleinfo
  ## remove duplcated cell line hybridization
  sampleinfo <- sampleinfo[!duplicated(sampleinfo[ ,"cellid"]), , drop=FALSE]
  rownames(sampleinfo) <- sampleinfo[ ,"cellid"]

  ## update of cgp cell line collection
  celline.cgp <- data.frame("cellid"=as.character(celline.cgp[ , "Sample.name"]), celline.cgp)
  celline.cgp[ , "cellid"] <- as.character(celline.cgp[ , "cellid"])
  ## add url based on COSMIC IDs
  uurl <- paste("http://cancer.sanger.ac.uk/cell_lines/sample/overview?id=", celline.cgp[ , "ID_sample"], sep="")
  uurl[is.na(celline.cgp[ , "ID_sample"])] <- NA
  celline.cgp <- data.frame("cellid"=celline.cgp[ , "cellid"], "link"=uurl, celline.cgp[ , !is.element(colnames(celline.cgp), "cellid")])

  ## drugpheno
  cellnall <- sort(unique(c(as.character(sampleinfo[ ,"cellid"]), as.character(drugpheno[ ,"cellid"]))))
  dd <- data.frame(matrix(NA, ncol=ncol(drugpheno), nrow=length(cellnall), dimnames=list(cellnall, colnames(drugpheno))))
  newlev <- sapply(drugpheno, levels)
  newlev$cellid <- cellnall
  dd <- setcolclass.df(df=dd, colclass=sapply(drugpheno, class), factor.levels=newlev)
  dd[rownames(drugpheno),colnames(drugpheno)] <- drugpheno
  dd[ ,"cellid"] <- cellnall
  drugpheno <- dd
  
  ## mutation
  dd <- matrix(NA, ncol=ncol(mutation), nrow=length(cellnall), dimnames=list(cellnall, colnames(mutation)))
  dd[rownames(mutation), colnames(mutation)] <- mutation
  rownames(dd) <- cellnall
  mutation <- dd

  ## reproducibility between different screening sites
  ## camptothecin was screened at MGH (drug id 195) and WTSI (drug id 1003)
  ## data only available in the supplementary infomration of the Nature website
  myfn2 <- file.path(saveres, "nature_supplinfo_drugpheno_cgp.RData")
  if(!file.exists(myfn2)) {
    drugpheno.nature <- gdata::read.xls(xls=file.path(rawpath, "nature_supplementary_information.xls"), sheet=2)
    drugpheno.nature[drugpheno.nature == "" | drugpheno.nature == " "] <- NA
    save(list="drugpheno.nature", compress=TRUE, file=myfn2)
  } else { load(myfn2) }
  ## format column names
  coln2 <- gsub(" ", "", sapply(drugpheno.nature[1,], as.character))
  coln2[coln2 == ""] <- NA
  drugpheno.nature <- drugpheno.nature[-1, ,drop=FALSE]
  coln <- colnames(drugpheno.nature)
  coln2[is.na(coln2)] <- coln[is.na(coln2)]
  coln2 <- genefu::rename.duplicate(x=coln2, sep="_dupl")$new.x
  myx <- sapply(sapply(strsplit(coln2, "_"), function(x) { return(x[[1]]) }), Hmisc::all.is.numeric)
  coln2[myx] <- paste("drugid", gsub(pattern=badchars, replacement="_", x=toupper(coln2[myx])), sep="_")
  colnames(drugpheno.nature) <- coln2
  myx <- sapply(strsplit(colnames(drugpheno.nature), "_"), function(x) { return(all(x[c(length(x)-1, length(x))] == c("IC", "50"))) })
  ic50 <- drugpheno.nature[ , myx, drop=FALSE]
  nn <- dimnames(ic50)
  nn[[2]] <- gsub("_IC_50", "", nn[[2]])
  ic50 <- apply(ic50, 2, as.numeric)
  dimnames(ic50) <- nn
  ic50 <- exp(ic50) / 10^6
  ## camptothecin
  pdf(file.path(saveres, "cgp_camptothecin_mgh_wtsi_paper.pdf"))
  yy <- -log10(ic50[ , "drugid_195", drop=FALSE])
  xx <- -log10(ic50[ , "drugid_1003", drop=FALSE])
  ccix <- complete.cases(xx, yy)
  nnn <- sum(ccix)
  cc <- cor.test(x=xx, y=yy, method="spearman", use="complete.obs")
  cci <- spearmanCI(x=cc$estimate, n=sum(ccix))
  par(mar=c(4, 4, 3, 1) + 0.1)
  llim <- round(range(c(xx, yy), na.rm=TRUE) * 10) / 10
  myScatterPlot(x=xx, y=yy, xlab="-log10 IC50 (WTSI)", ylab="-log10 IC50 (MGH)", main="CAMPTOTHECIN", pch=16, method="transparent", transparency=0.75)
  legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("Rs=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, nnn), text.font=2)
  dev.off()

  ## drug information
  message("Read drug information")
  myfn2 <- file.path(saveres, "nature_supplinfo_druginfo_cgp.RData")
  if(!file.exists(myfn2)) {
    druginfo <- gdata::read.xls(xls=file.path(rawpath, "nature_supplementary_information.xls"), sheet=4)
    druginfo[druginfo == "" | druginfo == " "] <- NA
    save(list="druginfo", compress=TRUE, file=myfn2)
  } else { load(myfn2) }
  druginfo <- data.frame("drugid"=gsub(pattern =badchars, replacement="", x=toupper(druginfo[ ,"Drug.ID"])), druginfo)
  rownames(druginfo) <- paste("drugid", as.character(druginfo[ ,"drugid"]), sep="_")

  ## drug concentration
  message("Read drug concentration")
  drugconc <- read.csv(file.path(rawpath, "cgp_drug_concentration.csv"), stringsAsFactors=FALSE)
  drugconc[drugconc == "" | drugconc == " "] <- NA
  drugconc <- data.frame("drug.name"=toupper(gsub(badchars, "", drugconc[ ,"Compound.Name"])), drugconc)
  if(all(!is.element(drugconc[ , "drug.name"], drugnid[ , "drug.name"]))) { stop("Screening concentration for drugs ithout identifiers!") }
  rownames(drugconc) <- rownames(drugnid)[match(drugconc[ , "drug.name"], drugnid[ , "drug.name"])]
  drugconc <- data.frame("drugid"=rownames(drugconc), drugconc)

  ## combine all drugs
  dix <- sort(unique(c(rownames(druginfo), rownames(drugconc), paste("drugid", sapply(strsplit(colnames(drugpheno)[grep("^drugid_", colnames(drugpheno))], "_"), function(x) { return(x[[2]]) }), sep="_"))))
  ## update druginfo
  druginfo2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(druginfo), dimnames=list(dix, colnames(druginfo))))
  newlev <- sapply(druginfo, levels)
  newlev$drugid <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
  druginfo2 <- setcolclass.df(df=druginfo2, colclass=sapply(druginfo, class), factor.levels=newlev)
  druginfo2[match(rownames(druginfo), dix), colnames(druginfo)] <- druginfo
  druginfo2[ , "drugid"] <- newlev$drugid
  druginfo <- druginfo2
  ## update drugconc
  drugconc2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(drugconc), dimnames=list(dix, colnames(drugconc))))
  newlev <- sapply(drugconc, levels)
  newlev$drugid <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
  drugconc2 <- setcolclass.df(df=drugconc2, colclass=sapply(drugconc, class), factor.levels=newlev)
  drugconc2[match(rownames(drugconc), dix), colnames(drugconc)] <- drugconc
  drugconc2[ , "drugid"] <- newlev$drugid
  drugconc <- drugconc2

  ## report concentrations per cell line and per drug
  drugconc2 <- data.frame(matrix(NA, nrow=nrow(drugconc) * length(cellnall), ncol=6, dimnames=list(paste(rep(cellnall, times=nrow(drugconc)), rep(rownames(drugconc), each=length(cellnall)), sep="..."), c("cellid", "drugid", "drug.name", "nbr.conc.tested", "min.Dose.uM", "max.Dose.uM"))))
  drugconc2[ , "cellid"] <- rep(cellnall, times=nrow(drugconc))
  drugconc2[ , "drugid"] <- rep(rownames(drugconc), each=length(cellnall))
  drugconc2[ , "drug.name"] <- rep(as.character(drugconc[ ,"drug.name"]), each=length(cellnall))
  ## as mentioned in the supplementary information of Garnett et al., a single cell line is used on each plate and treated with 28 different drugs over a 9-pt, 256-fold concentration range
  drugconc2[ , "nbr.conc.tested"] <- 9
  drugconc2[ , "min.Dose.uM"] <- rep(drugconc[ , "Min.Concentration.micromolar."], each=length(cellnall))
  drugconc2[ , "max.Dose.uM"] <- rep(drugconc[ , "Max.Concentration.micromolar."], each=length(cellnall))
  drugconc <- drugconc2
  
  # stop()
  # ## keep only the CEL files present in sampleinfo
  # myx <- is.element(names(celfns), sampleinfo[ ,"samplename"])
  # celfn <- celfn[myx]
  # celfns <- celfns[myx]

  ## normalization
  message("Normalize gene expression data")
  # rr <- just.rma(filenames=celfn)
  ## frma normalization using parallel
  splitix <- parallel::splitIndices(nx=length(celfn), ncl=nbcore)
  splitix <- splitix[sapply(splitix, length) > 0]
  res <- parallel::mclapply(splitix, function(x, celfn) {
    ## fRMA
    tt <- celfn[x]
    names(tt) <- NULL
    abatch <- affy::read.affybatch(filenames=tt)
    rr <- frma::frma(object=abatch, summarize="robust_weighted_average", verbose=TRUE, input.vecs=hthgu133afrmavecs)
    rr <- exprs(rr)
  }, celfn=celfn)
  datat <- t(do.call(cbind, res))

  ## build annotation matrix
  message("Build annotation matrix")
  require(biomaRt) || stop("Library biomaRt is not available!")
  mart.db <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="hsapiens_gene_ensembl")
  ## select the best probe for a single gene
  require(jetset) || stop("Library jetset is not available!")
  js <- jetset::jscores(chip="hgu133a", probeset=colnames(datat))
  js <- js[colnames(datat), , drop=FALSE]
  ## identify the best probeset for each entrez gene id
  geneid1 <- as.character(js[ ,"EntrezID"])
  names(geneid1) <- rownames(js)
  geneid2 <- sort(unique(geneid1))
  names(geneid2) <- paste("geneid", geneid2, sep=".")
  gix1 <- !is.na(geneid1)
  gix2 <- !is.na(geneid2)
  geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
  ## probes corresponding to common gene ids
  gg <- names(geneid1)[is.element(geneid1, geneid.common)]
  gid <- geneid1[is.element(geneid1, geneid.common)]
  ## duplicated gene ids
  gid.dupl <- unique(gid[duplicated(gid)])
  gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
  ## unique gene ids
  gid.uniq <- gid[!is.element(gid, gid.dupl)]
  gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
  ## which are the best probe for each gene
  js <- data.frame(js, "best"=FALSE)
  js[gg.uniq, "best"] <- TRUE
  ## data for duplicated gene ids
  if(length(gid.dupl) > 0) {	
  	library(jetset)
  	## use jetset oevrall score to select the best probeset
  	myscore <- js[gg.dupl,"overall"]
  	myscore <- cbind("probe"=gg.dupl, "gid"=geneid1[gg.dupl], "score"=myscore)
  	myscore <- myscore[order(as.numeric(myscore[ , "score"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
  	myscore <- myscore[!duplicated(myscore[ , "gid"]), , drop=FALSE]
  	js[myscore[ ,"probe"], "best"] <- TRUE
  }
  ## more annotations from biomart
  ugid <- sort(unique(js[ ,"EntrezID"]))
  ss <- "entrezgene"
  gene.an <- biomaRt::getBM(attributes=c(ss, "ensembl_gene_id", "hgnc_symbol", "unigene", "description", "chromosome_name", "start_position", "end_position", "strand", "band"), filters="entrezgene", values=ugid, mart=mart.db)
  gene.an[gene.an == "" | gene.an == " "] <- NA
  gene.an <- gene.an[!is.na(gene.an[ , ss]) & !duplicated(gene.an[ , ss]), , drop=FALSE]
  gene.an <- gene.an[is.element(gene.an[ ,ss], ugid), ,drop=FALSE]
  annot <- data.frame(matrix(NA, nrow=ncol(datat), ncol=ncol(gene.an)+1, dimnames=list(colnames(datat), c("probe", colnames(gene.an)))))
  annot[match(gene.an[ , ss], js[ ,"EntrezID"]), colnames(gene.an)] <- gene.an
  annot[ ,"probe"] <- colnames(datat)
  colnames(js)[colnames(js) != "best"] <- paste("jetset", colnames(js)[colnames(js) != "best"], sep=".")
  annot <- data.frame(annot, "EntrezGene.ID"=js[ ,"jetset.EntrezID"], js)
  
  ## save the full dataset, with duplicates
  sampleinfo.full.cgp <- sampleinfo2
  data.full.cgp <- datat
  rownames(data.full.cgp) <- gsub(".CEL.gz", "", rownames(data.full.cgp))
  data.full.cgp <- data.full.cgp[match(rownames(sampleinfo.full.cgp), rownames(data.full.cgp)), , drop=FALSE]
  annot.full.cgp <- annot
  save(list=c("data.full.cgp", "annot.full.cgp", "sampleinfo.full.cgp"), compress=TRUE, file=file.path(saveres, "cgp_full_frma.RData"))
  
  ## match the experiment labels
  myx <- rownames(sampleinfo)[match(rownames(datat), as.character(sampleinfo[ ,"filename"]))]
  datat <- datat[!is.na(myx), , drop=FALSE]
  myx <- myx[!is.na(myx)]
  rownames(datat) <- myx

  ## keep only experiments for which we have all the info
  myx <- fold(intersect, rownames(drugpheno), rownames(sampleinfo), rownames(datat))
  data.ge.cgp <- datat[myx, ,drop=FALSE]
  annot.ge.cgp <- annot
  sampleinfo.ge.cgp <- sampleinfo[myx, , drop=FALSE]
  drugpheno.ge.cgp <- drugpheno[myx, , drop=FALSE]
  druginfo.ge.cgp <- druginfo
  drugconc.ge.cgp <- drugconc[is.element(drugconc[ ,"cellid"], myx), , drop=FALSE]
  mutation.cgp <- mutation[myx, , drop=FALSE]

  ## make sure that cellid are not factors
  celline.cgp[, "cellid"] <- as.character(celline.cgp[, "cellid"])
  sampleinfo.ge.cgp[, "cellid"] <- as.character(sampleinfo.ge.cgp[, "cellid"])
  drugpheno.ge.cgp[, "cellid"] <- as.character(drugpheno.ge.cgp[, "cellid"])
  drugconc.ge.cgp[, "cellid"] <- as.character(drugconc.ge.cgp[, "cellid"])

  message("Save data")
  write.csv(celline.cgp, file=file.path(saveres, "cell_line_collection_cgp.csv"), row.names=FALSE)
  write.csv(annot.ge.cgp, file=file.path(saveres, "annot_ge_cgp.csv"), row.names=FALSE)
  write.csv(t(data.ge.cgp), file=file.path(saveres, "data_ge_cgp.csv"))
  write.csv(drugpheno.ge.cgp, file=file.path(saveres, "drugpheno_ge_cgp.csv"), row.names=FALSE)
  write.csv(druginfo.ge.cgp, file=file.path(saveres, "druginfo_ge_cgp.csv"), row.names=FALSE)
  write.csv(drugconc.ge.cgp, file=file.path(saveres, "drugconc_ge_cgp.csv"), row.names=FALSE)
  write.csv(sampleinfo.ge.cgp, file=file.path(saveres, "sampleinfo_ge_cgp.csv"), row.names=FALSE)
  write.csv(mutation.cgp, file=file.path(saveres, "mutation_cgp.csv"), row.names=TRUE)
  save(list=c("data.ge.cgp", "annot.ge.cgp", "sampleinfo.ge.cgp", "mutation.cgp", "drugpheno.ge.cgp", "druginfo.ge.cgp", "drugconc.ge.cgp", "celline.cgp"), compress=TRUE, file=myfn)
}





## end




