########################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## February 11, 2013
########################


# rm(list=ls())

datapath <- file.path("data", "CCLE")

rawpath <- file.path(datapath, "raw_ge")
if(!file.exists(rawpath)) { dir.create(rawpath, showWarnings=FALSE, recursive=TRUE) }

require(affy) || stop("Library affy is not available!")
require(R.utils) || stop("Library R.utils is not available!")
  
########################
## download data
########################
myfn <- file.path(rawpath, "celfile_timestamp.RData")
if(!file.exists(myfn)) {
  message("Download genomic data")
  dir.create(file.path(rawpath, "dwl"), showWarnings=FALSE)
  
  ## download and compress CEL files
  dwl.status <- download.file(url="http://www.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_21/CCLE_Expression.Arrays_2013-03-18.tar.gz?downloadff=true&fileId=9619", destfile=file.path(rawpath, "dwl", "CCLE_Expression.Arrays_2013-03-18.tar.gz"))
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  archivn <- "CCLE_Expression.Arrays_2013-03-18"
  ## untar
  res <- untar(tarfile=file.path(rawpath, "dwl", sprintf("%s.tar.gz", archivn)), exdir=rawpath)
  
  fff <- affy::list.celfiles(file.path(rawpath, archivn))
  res <- apply(cbind("from"=file.path(rawpath, archivn, fff), "to"=file.path(rawpath, fff)), 1, function(x) { return(file.rename(from=x[1], to=x[2])) })
  unlink(file.path(rawpath, archivn), recursive=TRUE)
  
  celfile.timestamp <- as.character(file.info(file.path(rawpath, fff))[, "mtime"])
  celfile.timestamp <- t(sapply(strsplit(celfile.timestamp, split=" "), function(x) { return(x) }))
  dimnames(celfile.timestamp) <- list(fff, c("file.day", "file.hour"))
  
  ## compress each CEL file individually using gzip
  rr <- sapply(file.path(rawpath, fff), R.utils::gzip, overwrite=TRUE)
   
  # unlink(file.path(rawpath, "dwl"), recursive=TRUE)
  write.csv(celfile.timestamp, file=file.path(rawpath, "celfile_timestamp.csv"))
  save(list=c("celfile.timestamp"), compress=TRUE, file=myfn)
}

## cell line annotations
message("Download cell line annotation")
dwl.status <- download.file(url="http://www.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_22/CCLE_sample_info_file_2012-10-18.txt?downloadff=true&fileId=6801", destfile=file.path(rawpath, "dwl", "CCLE_sample_info_file_2012-10-18.txt"))
if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
file.copy(from=file.path(rawpath, "dwl", "CCLE_sample_info_file_2012-10-18.txt"), to=file.path(rawpath, "ccle_sample_info_file.txt"))
  
 ## drug info
message("Download drug information")
dwl.status <- download.file(url="http://www.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_27/CCLE_NP24.2009_profiling_2012.02.20.csv?downloadff=true&fileId=3422", destfile=file.path(rawpath, "dwl", "CCLE_NP24.2009_profiling_2012.02.20.csv"))
if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
file.copy(from=file.path(rawpath, "dwl", "CCLE_NP24.2009_profiling_2012.02.20.csv"), to=file.path(rawpath, "ccle_drug_info_file.csv"))
  
 ## drug pheno
message("Download drug sensitivity measurements")
## drug sensitivity data from the addendum in Nature
dwl.status <- download.file(url="http://www.nature.com/nature/journal/v492/n7428/extref/nature11735-s2.xls", destfile=file.path(rawpath, "dwl", "nature11735-s2.xls"))
if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
file.copy(from=file.path(rawpath, "dwl", "nature11735-s2.xls"), to=file.path(rawpath, "ccle_drug_pheno_file.xls"))
## drug sensitivity data from the CCLE website
# dwl.status <- download.file(url="http://www.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_27/CCLE_NP24.2009_Drug_data_2012.02.20.csv?downloadff=true&fileId=2114", destfile=file.path(rawpath, "dwl", "CCLE_NP24.2009_Drug_data_2012.02.20.csv"))
# if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
# file.copy(from=file.path(rawpath, "dwl", "CCLE_NP24.2009_Drug_data_2012.02.20.csv"), to=file.path(rawpath, "ccle_drug_pheno_file.csv"))

## mutations
message("Download mutation data")
## hybrid capture
dwl.status <- download.file(url="http://www.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_26/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf.gz?downloadff=true&fileId=6873", destfile=file.path(rawpath, "dwl", "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf.gz"))
if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
gunzip(filename=file.path(rawpath, "dwl", "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf.gz"), overwrite=TRUE)
file.copy(from=file.path(rawpath, "dwl", "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf"), to=file.path(rawpath, "ccle_mutations_hybrid.maf"))
## oncomap3
dwl.status <- download.file(url="http://www.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_23/CCLE_Oncomap3_2012-04-09.maf?downloadff=true&fileId=3000", destfile=file.path(rawpath, "dwl", "CCLE_Oncomap3_2012-04-09.maf"))
if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
file.copy(from=file.path(rawpath, "dwl", "CCLE_Oncomap3_2012-04-09.maf"), to=file.path(rawpath, "ccle_mutations_oncomap.maf"))

########################
## normalizae and format data
########################
myfn <- file.path(saveres, "ccle_frma.RData")
if(!file.exists(myfn)) {
  
  require(gdata) || stop("Library gdata is not available!")
  require(affy) || stop("Library affy is not available!")
  require(Hmisc) || stop("Library Hmisc is not available!")
  require(genefu) || stop("Library genefu is not available!")
  require(frma) || stop("Library frma is not available!")
  require(hgu133plus2frmavecs) || stop("Library hgu133plus2frmavecs is not available!")
  data(hgu133plus2frmavecs)
  require(hgu133plus2cdf) || stop("Library hgu133plus2cdf is not available!")
  data(hgu133plus2cdf)

  ## CEL file names
  celfn <- list.celfiles(rawpath, full.names=TRUE)
  celfns <- list.celfiles(rawpath, full.names=FALSE)
  ## experiments' names
  names(celfn) <- names(celfns) <- gsub(".CEL.gz", "", celfns)
  ## chip type and date
  chipt <- sapply(celfn, celfileChip)
  chipd <- t(sapply(celfn, celfileDateHour))
  ## reorder CEL files by hybridization time or timestamp
  load(file.path(rawpath, "celfile_timestamp.RData"))
  myx <- NULL
  if(any(!complete.cases(chipd))) {
    ## all hybridization dates are not available
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

  ## drug information
  message("Read drug information")
  druginfo <- read.csv(file.path(rawpath, "ccle_drug_info_file.csv"), stringsAsFactors=FALSE)
  druginfo[druginfo == "" | druginfo == " "] <- NA
  ## manual drug name curation
  druginfo[druginfo[ ,"Compound..code.or.generic.name."] == "Panobinostat\xa0\xa0","Compound..code.or.generic.name."] <- "Panobinostat"
  druginfo[druginfo[ ,"Compound..code.or.generic.name."] == "PF-2341066","Compound..code.or.generic.name."] <- "Crizotinib"
  druginfo[druginfo[ ,"Compound..code.or.generic.name."] == "PD-0332991 ","Compound..code.or.generic.name."] <- "PD-0332991"
  druginfo <- data.frame(druginfo, "drugid"=paste("drugid", toupper(gsub(pattern=badchars, replacement="", x=toupper(druginfo[ ,"Compound..code.or.generic.name."]))), sep="_"))
  rownames(druginfo) <- druginfo[ ,"drugid"]

  ## phenotype for the drugs
  message("Read drug sensitivity measurements")
  ## drug sensitivity data from the addendum in Nature
  myfn2 <- file.path(saveres, "ccle_drug_pheno_file.RData")
  if(!file.exists(myfn2)) {
    drugpheno <- gdata::read.xls(xls=file.path(rawpath, "ccle_drug_pheno_file.xls"), sheet=12)
    drugpheno[drugpheno == "" | drugpheno == " "] <- NA
    save(list="drugpheno", compress=TRUE, file=myfn2)
  } else { load(myfn2) }
  nn <- apply(drugpheno[1, ], 2, as.character)
  nn <- gsub(pattern=badchars, replacement="_", x=nn)
  drugpheno <- drugpheno[-c(1),!is.na(nn), drop=FALSE]
  nn <- nn[!is.na(nn)]
  nn[match(c("Primary_Cell_Line_Name", "IC50_(_M)_(norm)", "EC50_(_M)_(norm)", "Amax_(norm)" ,"ActArea_(norm)", "Doses_(uM)", "Activity_Data_(normalized_data)", "Activity_Error_(SD)_(normalized_data)", "FitType_(norm)"), nn)] <- c("Primary.Cell.Line.Name", "IC50..uM.",  "EC50..uM.", "Amax", "ActArea", "Doses..uM.", "Activity.Data..median.", "Activity.SD", "FitType")
  colnames(drugpheno) <- nn
  drugpheno <- setcolclass.df(df=drugpheno, colclass=rep("character", ncol(drugpheno)), factor.levels=sapply(drugpheno, levels))
  drugpheno <- setcolclass.df(df=drugpheno, colclass=c(rep("character", 8), rep("numeric", 4), "character", "numeric", "numeric"))
  ## drug sensitivity data from the CCLE website
  # drugpheno <- read.csv(file.path(rawpath, "ccle_drug_pheno_file.csv"), stringsAsFactors=FALSE)
  # colnames(drugpheno) <- gsub(pattern=badchars, replacement="_", x=colnames(drugpheno))
  # drugpheno[drugpheno == "" | drugpheno == " "] <- NA
  ## manual curation of drug names
  drugpheno[drugpheno[ ,"Compound"] == "ZD-6474", "Compound"] <- "Vandetanib"
  drugpheno[drugpheno[ ,"Compound"] == "PF2341066", "Compound"] <- "Crizotinib"
  drugpheno[ ,"Compound"] <- toupper(gsub(pattern=badchars, replacement="", drugpheno[ ,"Compound"]))

  drugpheno <- data.frame(drugpheno, "drugid"=paste("drugid", toupper(gsub(pattern=badchars, replacement="", x=drugpheno[ ,"Compound"])), sep="_"), "cellid"=drugpheno[ ,"Primary.Cell.Line.Name"])
  ## drug screening concentrations
  ll <- sapply(strsplit(drugpheno[ , "Doses..uM."], ","), length)
  drugconc <- lapply(strsplit(drugpheno[ , "Doses..uM."], ","), function(x) {
    xx <- as.numeric(x)
    xx2 <- c(0.0025, 0.0080, 0.0250, 0.0800, 0.2500, 0.8000, 2.5300, 8.0000)
    if(any(!is.element(xx, xx2))) { stop("Unexpected drug screening concentrations!") }
    xx3 <- rep(NA, length(xx2))
    names(xx3) <- xx2
    xx3[match(xx, names(xx3))] <- xx
    return(xx3)
    })
  drugconc <- do.call(rbind, drugconc)
  drugconc <- data.frame("cellid"=as.character(drugpheno[ , "cellid"]), "drugid"=as.character(drugpheno[ , "drugid"]), "nbr.conc.tested"=ll, drugconc)
  dimnames(drugconc) <- list(paste(drugconc[ , "cellid"], drugconc[ , "drugid"], sep="..."), c("cellid", "drugid", "nbr.conc.tested", sprintf("Dose%i.uM", 1:(ncol(drugconc) - 3))))

  ## combine all drugs
  dix <- sort(unique(c(as.character(druginfo[ ,"drugid"]), as.character(drugpheno[ ,"drugid"]))))
  ## update druginfo
  druginfo2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(druginfo), dimnames=list(dix, colnames(druginfo))))
  newlev <- sapply(druginfo, levels)
  newlev$drugid <- dix
  druginfo2 <- setcolclass.df(df=druginfo2, colclass=sapply(druginfo, class), factor.levels=newlev)
  druginfo2[match(as.character(druginfo[ ,"drugid"]), dix), colnames(druginfo)] <- druginfo
  druginfo2[ , "drugid"] <- newlev$drugid
  druginfo <- druginfo2

  ## reorder the drug phenotypes per cell lines
  ## should be done for Doses..uM., Activity.Data..median., Activity.SD, Num.Data, FitType, EC50..uM., IC50..uM., Amax, ActArea
  celln <- unique(as.character(drugpheno[ , "cellid"]))
  drugpheno.all <- NULL
  ## IC50
  drugphenot <- matrix(NA, nrow=length(celln), ncol=length(dix), dimnames=list(as.character(celln), dix))
  for(i in 1:length(celln)) {
      tt <- as.numeric(drugpheno[drugpheno[ ,"cellid"] == celln[i] ,"IC50..uM."])
      names(tt) <- as.character(drugpheno[drugpheno[ ,"cellid"] == celln[i] ,"drugid"])
      drugphenot[celln[i],names(tt)] <- tt
  }
  colnames(drugphenot) <- paste(colnames(drugphenot), "IC50", sep="_")
  if(is.null(drugpheno.all)) { drugpheno.all <- drugphenot } else { drugpheno.all <- data.frame(drugpheno.all, drugphenot) }
  ## EC50..uM.
  drugphenot <- matrix(NA, nrow=length(celln), ncol=length(dix), dimnames=list(as.character(celln), dix))
  for(i in 1:length(celln)) {
      tt <- as.numeric(drugpheno[drugpheno[ ,"cellid"] == celln[i] ,"EC50..uM."])
      names(tt) <- as.character(drugpheno[drugpheno[ ,"cellid"] == celln[i] ,"drugid"])
      drugphenot[celln[i],names(tt)] <- tt
  }
  colnames(drugphenot) <- paste(colnames(drugphenot), "EC50", sep="_")
  if(is.null(drugpheno.all)) { drugpheno.all <- drugphenot } else { drugpheno.all <- data.frame(drugpheno.all, drugphenot) }
  ## Amax
  drugphenot <- matrix(NA, nrow=length(celln), ncol=length(dix), dimnames=list(as.character(celln), dix))
  for(i in 1:length(celln)) {
      tt <- as.numeric(drugpheno[drugpheno[ ,"cellid"] == celln[i] ,"Amax"])
      names(tt) <- as.character(drugpheno[drugpheno[ ,"cellid"] == celln[i] ,"drugid"])
      drugphenot[celln[i],names(tt)] <- tt
  }
  colnames(drugphenot) <- paste(colnames(drugphenot), "Amax", sep="_")
  if(is.null(drugpheno.all)) { drugpheno.all <- drugphenot } else { drugpheno.all <- data.frame(drugpheno.all, drugphenot) }
  ## ActArea
  drugphenot <- matrix(NA, nrow=length(celln), ncol=length(dix), dimnames=list(as.character(celln), dix))
  for(i in 1:length(celln)) {
      tt <- as.numeric(drugpheno[drugpheno[ ,"cellid"] == celln[i] ,"ActArea"])
      names(tt) <- as.character(drugpheno[drugpheno[ ,"cellid"] == celln[i] ,"drugid"])
      drugphenot[celln[i],names(tt)] <- tt
  }
  colnames(drugphenot) <- paste(colnames(drugphenot), "ActivityArea", sep="_")
  if(is.null(drugpheno.all)) { drugpheno.all <- drugphenot } else { drugpheno.all <- data.frame(drugpheno.all, drugphenot) }
  ## Dose
  drugphenot <- matrix(NA, nrow=length(celln), ncol=length(dix), dimnames=list(as.character(celln), dix))
  for(i in 1:length(celln)) {
      tt <- drugpheno[drugpheno[ ,"cellid"] == celln[i] ,"Doses..uM."]
      names(tt) <- as.character(drugpheno[drugpheno[ ,"cellid"] == celln[i] ,"drugid"])
      drugphenot[celln[i],names(tt)] <- tt
  }
  colnames(drugphenot) <- paste(colnames(drugphenot), "Doses", sep="_")
  if(is.null(drugpheno.all)) { drugpheno.all <- drugphenot } else { drugpheno.all <- data.frame(drugpheno.all, drugphenot) }
  ## Activity.Data..median
  drugphenot <- matrix(NA, nrow=length(celln), ncol=length(dix), dimnames=list(as.character(celln), dix))
  for(i in 1:length(celln)) {
      tt <- drugpheno[drugpheno[ ,"cellid"] == celln[i] ,"Activity.Data..median."]
      names(tt) <- as.character(drugpheno[drugpheno[ ,"cellid"] == celln[i] ,"drugid"])
      drugphenot[celln[i],names(tt)] <- tt
  }
  colnames(drugphenot) <- paste(colnames(drugphenot), "ActivityDataMedian", sep="_")
  if(is.null(drugpheno.all)) { drugpheno.all <- drugphenot } else { drugpheno.all <- data.frame(drugpheno.all, drugphenot) }
  ## Activity.SD
  drugphenot <- matrix(NA, nrow=length(celln), ncol=length(dix), dimnames=list(as.character(celln), dix))
  for(i in 1:length(celln)) {
      tt <- drugpheno[drugpheno[ ,"cellid"] == celln[i] ,"Activity.SD"]
      names(tt) <- as.character(drugpheno[drugpheno[ ,"cellid"] == celln[i] ,"drugid"])
      drugphenot[celln[i],names(tt)] <- tt
  }
  colnames(drugphenot) <- paste(colnames(drugphenot), "ActivitySD", sep="_")
  if(is.null(drugpheno.all)) { drugpheno.all <- drugphenot } else { drugpheno.all <- data.frame(drugpheno.all, drugphenot) }
  ## FitType
  drugphenot <- matrix(NA, nrow=length(celln), ncol=length(dix), dimnames=list(as.character(celln), dix))
  for(i in 1:length(celln)) {
      tt <- drugpheno[drugpheno[ ,"cellid"] == celln[i] ,"FitType"]
      names(tt) <- as.character(drugpheno[drugpheno[ ,"cellid"] == celln[i] ,"drugid"])
      drugphenot[celln[i],names(tt)] <- tt
  }
  colnames(drugphenot) <- paste(colnames(drugphenot), "FitType", sep="_")
  if(is.null(drugpheno.all)) { drugpheno.all <- drugphenot } else { drugpheno.all <- data.frame(drugpheno.all, drugphenot) }
  ## save the new spreadsheet
  drugpheno <- data.frame("cellid"=rownames(drugpheno.all), drugpheno.all)
  
  ## info about each experiment
  message("Read sample information")
  sampleinfo <- read.csv(file.path(rawpath, "ccle_sample_info_file.txt"), sep="\t", stringsAsFactors=FALSE)
  sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
  sampleinfo <- data.frame("cellid"=as.character(sampleinfo[ ,"Cell.line.primary.name"]), sampleinfo)
  ## remove duplicated cell line hybridization
  ## only the most recent experiment (as determine by hyridization date or CEL file timestamp) will be kept for each cell line
  sampleinfo <- sampleinfo[!duplicated(sampleinfo[ ,"cellid"]), ,drop=FALSE]
  sampleinfo[ , "cellid"] <- as.character(sampleinfo[ , "cellid"])
  rownames(sampleinfo) <- as.character(sampleinfo[ , "cellid"])

  ## extract protein coding variants from mutation data
  myfn2 <- file.path(saveres, "ccle_mutations.RData")
  if(!file.exists(myfn2)) {
    message("Read (missense) mutation data")
    ## from hybrid capture
    mut <- read.csv(file.path(rawpath, "ccle_mutations_hybrid.maf"), stringsAsFactors=FALSE, sep="\t")
    mut <- mut[ , c("Hugo_Symbol", "Tumor_Sample_Barcode", "Protein_Change"), drop=FALSE]
    mut[!is.na(mut) & mut == ""] <- NA
    mut[is.na(mut[ , "Protein_Change"]) | mut[ , "Protein_Change"] == "", "Protein_Change"] <- "wt"
    mut[!is.na(mut[ , "Protein_Change"]) & (mut[ , "Protein_Change"] == "p.?" | mut[ , "Protein_Change"] == "p.0?"), "Protein_Change"] <- NA
    mut <- mut[complete.cases(mut), , drop=FALSE]
    myx <- !duplicated(paste(mut[ , c("Tumor_Sample_Barcode")], mut[ , c("Hugo_Symbol")], mut[ , c("Protein_Change")], sep="///"))
    mut <- mut[myx, , drop=FALSE]
    ## from oncomap
    mut2 <- read.csv(file.path(rawpath, "ccle_mutations_oncomap.maf"), stringsAsFactors=FALSE, sep="\t")
    mut2 <- mut2[ , c("Hugo_Symbol", "Tumor_Sample_Barcode", "Protein_Change"), drop=FALSE]
    mut2[!is.na(mut2) & mut2 == ""] <- NA
    # mut2[is.na(mut2[ , "Protein_Change"]) | mut2[ , "Protein_Change"] == "", "Protein_Change"] <- "wt"
    mut2[!is.na(mut2[ , "Protein_Change"]) & (mut2[ , "Protein_Change"] == "p.?" | mut2[ , "Protein_Change"] == "p.0?"), "Protein_Change"] <- NA  
    mut2 <- mut2[complete.cases(mut2), , drop=FALSE]
    myx <- !duplicated(paste(mut2[ , c("Tumor_Sample_Barcode")], mut2[ , c("Hugo_Symbol")], mut2[ , c("Protein_Change")], sep="///"))
    mut2 <- mut2[myx, , drop=FALSE]
    ## merge mutations discovered by hybrid capture and oncomap
    mutation <- rbind(mut, mut2)
    ucell <- sort(unique(mutation[ , "Tumor_Sample_Barcode"]))
    ugene <- sort(unique(mutation[ , "Hugo_Symbol"]))

    # ## create a list (cell lines) of lists (genes) of vectors (mutations)
    # splitix <- parallel::splitIndices(nx=length(ucell), ncl=nbcore)
    # system.time(res <- parallel::mclapply(splitix, function(x, ucell, ugene, mutation) {
    #   ucell2 <- ucell[x]
    #   dd <- lapply(ucell2, function(x, y) {
    #     ll <- lapply(y, function(x) { return("wt") })
    #     names(ll) <- y
    #     return(ll)
    #   }, y=ugene)
    #   names(dd) <- ucell2
    #   for(uc in 1:length(ucell2)) {
    #     # message(sprintf("Process mutation for cell line %s", ucell2[uc]))
    #     ucix <- mutation[ , "Tumor_Sample_Barcode"] == ucell2[uc]
    #     for(ug in 1:length(ugene)) {
    #       myx <- ucix & mutation[ , "Hugo_Symbol"] == ugene[ug]
    #       if(any(myx)) { dd[[ucell2[uc]]][[ugene[ug]]] <- unique(mutation[myx, "Protein_Change"]) }
    #     }
    #   }
    #   return(dd)
    # }, ucell=ucell, ugene=ugene, mutation=mutation))
    # res <- do.call(c, res)
    # mutation <- res
  
    dd <- matrix("wt", nrow=length(ucell), ncol=length(ugene), dimnames=list(ucell, ugene))
    ## maximum number of mutations per gene
    # mutmax <- max(table(paste(mutation[, c("Tumor_Sample_Barcode")], mutation[, c("Hugo_Symbol")], sep="::")))
    mm <- 1:nrow(mutation)
    ff <- TRUE
    while(length(mm) > 1) {
      myx <- !duplicated(paste(mutation[mm, c("Tumor_Sample_Barcode")], mutation[mm, c("Hugo_Symbol")], sep="///"))
      if(ff) {
        dd[as.matrix(mutation[mm[myx], c("Tumor_Sample_Barcode", "Hugo_Symbol")])] <- mutation[mm[myx], "Protein_Change"]
        ff <- FALSE
      } else {
          dd[as.matrix(mutation[mm[myx], c("Tumor_Sample_Barcode", "Hugo_Symbol")])] <- paste(dd[as.matrix(mutation[mm[myx], c("Tumor_Sample_Barcode", "Hugo_Symbol")])], mutation[mm[myx], "Protein_Change"], sep="///")
      }
      mm <- mm[!myx]
    }
    ## check for inconsistencies (wt + mutations)
    iix <- grep("///", dd)
    for(iii in iix) {
      x <- sort(unique(unlist(strsplit(dd[iii], split="///"))))
      if(length(x) > 1) { s <- x[!is.element(x, "wt")] }
      dd[iii] <- paste(x, collapse="///")
    }
    nn <- sampleinfo[match(rownames(dd), sampleinfo[ , "CCLE.name"]), "cellid"]
    ## remove if we do not have cell line identifier
    dd <- dd[!is.na(nn), , drop=FALSE]
    rownames(dd) <- nn[!is.na(nn)]
    mutation <- dd
    save(list=c("mutation"), compress=TRUE, file=myfn2)
  } else { load(myfn2) }

  iix <- which(!is.element(rownames(drugpheno), rownames(sampleinfo)))
  celline.ccle <- data.frame(matrix(NA, nrow=nrow(sampleinfo) + length(iix), ncol=ncol(sampleinfo), dimnames=list(c(rownames(sampleinfo), rownames(drugpheno)[iix]), colnames(sampleinfo))))
  celline.ccle[rownames(sampleinfo), colnames(sampleinfo)] <- sampleinfo
  celline.ccle[rownames(drugpheno)[iix], "Cell.line.primary.name"] <- rownames(drugpheno)[iix]
  celline.ccle[rownames(drugpheno)[iix], "cellid"] <- rownames(drugpheno)[iix]
  ## add url based on CCLE name
  uurl <- paste("http://www.broadinstitute.org/ccle/cell%20lines/", celline.ccle[ , "CCLE.name"], sep="")
  uurl[is.na(celline.ccle[ , "CCLE.name"])] <- NA
  celline.ccle <- data.frame("cellid"=celline.ccle[ , "cellid"], "link"=uurl, celline.ccle[ , !is.element(colnames(celline.ccle), "cellid")])

  # ## rename a problematic cell line names
  # sampleinfo[sampleinfo[ ,"Cell.line.primary.name"] == "T.T", "Cell.line.primary.name"] <- "TT_OESOPHAGUS"
  # ## cell line identifiers
  # dupln <- sum(duplicated(sampleinfo[ ,"Cell.line.primary.name"]))
  # celln2 <- gsub(pattern =badchars, replacement="", x=toupper(as.character(sampleinfo[ ,"Cell.line.primary.name"])))
  # dupln <- c(dupln, sum(duplicated(celln2)))
  # ## remove specific strings from cell line names
  # celln2 <- gsub(pattern ="NCI", replacement="", x=celln2)
  # dupln <- c(dupln, sum(duplicated(celln2)))
  # celln2 <- gsub(pattern ="NIH", replacement="", x=celln2)
  # dupln <- c(dupln, sum(duplicated(celln2)))
  # if(any(dupln[-1] > dupln[1])) { stop("Duplicated cell line identifiers due to curation of sample information!") }


  ## keep only the CEL files present in sampleinfo
  fn <- sampleinfo[ ,"Expression.arrays"]
  # if(any(!is.element(fn[!is.na(fn)], names(celfns)))) { stop("some CEL files are missing for the CCLE project") }
  myx <- intersect(fn[!is.na(fn)], names(celfns))
  celfn <- celfn[myx]
  celfns <- celfns[myx]
  chipt <- chipt[myx]
  chipd <- chipd[myx, , drop=FALSE]
  celfile.timestamp <- celfile.timestamp[paste(myx, ".CEL", sep=""), , drop=FALSE]
  sampleinfo <- sampleinfo[match(myx, sampleinfo[ ,"Expression.arrays"]), ,drop=FALSE]
  sampleinfo <- data.frame("samplename"=names(celfns), "filename"=celfns, "chiptype"=chipt, "hybridization.date"=chipd[ ,"day"], "hybridization.hour"=chipd[ ,"hour"], "file.day"=celfile.timestamp[ ,"file.day"], "file.hour"=celfile.timestamp[ ,"file.hour"], "batch"=NA, sampleinfo)
  rownames(sampleinfo) <- as.character(sampleinfo[ ,"cellid"])

  ## use cell lines as identifiers
  cellnall <- sort(unique(c(as.character(sampleinfo[ ,"cellid"]), as.character(drugpheno[ ,"cellid"]), rownames(mutation))))
  ## update drugpheno
  dd <- data.frame(matrix(NA, ncol=ncol(drugpheno), nrow=length(cellnall), dimnames=list(cellnall, colnames(drugpheno))))
  new_levels <- sapply(drugpheno, levels)
  new_levels$cellid <- cellnall
  dd <- setcolclass.df(df=dd, colclass=sapply(drugpheno, class), factor.levels=new_levels)
  dd[rownames(drugpheno),colnames(drugpheno)] <- drugpheno
  dd[ ,"cellid"] <- cellnall
  drugpheno <- dd
  ## update mutations
  dd <- matrix(NA, ncol=ncol(mutation), nrow=length(cellnall), dimnames=list(cellnall, colnames(mutation)))
  dd[rownames(mutation), colnames(mutation)] <- mutation
  mutation <- dd
  
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
    rr <- frma::frma(object=abatch, summarize="robust_weighted_average", verbose=TRUE, input.vecs=hgu133plus2frmavecs)
    return(exprs(rr))
  }, celfn=celfn)
  datat <- t(do.call(cbind, res))
  ## match the experiment labels
  rownames(datat) <- rownames(sampleinfo)[match(rownames(datat), as.character(sampleinfo[ ,"filename"]))]

  ## build annotation matrix
  message("Build annotation matrix")
  library(biomaRt)
  mart.db <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="hsapiens_gene_ensembl")
  ## select the best probe for a single gene
  library(jetset)
  js <- jetset::jscores(chip="hgu133plus2", probeset=colnames(datat))
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

  ## keep only experiments for which we have all the info
  myx <- fold(intersect, rownames(drugpheno), rownames(sampleinfo), rownames(datat))
  data.ge.ccle <- datat[myx, , drop=FALSE]
  annot.ge.ccle <- annot
  sampleinfo.ge.ccle <- sampleinfo[myx, , drop=FALSE]
  drugpheno.ge.ccle <- drugpheno[myx, , drop=FALSE]
  druginfo.ge.ccle <- druginfo
  drugconc.ge.ccle <- drugconc
  mutation.ccle <- mutation[myx, , drop=FALSE]

  ## make sure that cellid are not factors
  celline.ccle[, "cellid"] <- as.character(celline.ccle[, "cellid"])
  sampleinfo.ge.ccle[, "cellid"] <- as.character(sampleinfo.ge.ccle[, "cellid"])
  drugpheno.ge.ccle[, "cellid"] <- as.character(drugpheno.ge.ccle[, "cellid"])
  drugconc.ge.ccle[, "cellid"] <- as.character(drugconc.ge.ccle[, "cellid"])

  message("Save data")
  write.csv(celline.ccle, file=file.path(saveres, "cell_line_collection_ccle.csv"), row.names=FALSE)
  write.csv(annot.ge.ccle, file=file.path(saveres, "annot_ge_ccle.csv"), row.names=FALSE)
  write.csv(t(data.ge.ccle), file=file.path(saveres, "data_ge_ccle.csv"))
  write.csv(drugpheno.ge.ccle, file=file.path(saveres, "drugpheno_ge_ccle.csv"), row.names=FALSE)
  write.csv(druginfo.ge.ccle, file=file.path(saveres, "druginfo_ge_ccle.csv"), row.names=FALSE)
  write.csv(drugconc.ge.ccle, file=file.path(saveres, "drugconc_ge_ccle.csv"), row.names=FALSE)
  write.csv(sampleinfo.ge.ccle, file=file.path(saveres, "sampleinfo_ge_ccle.csv"), row.names=FALSE)
  write.csv(mutation.ccle, file=file.path(saveres, "mutation_ccle.csv"), row.names=TRUE)
  save(list=c("data.ge.ccle", "annot.ge.ccle", "sampleinfo.ge.ccle", "drugpheno.ge.ccle", "druginfo.ge.ccle", "drugconc.ge.ccle", "mutation.ccle", "celline.ccle"), compress=TRUE, file=myfn)
}




## end
