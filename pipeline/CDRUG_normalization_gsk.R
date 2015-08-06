########################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## February 11, 2013
########################

# rm(list=ls())

datapath <- file.path("data", "GSKCELLINES")
rawpath <- file.path(datapath, "raw_ge")
if(!file.exists(rawpath)) { dir.create(rawpath, showWarnings=FALSE, recursive=TRUE) }
  
require(gdata) || stop("Library gdata is not available!")
require(R.utils) || stop("Library R.utils is not available!")
require(genefu) || stop("Library genefu is not available!")

########################
## download data
########################
ftpdir <- "ftp://caftpd.nci.nih.gov/pub/caARRAY/transcript_profiling"
myfn <- file.path(rawpath, "celfile_timestamp.RData")
if(!file.exists(myfn)) {
  message("Download genomic data")
  
  dir.create(file.path(rawpath, "dwl"), showWarnings=FALSE)
  
  ## download and compress CEL files
  celfile.timestamp <- celfn <- NULL
  i <- 1
  while(i <= 8) {
    ## assuming there are only 9 zip archives (need to check if the update version has more)
   dwl.status <- download.file(url=sprintf("%s/cel/cel_0%i.zip", ftpdir, i), destfile=file.path(rawpath, "dwl", sprintf("cel_0%i.zip", i)))
   if(dwl.status != 0) {
     message("\t-> download failed, let's try again ...")
     file.remove(file.path(rawpath, "dwl", sprintf("cel_0%i.zip", i)))
     i <- i - 1
    } else {
       ## unzip archive
       fff <- unzip(zipfile=file.path(rawpath, "dwl", sprintf("cel_0%i.zip", i)), list=TRUE)
       celfile.timestamp <- c(celfile.timestamp, as.character(fff[ ,"Date"]))
       celfn <- c(celfn, as.character(fff[ ,"Name"]))
       res <- unzip(zipfile=file.path(rawpath, "dwl", sprintf("cel_0%i.zip", i)), exdir=rawpath)
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
dwl.status <- download.file(url=sprintf("%s/GSK_RNA.sdrf", ftpdir), destfile=file.path(rawpath, "dwl", "GSK_RNA.sdrf"))
if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
file.copy(from=file.path(rawpath, "dwl", "GSK_RNA.sdrf"), to=file.path(rawpath, "GSK_RNA.sdrf"))
  
## download drug sensitivity (release 2)
message("Download drug sensitivity measurements")
dwl.status <- download.file(url="http://cancerres.aacrjournals.org/content/suppl/2010/04/19/0008-5472.CAN-09-3788.DC1/stab_2.xls", destfile=file.path(rawpath, "dwl", "stab_2.xls"))
if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
file.copy(from=file.path(rawpath, "dwl", "stab_2.xls"), to=file.path(rawpath, "gsk_drug_sensitivity.xls"))


########################
## normalize and format data
########################
myfn <- file.path(saveres, "gskcellines_frma.RData")
if(!file.exists(myfn)) {

  require(affy) || stop("Library affy is not available!")
  require(Hmisc) || stop("Library Hmisc is not available!")
  require(genefu) || stop("Library genefu is not available!")
  require(frma) || stop("Library frma is not available!")
  require(hgu133plus2frmavecs) || stop("Library hgu133plus2frmavecs is not available!")
  data(hgu133plus2frmavecs)
  require(hgu133plus2cdf) || stop("Library hgu133plus2cdf is not available!")
  data(hgu133plus2cdf)
  
  ## load timestamp
  load(file.path(rawpath, "celfile_timestamp.RData"))

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

  ## info about each experiment
  message("Read sample information")
  sampleinfo <- read.csv(file.path(rawpath, "GSK_RNA.sdrf"), sep="\t", stringsAsFactors=FALSE, comment.char="#")
  sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
  ## curate cell line names
  rownames(sampleinfo) <- gsub(" - Replicate ", "_rep", sampleinfo[ , "Source.Name"])
  sampleinfo <- data.frame("samplename"=rownames(sampleinfo), "cellid"=sampleinfo[ , "Characteristics.Cell.Line.Name."], "filename"=sprintf("%s.gz", sampleinfo[ , "Array.Data.File"]), "tissue.type"=tolower(gsub("_$", "", gsub(badchars, "_", sampleinfo[ , "Characteristics.OrganismPart."]))), sampleinfo)
  
  ## read drug phenotypes
  drugpheno <- read.xls(xls=file.path(rawpath, "gsk_drug_sensitivity.xls"), sheet=1, stringsAsFactors=FALSE)
  drugpheno[drugpheno == ""] <- NA
  cn <- gsub(badchars, ".", drugpheno[6, ])
  cn <- gsub("CL.ID", "CL_ID", cn)
  drugpheno <- drugpheno[-(1:6), , drop=FALSE]
  drugpheno <- drugpheno[!is.na(drugpheno[ , 2]), , drop=FALSE]
  dimnames(drugpheno) <- list(drugpheno[ , 2], cn)
  celline.info <- drugpheno[ , c("CL_ID", "Cell.Line", "Site", "Dx"), drop=FALSE]
  mutation <- drugpheno[ , c("HRAS.", "KRAS.", "NRAS.", "BRAF.", "PIK3CA.", "PTEN."), drop=FALSE]
  drugpheno <- data.matrix(drugpheno[ , grep("IC50", cn), drop=FALSE])
  ## IC50 in nano molar

  ## keep only the CEL files present in sampleinfo
  myx <- intersect(sampleinfo[ , "filename"], celfns)
  celfn <- celfn[match(myx, celfns)]
  celfns <- celfns[match(myx, celfns)]
  sampelinfo <- sampleinfo[match(myx, sampleinfo[ , "filename"]), , drop=FALSE]

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
    rr <- exprs(rr)
  }, celfn=celfn)
  datat <- t(do.call(cbind, res))
  ## match the experiment labels
  rownames(datat) <- rownames(sampleinfo)[match(rownames(datat), as.character(sampleinfo[ ,"filename"]))]

  ## build annotation matrix
  message("Build annotation matrix")
  require(biomaRt) || stop("Library biomaRt is not available!")
  mart.db <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="hsapiens_gene_ensembl")
  ## select the best probe for a single gene
  require(jetset) || stop("Library jetset is not available!")
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

  sampleinfo.ge.gsk <- sampleinfo
  data.ge.gsk <- datat
  annot.ge.gsk <- annot
  drugpheno.gsk <- drugpheno
  mutation.gsk <- mutation
  
  save(list=c("data.ge.gsk", "annot.ge.gsk", "sampleinfo.ge.gsk", "drugpheno.gsk", "mutation.gsk"), compress=TRUE, file=file.path(saveres, "gskcellines_complete_frma.RData"))
  
  ## average replicates
  message("Averaging replicates")
  ucell <- as.character(unique(sampleinfo.ge.gsk[ , "cellid"]))
  splitix <- parallel::splitIndices(nx=length(ucell), ncl=nbcore)
  res <- parallel::mclapply(splitix, function(x, ucell, data) {
    dd <- t(sapply(ucell[x], function(x, y) {
      iix <- grep(x, rownames(y))
      dd <- apply(y[iix, , drop=FALSE], 2, mean, na.rm=FALSE)
      return(dd)
    }, y=data))
    # dd <- matrix(NA, ncol=ncol(data.ge.gsk), nrow=length(ucell), dimnames=list(ucell, colnames(data.ge.gsk)))
    # for(i in 1:length(ucell)) {
    #   iix <- grep(ucell[i], rownames(data.gsk))
    #   dd[ucell[i], ] <- apply(data.ge.gsk[iix, , drop=FALSE], 2, mean, na.rm=FALSE)
    # }
    return(dd)
  }, ucell=ucell, data=data.ge.gsk)
  data.ge.gsk <- do.call(rbind, res)
  sampleinfo.ge.gsk <- sampleinfo.ge.gsk[!duplicated(sampleinfo.ge.gsk[ , "cellid"]), , drop=FALSE]
  rownames(sampleinfo.ge.gsk) <- as.character(sampleinfo.ge.gsk[ , "cellid"])
  sampleinfo.ge.gsk <- sampleinfo.ge.gsk[rownames(data.ge.gsk), , drop=FALSE]
  
  ## match cell line names
  ## no hits with partial matching
  # myx <- !is.element(rownames(sampleinfo.ge.gsk), gsub(badchars, "", rownames(drugpheno.gsk)))
  # myx2 <- !is.element(gsub(badchars, "", rownames(drugpheno.gsk)), rownames(sampleinfo.ge.gsk))
  # tt <- lapply(rownames(sampleinfo.ge.gsk)[myx], agrep, x=rownames(drugpheno.gsk)[myx2], value=TRUE)
  # tt <- mapply(function(x, y) { return(sprintf("%s -> %s", x, paste(y, collapse=":::"))) }, x=as.list(rownames(sampleinfo.ge.gsk)[myx]), y=tt)
  # write.csv(tt, "temp.csv")
  drugpheno.ge.gsk <- matrix(NA, ncol=ncol(drugpheno.gsk), nrow=nrow(sampleinfo.ge.gsk), dimnames=list(rownames(sampleinfo.ge.gsk), colnames(drugpheno.gsk)))
  matchix <- match(rownames(drugpheno.ge.gsk), gsub(badchars, "", rownames(drugpheno.gsk)))
  matchix.na <- is.na(matchix)
  drugpheno.ge.gsk[!matchix.na, ] <- drugpheno.gsk[matchix[!matchix.na], , drop=FALSE]
  colnames(drugpheno.ge.gsk) <- paste("drugid", toupper(sapply(strsplit(gsub(badchars, "", colnames(drugpheno.ge.gsk)), split="[(]gIC50"), function(x) { return(x[[1]]) })), sep="_")
  ## transofrm into microM
  drugpheno.ge.gsk <- drugpheno.ge.gsk / 1000
  ## drug concentrations
  iix <- which(!is.na(drugpheno.ge.gsk))
  posix <- pos2coord(pos=iix, dim.mat=dim(drugpheno.ge.gsk))
  drugconc.ge.gsk <- t(apply(posix, 1, function(x, y, z) {
    cc <- c(0.0003, 0.0032, 0.01, 0.032, 0.1, 0.317, 1, 3.16, 10)
    names(cc) <- sprintf("Dose%i.uM", 1:length(cc))
    return(c("cellid"=y[x[1]], "drugid"=z[x[2]], "nbr.conc.tested"=length(cc), cc))
  }, y=rownames(drugpheno.ge.gsk), z=colnames(drugpheno.ge.gsk)))
  rownames(drugconc.ge.gsk) <- paste(drugconc.ge.gsk[ , "cellid"], drugconc.ge.gsk[ , "drugid"], sep="...")
  
  message("Save GSK data")
  write.csv(annot.ge.gsk, file=file.path(saveres, "annot_ge_gsk.csv"), row.names=FALSE)
  write.csv(t(data.ge.gsk), file=file.path(saveres, "data_ge_gsk.csv"))
  write.csv(sampleinfo.ge.gsk, file=file.path(saveres, "sampleinfo_ge_gsk.csv"), row.names=FALSE)
  write.csv(drugpheno.ge.gsk, file=file.path(saveres, "drugpheno_ge_gsk.csv"))
  write.csv(drugconc.ge.gsk, file=file.path(saveres, "drugconc_ge_gsk.csv"))
  save(list=c("data.ge.gsk", "annot.ge.gsk", "sampleinfo.ge.gsk", "drugconc.ge.gsk", "drugpheno.ge.gsk", "mutation.gsk"), compress=TRUE, file=myfn)
}





## end




