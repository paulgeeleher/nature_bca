########################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## February 11, 2013
########################


###############################################################################
###############################################################################
## compute gene drug-associations using a linear model controlled for tissue types
###############################################################################
###############################################################################

require(amap) || stop("Library amap is not available!")
require(survcomp) || stop("Library survcomp is not available!")
require(vcd) || stop("Library vcd is not available!")
require(epibasix) || stop("Library gplots is not available!")
require(plotrix) || stop("Library plotrix is not available!")
require(WriteXLS) || stop("Library WriteXLS is not available!")
require(xtable) || stop("Library xtable is not available!")
require(gplots) || stop("Library gplots is not available!")

###############################################################################
###############################################################################
## load data
###############################################################################
###############################################################################

load(file.path(saveres, "CDRUG_cgp_ccle.RData"))

tissue.cgp <- as.character(sampleinfo.cgp[ , "tissue.type"])
tissue.ccle <- as.character(sampleinfo.ccle[ , "tissue.type"])
tissue <- tissue.cgp
tissuen <- sort(unique(as.character(tissue)))
ic50f.ccle <- ic50.ccle <- -log10(ic50.ccle / 10^6)
ic50f.cgp <- ic50.cgp <- -log10(ic50.cgp / 10^6)
## filter out low quality ic50 at high concentrations
ic50f.ccle[!ic50.filt.ccle] <- NA
ic50f.cgp[!ic50.filt.cgp] <- NA
drugn <- colnames(ic50.cgp)


########################
## IC50 all tissues
########################

## use only gene epxressions from CGP
myfn <- file.path(saveres, "assoc_genecgp_ic50_cgp_ccle.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on IC50 using gene expressions from CGP in CCLE")
  assoc.ic50.ccle <- NULL
  for(i in 1:ncol(ic50.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(ic50.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.cgp), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, ic50, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=ic50, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.cgp, ic50=ic50.ccle[ ,i], tissue=tissue.ccle)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.cgp), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.ic50.ccle <- c(assoc.ic50.ccle, list(mcres))
  }
  message("")
  names(assoc.ic50.ccle) <- colnames(ic50.ccle)
  ## CGP
  message("Gene-drug association based on IC50 using gene expressions from CGP in CGP")
  assoc.ic50.cgp <- NULL
  for(i in 1:ncol(ic50.cgp)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(ic50.cgp)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.cgp), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, ic50, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=ic50, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.cgp, ic50=ic50.cgp[ ,i], tissue=tissue.cgp)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.cgp), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.ic50.cgp <- c(assoc.ic50.cgp, list(mcres))
  }
  message("")
  names(assoc.ic50.cgp) <- colnames(ic50.cgp)
  save(list=c("assoc.ic50.cgp", "assoc.ic50.ccle"), compress=TRUE, file=myfn)
}

## use only gene epxressions from CCLE
myfn <- file.path(saveres, "assoc_geneccle_ic50_cgp_ccle.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on IC50 using gene expressions from CCLE in CCLE")
  assoc.ic50.ccle <- NULL
  for(i in 1:ncol(ic50.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(ic50.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, ic50, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=ic50, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.ccle, ic50=ic50.ccle[ ,i], tissue=tissue.ccle)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ccle), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.ic50.ccle <- c(assoc.ic50.ccle, list(mcres))
  }
  message("")
  names(assoc.ic50.ccle) <- colnames(ic50.ccle)
  ## CGP
  message("Gene-drug association based on IC50 using gene expressions from CCLE in CGP")
  assoc.ic50.cgp <- NULL
  for(i in 1:ncol(ic50.cgp)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(ic50.cgp)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, ic50, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=ic50, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.ccle, ic50=ic50.cgp[ ,i], tissue=tissue.cgp)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ccle), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.ic50.cgp <- c(assoc.ic50.cgp, list(mcres))
  }
  message("")
  names(assoc.ic50.cgp) <- colnames(ic50.cgp)
  save(list=c("assoc.ic50.cgp", "assoc.ic50.ccle"), compress=TRUE, file=myfn)
}

## use only ic50 from CGP
myfn <- file.path(saveres, "assoc_gene_ic50cgp_cgp_ccle.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on IC50 using IC50 from CGP in CCLE")
  assoc.ic50.ccle <- NULL
  for(i in 1:ncol(ic50.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(ic50.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, ic50, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=ic50, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.ccle, ic50=ic50.cgp[ ,i], tissue=tissue.ccle)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ccle), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.ic50.ccle <- c(assoc.ic50.ccle, list(mcres))
  }
  message("")
  names(assoc.ic50.ccle) <- colnames(ic50.ccle)
  ## CGP
  message("Gene-drug association based on IC50 using IC50 from CGP in CGP")
  assoc.ic50.cgp <- NULL
  for(i in 1:ncol(ic50.cgp)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(ic50.cgp)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.cgp), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, ic50, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=ic50, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.cgp, ic50=ic50.cgp[ ,i], tissue=tissue.cgp)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.cgp), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.ic50.cgp <- c(assoc.ic50.cgp, list(mcres))
  }
  message("")
  names(assoc.ic50.cgp) <- colnames(ic50.cgp)
  save(list=c("assoc.ic50.cgp", "assoc.ic50.ccle"), compress=TRUE, file=myfn)
}

## use only ic50 from CCLE
myfn <- file.path(saveres, "assoc_gene_ic50ccle_cgp_ccle.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on IC50 using IC50 from CCLE in CCLE")
  assoc.ic50.ccle <- NULL
  for(i in 1:ncol(ic50.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(ic50.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, ic50, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=ic50, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.ccle, ic50=ic50.ccle[ ,i], tissue=tissue.ccle)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ccle), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.ic50.ccle <- c(assoc.ic50.ccle, list(mcres))
  }
  message("")
  names(assoc.ic50.ccle) <- colnames(ic50.ccle)
  ## CGP
  message("Gene-drug association based on IC50 using IC50 from CCLE in CGP")
  assoc.ic50.cgp <- NULL
  for(i in 1:ncol(ic50.cgp)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(ic50.cgp)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.cgp), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, ic50, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=ic50, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.cgp, ic50=ic50.ccle[ ,i], tissue=tissue.cgp)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.cgp), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.ic50.cgp <- c(assoc.ic50.cgp, list(mcres))
  }
  message("")
  names(assoc.ic50.cgp) <- colnames(ic50.cgp)
  save(list=c("assoc.ic50.cgp", "assoc.ic50.ccle"), compress=TRUE, file=myfn)
}

tt <- matrix(NA, nrow=length(drugn), ncol=5, dimnames=list(drugn, c("rho", "lower", "upper", "p", "n")))
correlations.assoc.stats <- list("auc"=tt, "auc.filt"=tt, "auc.call"=tt, "auc.call.filt"=tt, "auc"=tt, "auc.filt"=tt, "auc.call"=tt, "auc.call.filt"=tt)

########################
## AUC all tissues
########################

## use only gene epxressions from CGP
myfn <- file.path(saveres, "assoc_genecgp_auc_cgp_ccle.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on AUC using gene expressions from CGP in CCLE")
  assoc.auc.ccle <- NULL
  for(i in 1:ncol(auc.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.cgp), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.cgp, auc=auc.ccle[ ,i], tissue=tissue.ccle)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.cgp), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.ccle <- c(assoc.auc.ccle, list(mcres))
  }
  message("")
  names(assoc.auc.ccle) <- colnames(auc.ccle)
  ## CGP
  message("Gene-drug association based on AUC using gene expressions from CGP in CGP")
  assoc.auc.cgp <- NULL
  for(i in 1:ncol(auc.cgp)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.cgp)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.cgp), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.cgp, auc=auc.cgp[ ,i], tissue=tissue.cgp)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.cgp), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.cgp <- c(assoc.auc.cgp, list(mcres))
  }
  message("")
  names(assoc.auc.cgp) <- colnames(auc.cgp)
  save(list=c("assoc.auc.cgp", "assoc.auc.ccle"), compress=TRUE, file=myfn)
}

## use only gene epxressions from CCLE
myfn <- file.path(saveres, "assoc_geneccle_auc_cgp_ccle.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on AUC using gene expressions from CCLE in CCLE")
  assoc.auc.ccle <- NULL
  for(i in 1:ncol(auc.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.ccle, auc=auc.ccle[ ,i], tissue=tissue.ccle)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ccle), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.ccle <- c(assoc.auc.ccle, list(mcres))
  }
  message("")
  names(assoc.auc.ccle) <- colnames(auc.ccle)
  ## CGP
  message("Gene-drug association based on AUC using gene expressions from CCLE in CGP")
  assoc.auc.cgp <- NULL
  for(i in 1:ncol(auc.cgp)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.cgp)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.ccle, auc=auc.cgp[ ,i], tissue=tissue.cgp)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ccle), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.cgp <- c(assoc.auc.cgp, list(mcres))
  }
  message("")
  names(assoc.auc.cgp) <- colnames(auc.cgp)
  save(list=c("assoc.auc.cgp", "assoc.auc.ccle"), compress=TRUE, file=myfn)
}

## use only auc from CGP
myfn <- file.path(saveres, "assoc_gene_auccgp_cgp_ccle.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on AUC using AUC from CGP in CCLE")
  assoc.auc.ccle <- NULL
  for(i in 1:ncol(auc.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.ccle, auc=auc.cgp[ ,i], tissue=tissue.ccle)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ccle), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.ccle <- c(assoc.auc.ccle, list(mcres))
  }
  message("")
  names(assoc.auc.ccle) <- colnames(auc.ccle)
  ## CGP
  message("Gene-drug association based on AUC using AUC from CGP in CGP")
  assoc.auc.cgp <- NULL
  for(i in 1:ncol(auc.cgp)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.cgp)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.cgp), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.cgp, auc=auc.cgp[ ,i], tissue=tissue.cgp)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.cgp), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.cgp <- c(assoc.auc.cgp, list(mcres))
  }
  message("")
  names(assoc.auc.cgp) <- colnames(auc.cgp)
  save(list=c("assoc.auc.cgp", "assoc.auc.ccle"), compress=TRUE, file=myfn)
}

## use only auc from CCLE
myfn <- file.path(saveres, "assoc_gene_aucccle_cgp_ccle.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on AUC using AUC from CCLE in CCLE")
  assoc.auc.ccle <- NULL
  for(i in 1:ncol(auc.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.ccle, auc=auc.ccle[ ,i], tissue=tissue.ccle)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ccle), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.ccle <- c(assoc.auc.ccle, list(mcres))
  }
  message("")
  names(assoc.auc.ccle) <- colnames(auc.ccle)
  ## CGP
  message("Gene-drug association based on AUC using AUC from CCLE in CGP")
  assoc.auc.cgp <- NULL
  for(i in 1:ncol(auc.cgp)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.cgp)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.cgp), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.cgp, auc=auc.ccle[ ,i], tissue=tissue.cgp)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.cgp), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.cgp <- c(assoc.auc.cgp, list(mcres))
  }
  message("")
  names(assoc.auc.cgp) <- colnames(auc.cgp)
  save(list=c("assoc.auc.cgp", "assoc.auc.ccle"), compress=TRUE, file=myfn)
}

myfn <- file.path(saveres, "correlations_assoc_ic50_auc_swap_cgp_ccle.RData")
if(!file.exists(myfn)) {
  ## gene-drug association using IC50
  filel <- list("Original"="assoc_ic50_cgp_ccle.RData", "GeneCGP.fixed"="assoc_genecgp_ic50_cgp_ccle.RData", "GeneCCLE.fixed"="assoc_geneccle_ic50_cgp_ccle.RData", "DrugCGP.fixed"="assoc_gene_ic50cgp_cgp_ccle.RData", "DrugCCLE.fixed"="assoc_gene_ic50ccle_cgp_ccle.RData")
  ## correlations for each drug in each setting
  ll <- ll.filt <- NULL
  for(i in 1:length(filel)) {
    ## load gene-drug associations
    load(file.path(saveres, filel[[i]]))
    ## all genes
    rr <- t(mapply(FUN=function(x, y) {
      myx <- complete.cases(x[ , "estimate"], y[ , "estimate"]) 
      if(sum(myx) < minsample) {
       cc <- list("estimate"=NA, "p.value"=NA) 
      } else {
        cc <- cor.test(x[myx, "estimate"], y[myx, "estimate"], method="spearman", use="complete.obs", alternative="greater")
      }
      return(c("rho"=cc$estimate, "p"=cc$p.value))
    }, x=assoc.ic50.cgp, y=assoc.ic50.ccle))
    colnames(rr) <- c("rho", "p")
    ll <- c(ll, list(rr))
    ## only genes significantly associated woith drug response
    rr <- t(mapply(FUN=function(x, y) {
      myx <- complete.cases(x[ , c("estimate", "fdr")], y[ , c("estimate", "fdr")]) & (x[ , "fdr"] < myfdr | x[ , "fdr"] < myfdr)
      if(sum(myx) < minsample) {
       cc <- list("estimate"=NA, "p.value"=NA) 
      } else {
        cc <- cor.test(x[myx, "estimate"], y[myx, "estimate"], method="spearman", use="complete.obs", alternative="greater")
      }
      return(c("rho"=cc$estimate, "p"=cc$p.value))
    }, x=assoc.ic50.cgp, y=assoc.ic50.ccle))
    colnames(rr) <- c("rho", "p")
    ll.filt <- c(ll.filt, list(rr))
  }
  names(ll) <- names(ll.filt) <- names(filel)
  ll.ic50 <- ll
  ll.filt.ic50 <- ll.filt

  ## gene-drug association using AUC
  filel <- list("Original"="assoc_auc_cgp_ccle.RData", "GeneCGP.fixed"="assoc_genecgp_auc_cgp_ccle.RData", "GeneCCLE.fixed"="assoc_geneccle_auc_cgp_ccle.RData", "DrugCGP.fixed"="assoc_gene_auccgp_cgp_ccle.RData", "DrugCCLE.fixed"="assoc_gene_aucccle_cgp_ccle.RData")
  ## correlations for each drug in each setting
  ll <- ll.filt <- NULL
  for(i in 1:length(filel)) {
    ## load gene-drug associations
    load(file.path(saveres, filel[[i]]))
    ## all genes
    rr <- t(mapply(FUN=function(x, y) {
      myx <- complete.cases(x[ , "estimate"], y[ , "estimate"]) 
      if(sum(myx) < minsample) {
       cc <- list("estimate"=NA, "p.value"=NA) 
      } else {
        cc <- cor.test(x[myx, "estimate"], y[myx, "estimate"], method="spearman", use="complete.obs", alternative="greater")
      }
      return(c("rho"=cc$estimate, "p"=cc$p.value))
    }, x=assoc.auc.cgp, y=assoc.auc.ccle))
    colnames(rr) <- c("rho", "p")
    ll <- c(ll, list(rr))
    ## only genes significantly associated woith drug response
    rr <- t(mapply(FUN=function(x, y) {
      myx <- complete.cases(x[ , c("estimate", "fdr")], y[ , c("estimate", "fdr")]) & (x[ , "fdr"] < myfdr | x[ , "fdr"] < myfdr)
      if(sum(myx) < minsample) {
       cc <- list("estimate"=NA, "p.value"=NA) 
      } else {
        cc <- cor.test(x[myx, "estimate"], y[myx, "estimate"], method="spearman", use="complete.obs", alternative="greater")
      }
      return(c("rho"=cc$estimate, "p"=cc$p.value))
    }, x=assoc.auc.cgp, y=assoc.auc.ccle))
    colnames(rr) <- c("rho", "p")
    ll.filt <- c(ll.filt, list(rr))
  }
  names(ll) <- names(ll.filt) <- names(filel)
  ll.auc <- ll
  ll.filt.auc <- ll.filt

  save(list=c("ll.ic50", "ll.filt.ic50", "ll.auc", "ll.filt.auc"), compress=TRUE, file=myfn)
} else { load(myfn) }

## boxplot for gene-drug associations based on ic50 in different settings
pdf(file.path(saveres, sprintf("boxplot2_assoc_gene_ic50_auc_cgp_ccle_swap_paper.pdf")), height=12, width=14)
par(mfrow=c(2, 2), las=1, mar=c(6, 4, 4, 2) + 0.1, xaxt="n")

## ic50, all genes
tt <- lapply(ll.ic50, function(x) { return(x[ , "rho"]) })
wt <- kruskal.test(x=tt[sapply(tt, length) > 0])
mp <- boxplot(tt, outline=FALSE, ylab="Rs", main="Correlations of gene-drug associations with IC50 (all)", border="grey50", ylim=c(-0.1, 1), col="white", range=0)
# # myccol <- gplots::rich.colors(length(tt), alpha=0.75)[-1]
myccol <- myBluePalette(length(tt))[-1]
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=c("lightgrey", myccol))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=1:length(tt))
points(x=jitter(unlist(xx), 1), y=unlist(tt), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", myccol))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")

## ic50, significant genes
# par(las=1, mar=c(6, 4, 4, 2) + 0.1, xaxt="n")
tt <- lapply(ll.filt.ic50, function(x) { return(x[ , "rho"]) })
wt <- kruskal.test(x=tt[sapply(tt, length) > 0])
mp <- boxplot(tt, outline=FALSE, ylab="Rs", main=sprintf("Correlations of gene-drug associations with IC50 (FDR < %i%%)", round(myfdr * 100)), border="grey50", ylim=c(-0.1, 1), col="white", range=0)
# myccol <- gplots::rich.colors(length(tt), alpha=0.75)[-1]
myccol <- myBluePalette(length(tt))[-1]
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=c("lightgrey", myccol))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=1:length(tt))
points(x=jitter(unlist(xx), 1), y=unlist(tt), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", myccol))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")

## auc, all genes
tt <- lapply(ll.auc, function(x) { return(x[ , "rho"]) })
wt <- kruskal.test(x=tt[sapply(tt, length) > 0])
mp <- boxplot(tt, outline=FALSE, ylab="Rs", main="Correlations of gene-drug associations with AUC (all)", border="grey50", ylim=c(-0.1, 1), col="white", range=0)
# myccol <- gplots::rich.colors(length(tt), alpha=0.75)[-1]
myccol <- myBluePalette(length(tt))[-1]
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=c("lightgrey", myccol))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=1:length(tt))
points(x=jitter(unlist(xx), 1), y=unlist(tt), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", myccol))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")

## auc, significant genes
# par(las=1, mar=c(6, 4, 4, 2) + 0.1, xaxt="n")
tt <- lapply(ll.filt.auc, function(x) { return(x[ , "rho"]) })
wt <- kruskal.test(x=tt[sapply(tt, length) > 0])
mp <- boxplot(tt, outline=FALSE, ylab="Rs", main=sprintf("Correlations of gene-drug associations with AUC (FDR < %i%%)", round(myfdr * 100)), border="grey50", ylim=c(-0.1, 1), col="white", range=0)
# myccol <- gplots::rich.colors(length(tt), alpha=0.75)[-1]
myccol <- myBluePalette(length(tt))[-1]
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=c("lightgrey", myccol))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=1:length(tt))
points(x=jitter(unlist(xx), 1), y=unlist(tt), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", myccol))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")

dev.off()


###############################################################################
###############################################################################
## compute pathway drug-associations using GSEA
###############################################################################
###############################################################################

########################
## IC50 all tissues
########################

myGSEAbis.ic50 <- function(myfn, myfn2) {
  if(!file.exists(myfn2)) {
    load(myfn)
    ## CCLE
    ## create rankings for each drug
    dir.create(file.path(saveres, "GSEAbis", "rankings", "CCLE", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
    rank.files <- NULL
    for(i in 1:length(assoc.ic50.ccle)) {
      ss <- assoc.ic50.ccle[[i]][ , "stat"]
      ss <- sort(ss, decreasing=TRUE, na.last=NA)
      ss[ss == Inf] <- .Machine$double.xmax
      ss[ss == -Inf] <- -.Machine$double.xmax
      rankg <- cbind(annot[names(ss), "jetset.EntrezID"], ss)
      fff <- file.path(saveres, "GSEAbis", "rankings", "CCLE", "ALLTISSUES", sprintf("ic50_%s.rnk", gsub("drugid_", "", names(assoc.ic50.ccle)[i])))
      write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
      rank.files <- c(rank.files, fff)  
    }
    names(rank.files) <- gsub("drugid_", "", names(assoc.ic50.ccle))
    ## GSEA
    dir.create(file.path(saveres, "GSEAbis", "reports", "CCLE", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
    gsea.res <- mclapply(1:length(rank.files), function(x, y, z, ...) { 
    	tt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=y[x], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=file.path(saveres, "GSEAbis", "reports", "CCLE", "ALLTISSUES"), replace.res=FALSE, gsea.seed=987654321)
      tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)], "description"=z[rownames(tt)])
      tt[tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (gsea.nperm + 1)
      return(tt)
    }, y=rank.files, z=genesets$description, gsea.exec=gsea.exec, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, saveres=saveres)
    names(gsea.res) <- names(rank.files)
    ## save all class comparisons into a multisheets xls file
    gsea.ic50.ccle <- lapply(gsea.res, function(x, y) {
      rr <- matrix(NA, nrow=length(y), ncol=3, dimnames=list(y, c("nes", "pvalue", "fdr")))
      rr[as.character(x[ , "geneset"]), ] <- data.matrix(x[ , c("NES", "NOM.p.val", "FDR.q.val")])
      return(rr)
    }, y=names(genesets$geneset))
    ## CGP
    ## create rankings for each drug
    dir.create(file.path(saveres, "GSEAbis", "rankings", "CGP", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
    rank.files <- NULL
    for(i in 1:length(assoc.ic50.cgp)) {
      ss <- assoc.ic50.cgp[[i]][ , "stat"]
      ss <- sort(ss, decreasing=TRUE, na.last=NA)
      ss[ss == Inf] <- .Machine$double.xmax
      ss[ss == -Inf] <- -.Machine$double.xmax
      rankg <- cbind(annot[names(ss), "jetset.EntrezID"], ss)
      fff <- file.path(saveres, "GSEAbis", "rankings", "CGP", "ALLTISSUES", sprintf("ic50_%s.rnk", gsub("drugid_", "", names(assoc.ic50.cgp)[i])))
      write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
      rank.files <- c(rank.files, fff)  
    }
    names(rank.files) <- gsub("drugid_", "", names(assoc.ic50.cgp))
    ## GSEA
    dir.create(file.path(saveres, "GSEAbis", "reports", "CGP", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
    gsea.res <- mclapply(1:length(rank.files), function(x, y, z, ...) { 
    	tt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=y[x], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=file.path(saveres, "GSEAbis", "reports", "CGP", "ALLTISSUES"), replace.res=FALSE, gsea.seed=987654321)
      tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)], "description"=z[rownames(tt)])
      tt[tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (gsea.nperm + 1)
      return(tt)
    }, y=rank.files, z=genesets$description, gsea.exec=gsea.exec, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, saveres=saveres)
    names(gsea.res) <- names(rank.files)
    ## save all class comparisons into a multisheets xls file
    gsea.ic50.cgp <- lapply(gsea.res, function(x, y) {
      rr <- matrix(NA, nrow=length(y), ncol=3, dimnames=list(y, c("nes", "pvalue", "fdr")))
      rr[as.character(x[ , "geneset"]), ] <- data.matrix(x[ , c("NES", "NOM.p.val", "FDR.q.val")])
      return(rr)
    }, y=names(genesets$geneset))
    unlink(file.path(saveres, "GSEAbis"), recursive=TRUE)
    save(list=c("gsea.ic50.ccle", "gsea.ic50.cgp"), compress=TRUE, file=myfn2)
  }
}

message("Pathway-drug association based on IC50 using gene expression from CGP")
myfn2 <- file.path(saveres, "gsea_genecgp_ic50_cgp_ccle.RData")
myfn <- file.path(saveres, "assoc_genecgp_ic50_cgp_ccle.RData")
myGSEAbis.ic50(myfn=myfn, myfn2=myfn2)

message("Pathway-drug association based on IC50 using gene expression from CCLE")
myfn2 <- file.path(saveres, "gsea_geneccle_ic50_cgp_ccle.RData")
myfn <- file.path(saveres, "assoc_geneccle_ic50_cgp_ccle.RData")
myGSEAbis.ic50(myfn=myfn, myfn2=myfn2)

message("Pathway-drug association based on IC50 using IC50 from CGP")
myfn2 <- file.path(saveres, "gsea_gene_ic50cgp_cgp_ccle.RData")
myfn <- file.path(saveres, "assoc_gene_ic50cgp_cgp_ccle.RData")
myGSEAbis.ic50(myfn=myfn, myfn2=myfn2)

message("Pathway-drug association based on IC50 using IC50 from CCLE")
myfn2 <- file.path(saveres, "gsea_gene_ic50ccle_cgp_ccle.RData")
myfn <- file.path(saveres, "assoc_gene_ic50ccle_cgp_ccle.RData")
myGSEAbis.ic50(myfn=myfn, myfn2=myfn2)


########################
## AUC all tissues
########################

myGSEAbis.auc <- function(myfn, myfn2) {
  if(!file.exists(myfn2)) {
    load(myfn)
    ## CCLE
    ## create rankings for each drug
    dir.create(file.path(saveres, "GSEAbis", "rankings", "CCLE", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
    rank.files <- NULL
    for(i in 1:length(assoc.auc.ccle)) {
      ss <- assoc.auc.ccle[[i]][ , "stat"]
      ss <- sort(ss, decreasing=TRUE, na.last=NA)
      ss[ss == Inf] <- .Machine$double.xmax
      ss[ss == -Inf] <- -.Machine$double.xmax
      rankg <- cbind(annot[names(ss), "jetset.EntrezID"], ss)
      fff <- file.path(saveres, "GSEAbis", "rankings", "CCLE", "ALLTISSUES", sprintf("auc_%s.rnk", gsub("drugid_", "", names(assoc.auc.ccle)[i])))
      write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
      rank.files <- c(rank.files, fff)  
    }
    names(rank.files) <- gsub("drugid_", "", names(assoc.auc.ccle))
    ## GSEA
    dir.create(file.path(saveres, "GSEAbis", "reports", "CCLE", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
    gsea.res <- mclapply(1:length(rank.files), function(x, y, z, ...) { 
    	tt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=y[x], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=file.path(saveres, "GSEAbis", "reports", "CCLE", "ALLTISSUES"), replace.res=FALSE, gsea.seed=987654321)
      tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)], "description"=z[rownames(tt)])
      tt[tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (gsea.nperm + 1)
      return(tt)
    }, y=rank.files, z=genesets$description, gsea.exec=gsea.exec, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, saveres=saveres)
    names(gsea.res) <- names(rank.files)
    ## save all class comparisons into a multisheets xls file
    gsea.auc.ccle <- lapply(gsea.res, function(x, y) {
      rr <- matrix(NA, nrow=length(y), ncol=3, dimnames=list(y, c("nes", "pvalue", "fdr")))
      rr[as.character(x[ , "geneset"]), ] <- data.matrix(x[ , c("NES", "NOM.p.val", "FDR.q.val")])
      return(rr)
    }, y=names(genesets$geneset))
    ## CGP
    ## create rankings for each drug
    dir.create(file.path(saveres, "GSEAbis", "rankings", "CGP", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
    rank.files <- NULL
    for(i in 1:length(assoc.auc.cgp)) {
      ss <- assoc.auc.cgp[[i]][ , "stat"]
      ss <- sort(ss, decreasing=TRUE, na.last=NA)
      ss[ss == Inf] <- .Machine$double.xmax
      ss[ss == -Inf] <- -.Machine$double.xmax
      rankg <- cbind(annot[names(ss), "jetset.EntrezID"], ss)
      fff <- file.path(saveres, "GSEAbis", "rankings", "CGP", "ALLTISSUES", sprintf("auc_%s.rnk", gsub("drugid_", "", names(assoc.auc.cgp)[i])))
      write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
      rank.files <- c(rank.files, fff)  
    }
    names(rank.files) <- gsub("drugid_", "", names(assoc.auc.cgp))
    ## GSEA
    dir.create(file.path(saveres, "GSEAbis", "reports", "CGP", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
    gsea.res <- mclapply(1:length(rank.files), function(x, y, z, ...) { 
    	tt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=y[x], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=file.path(saveres, "GSEAbis", "reports", "CGP", "ALLTISSUES"), replace.res=FALSE, gsea.seed=987654321)
      tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)], "description"=z[rownames(tt)])
      tt[tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (gsea.nperm + 1)
      return(tt)
    }, y=rank.files, z=genesets$description, gsea.exec=gsea.exec, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, saveres=saveres)
    names(gsea.res) <- names(rank.files)
    ## save all class comparisons into a multisheets xls file
    gsea.auc.cgp <- lapply(gsea.res, function(x, y) {
      rr <- matrix(NA, nrow=length(y), ncol=3, dimnames=list(y, c("nes", "pvalue", "fdr")))
      rr[as.character(x[ , "geneset"]), ] <- data.matrix(x[ , c("NES", "NOM.p.val", "FDR.q.val")])
      return(rr)
    }, y=names(genesets$geneset))
    unlink(file.path(saveres, "GSEAbis"), recursive=TRUE)
    save(list=c("gsea.auc.ccle", "gsea.auc.cgp"), compress=TRUE, file=myfn2)
  }
}

message("Pathway-drug association based on AUC using gene expression from CGP")
myfn2 <- file.path(saveres, "gsea_genecgp_auc_cgp_ccle.RData")
myfn <- file.path(saveres, "assoc_genecgp_auc_cgp_ccle.RData")
myGSEAbis.auc(myfn=myfn, myfn2=myfn2)

message("Pathway-drug association based on AUC using gene expression from CCLE")
myfn2 <- file.path(saveres, "gsea_geneccle_auc_cgp_ccle.RData")
myfn <- file.path(saveres, "assoc_geneccle_auc_cgp_ccle.RData")
myGSEAbis.auc(myfn=myfn, myfn2=myfn2)

message("Pathway-drug association based on AUC using AUC from CGP")
myfn2 <- file.path(saveres, "gsea_gene_auccgp_cgp_ccle.RData")
myfn <- file.path(saveres, "assoc_gene_auccgp_cgp_ccle.RData")
myGSEAbis.auc(myfn=myfn, myfn2=myfn2)

message("Pathway-drug association based on AUC using AUC from CCLE")
myfn2 <- file.path(saveres, "gsea_gene_aucccle_cgp_ccle.RData")
myfn <- file.path(saveres, "assoc_gene_aucccle_cgp_ccle.RData")
myGSEAbis.auc(myfn=myfn, myfn2=myfn2)

myfn <- file.path(saveres, "correlations_gsea_ic50_auc_swap_cgp_ccle.RData")
if(!file.exists(myfn)) {
  ## pathway-drug association using IC50
  filel <- list("Original"="gsea_ic50_cgp_ccle.RData", "GeneCGP.fixed"="gsea_genecgp_ic50_cgp_ccle.RData", "GeneCCLE.fixed"="gsea_geneccle_ic50_cgp_ccle.RData", "DrugCGP.fixed"="gsea_gene_ic50cgp_cgp_ccle.RData", "DrugCCLE.fixed"="gsea_gene_ic50ccle_cgp_ccle.RData")
  ## correlations for each drug in each setting
  ll <- ll.filt <- NULL
  for(i in 1:length(filel)) {
    ## load gene-drug associations
    load(file.path(saveres, filel[[i]]))
    ## all genes
    rr <- t(mapply(FUN=function(x, y) {
      myx <- complete.cases(x[ , "nes"], y[ , "nes"]) 
      if(sum(myx) < minsample) {
       cc <- list("estimate"=NA, "p.value"=NA) 
      } else {
        cc <- cor.test(x[myx, "nes"], y[myx, "nes"], method="spearman", use="complete.obs", alternative="greater")
      }
      return(c("rho"=cc$estimate, "p"=cc$p.value))
    }, x=gsea.ic50.cgp, y=gsea.ic50.ccle))
    colnames(rr) <- c("rho", "p")
    ll <- c(ll, list(rr))
    ## only genes significantly associated woith drug response
    rr <- t(mapply(FUN=function(x, y) {
      myx <- complete.cases(x[ , c("nes", "fdr")], y[ , c("nes", "fdr")]) & (x[ , "fdr"] < myfdr | x[ , "fdr"] < myfdr)
      if(sum(myx) < minsample) {
       cc <- list("estimate"=NA, "p.value"=NA) 
      } else {
        cc <- cor.test(x[myx, "nes"], y[myx, "nes"], method="spearman", use="complete.obs", alternative="greater")
      }
      return(c("rho"=cc$estimate, "p"=cc$p.value))
    }, x=gsea.ic50.cgp, y=gsea.ic50.ccle))
    colnames(rr) <- c("rho", "p")
    ll.filt <- c(ll.filt, list(rr))
  }
  names(ll) <- names(ll.filt) <- names(filel)
  ll.ic50 <- ll
  ll.filt.ic50 <- ll.filt

  ## pathway-drug association using AUC
  filel <- list("Original"="gsea_auc_cgp_ccle.RData", "GeneCGP.fixed"="gsea_genecgp_auc_cgp_ccle.RData", "GeneCCLE.fixed"="gsea_geneccle_auc_cgp_ccle.RData", "DrugCGP.fixed"="gsea_gene_auccgp_cgp_ccle.RData", "DrugCCLE.fixed"="gsea_gene_aucccle_cgp_ccle.RData")
  ## correlations for each drug in each setting
  ll <- ll.filt <- NULL
  for(i in 1:length(filel)) {
    ## load gene-drug associations
    load(file.path(saveres, filel[[i]]))
    ## all genes
    rr <- t(mapply(FUN=function(x, y) {
      myx <- complete.cases(x[ , "nes"], y[ , "nes"]) 
      if(sum(myx) < minsample) {
       cc <- list("estimate"=NA, "p.value"=NA) 
      } else {
        cc <- cor.test(x[myx, "nes"], y[myx, "nes"], method="spearman", use="complete.obs", alternative="greater")
      }
      return(c("rho"=cc$estimate, "p"=cc$p.value))
    }, x=gsea.auc.cgp, y=gsea.auc.ccle))
    colnames(rr) <- c("rho", "p")
    ll <- c(ll, list(rr))
    ## only genes significantly associated woith drug response
    rr <- t(mapply(FUN=function(x, y) {
      myx <- complete.cases(x[ , c("nes", "fdr")], y[ , c("nes", "fdr")]) & (x[ , "fdr"] < myfdr | x[ , "fdr"] < myfdr)
      if(sum(myx) < minsample) {
       cc <- list("estimate"=NA, "p.value"=NA) 
      } else {
        cc <- cor.test(x[myx, "nes"], y[myx, "nes"], method="spearman", use="complete.obs", alternative="greater")
      }
      return(c("rho"=cc$estimate, "p"=cc$p.value))
    }, x=gsea.auc.cgp, y=gsea.auc.ccle))
    colnames(rr) <- c("rho", "p")
    ll.filt <- c(ll.filt, list(rr))
  }
  names(ll) <- names(ll.filt) <- names(filel)
  ll.auc <- ll
  ll.filt.auc <- ll.filt

  save(list=c("ll.ic50", "ll.filt.ic50", "ll.auc", "ll.filt.auc"), compress=TRUE, file=myfn)
} else { load(myfn) }

## boxplot for gene-drug associations based on ic50 in different settings
pdf(file.path(saveres, sprintf("boxplot2_gsea_gene_ic50_auc_cgp_ccle_swap_paper.pdf")), height=12, width=14)
par(mfrow=c(2, 2), las=1, mar=c(6, 4, 4, 2) + 0.1, xaxt="n")

## ic50, all genes
tt <- lapply(ll.ic50, function(x) { return(x[ , "rho"]) })
wt <- kruskal.test(x=tt[sapply(tt, length) > 0])
mp <- boxplot(tt, outline=FALSE, ylab="Rs", main="Correlations of pathway-drug associations with IC50 (all)", border="grey50", ylim=c(-0.3, 1), col="white", range=0)
# myccol <- gplots::rich.colors(length(tt), alpha=0.75)[-1]
myccol <- myBluePalette(length(tt))[-1]
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=c("lightgrey", myccol))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=1:length(tt))
points(x=jitter(unlist(xx), 1), y=unlist(tt), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", myccol))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")

## ic50, significant genes
# par(las=1, mar=c(6, 4, 4, 2) + 0.1, xaxt="n")
tt <- lapply(ll.filt.ic50, function(x) { return(x[ , "rho"]) })
wt <- kruskal.test(x=tt[sapply(tt, length) > 0])
mp <- boxplot(tt, outline=FALSE, ylab="Rs", main=sprintf("Correlations of pathway-drug associations with IC50 (FDR < %i%%)", round(myfdr * 100)), border="grey50", ylim=c(-0.3, 1), col="white", range=0)
# myccol <- gplots::rich.colors(length(tt), alpha=0.75)[-1]
myccol <- myBluePalette(length(tt))[-1]
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=c("lightgrey", myccol))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=1:length(tt))
points(x=jitter(unlist(xx), 1), y=unlist(tt), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", myccol))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")

## auc, all genes
tt <- lapply(ll.auc, function(x) { return(x[ , "rho"]) })
wt <- kruskal.test(x=tt[sapply(tt, length) > 0])
mp <- boxplot(tt, outline=FALSE, ylab="Rs", main="Correlations of pathway-drug associations with AUC (all)", border="grey50", ylim=c(-0.3, 1), col="white", range=0)
# myccol <- gplots::rich.colors(length(tt), alpha=0.75)[-1]
myccol <- myBluePalette(length(tt))[-1]
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=c("lightgrey", myccol))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=1:length(tt))
points(x=jitter(unlist(xx), 1), y=unlist(tt), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", myccol))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")

## auc, significant genes
# par(las=1, mar=c(6, 4, 4, 2) + 0.1, xaxt="n")
tt <- lapply(ll.filt.auc, function(x) { return(x[ , "rho"]) })
wt <- kruskal.test(x=tt[sapply(tt, length) > 0])
mp <- boxplot(tt, outline=FALSE, ylab="Rs", main=sprintf("Correlations of pathway-drug associations with AUC (FDR < %i%%)", round(myfdr * 100)), border="grey50", ylim=c(-0.3, 1), col="white", range=0)
# myccol <- gplots::rich.colors(length(tt), alpha=0.75)[-1]
myccol <- myBluePalette(length(tt))[-1]
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=c("lightgrey", myccol))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=1:length(tt))
points(x=jitter(unlist(xx), 1), y=unlist(tt), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", myccol))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")

dev.off()

###############################################################################
###############################################################################
## only results with all genes
###############################################################################
###############################################################################

## boxplot for gene-drug associations based on ic50 in different settings
pdf(file.path(saveres, sprintf("boxplot2_gsea_assoc_gene_ic50_auc_cgp_ccle_swapr.pdf")), height=12, width=14)
par(mfrow=c(2, 2), las=1, mar=c(6, 4, 4, 2) + 0.1, xaxt="n")

yylim <- c(-0.1, 1.1)
## gene-drug associations
load(file.path(saveres, "correlations_assoc_ic50_auc_swap_cgp_ccle.RData"))
## IC50
tt <- lapply(ll.ic50, function(x) { return(x[ , "rho"]) })
wt <- kruskal.test(x=tt[sapply(tt, length) > 0])
mp <- boxplot(tt, outline=FALSE, ylab="Rs", main="Correlations of gene-drug associations with IC50", border="grey50", ylim=yylim, col="white", range=0)
# myccol <- gplots::rich.colors(length(tt), alpha=0.75)[-1]
myccol <- myBluePalette(length(tt))[-1]
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=c("lightgrey", myccol))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=1:length(tt))
points(x=jitter(unlist(xx), 1), y=unlist(tt), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", myccol))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
## AUC
tt <- lapply(ll.auc, function(x) { return(x[ , "rho"]) })
wt <- kruskal.test(x=tt[sapply(tt, length) > 0])
mp <- boxplot(tt, outline=FALSE, ylab="Rs", main=sprintf("Correlations of gene-drug associations with AUC", round(myfdr * 100)), border="grey50", ylim=yylim, col="white", range=0)
# myccol <- gplots::rich.colors(length(tt), alpha=0.75)[-1]
myccol <- myBluePalette(length(tt))[-1]
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=c("lightgrey", myccol))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=1:length(tt))
points(x=jitter(unlist(xx), 1), y=unlist(tt), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", myccol))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")

yylim <- c(-0.3, 1.1)
## pathway-drug assocations
load(file.path(saveres, "correlations_gsea_ic50_auc_swap_cgp_ccle.RData"))
## auc, all genes
tt <- lapply(ll.ic50, function(x) { return(x[ , "rho"]) })
wt <- kruskal.test(x=tt[sapply(tt, length) > 0])
mp <- boxplot(tt, outline=FALSE, ylab="Rs", main="Correlations of pathway-drug associations with IC50", border="grey50", ylim=yylim, col="white", range=0)
# myccol <- gplots::rich.colors(length(tt), alpha=0.75)[-1]
myccol <- myBluePalette(length(tt))[-1]
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=c("lightgrey", myccol))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=1:length(tt))
points(x=jitter(unlist(xx), 1), y=unlist(tt), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", myccol))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")

## auc, significant genes
# par(las=1, mar=c(6, 4, 4, 2) + 0.1, xaxt="n")
tt <- lapply(ll.auc, function(x) { return(x[ , "rho"]) })
wt <- kruskal.test(x=tt[sapply(tt, length) > 0])
mp <- boxplot(tt, outline=FALSE, ylab="Rs", main=sprintf("Correlations of pathway-drug associations with AUC", round(myfdr * 100)), border="grey50", ylim=yylim, col="white", range=0)
# myccol <- gplots::rich.colors(length(tt), alpha=0.75)[-1]
myccol <- myBluePalette(length(tt))[-1]
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=c("lightgrey", myccol))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(tt, length), y=1:length(tt))
points(x=jitter(unlist(xx), 1), y=unlist(tt), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", myccol))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")

dev.off()




## end



