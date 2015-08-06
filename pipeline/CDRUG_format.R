########################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## June 13, 2013
########################

# rm(list=ls())

require(vcd) || stop("Library vcd is not available")
require(epibasix) || stop("Library epibasix is not available")
require(R.utils) || stop("Library R.utils is not available")
require(amap) || stop("Library amap is not available")
require(gplots) || stop("Library gplots is not available")
require(Vennerable) || stop("Library Vennerable is not available")

## CGP
load(file.path(saveres, "cgp_frma.RData"))
## gene centric data
myx <- which(annot.ge.cgp[ ,"best"])
myx <- myx[!duplicated(annot.ge.cgp[myx,"jetset.EntrezID"])]
data.cgp <- data.ge.cgp[ ,myx,drop=FALSE]
annot.cgp <- annot.ge.cgp[myx, ,drop=FALSE]
colnames(data.cgp) <- rownames(annot.cgp) <- paste("geneid", annot.cgp[ ,"jetset.EntrezID"], sep="_")
## drug information
druginfo.cgp <- druginfo.ge.cgp
## sample information
sampleinfo.cgp <- sampleinfo.ge.cgp
## drug sensitivity
drugpheno.cgp <- drugpheno.ge.cgp
## drug concentrations
drugconc.cgp <- drugconc.ge.cgp

## CCLE
load(file.path(saveres, "ccle_frma.RData"))
## gene centric data
myx <- which(annot.ge.ccle[ ,"best"])
myx <- myx[!duplicated(annot.ge.ccle[myx,"jetset.EntrezID"])]
data.ccle <- data.ge.ccle[ ,myx,drop=FALSE]
annot.ccle <- annot.ge.ccle[myx, ,drop=FALSE]
colnames(data.ccle) <- rownames(annot.ccle) <- paste("geneid", annot.ccle[ ,"jetset.EntrezID"], sep="_")
## drug information
druginfo.ccle <- druginfo.ge.ccle
## sample information
sampleinfo.ccle <- sampleinfo.ge.ccle
## drug sensitivity
drugpheno.ccle <- drugpheno.ge.ccle
## drug concentrations
drugconc.ccle <- drugconc.ge.ccle

## gsk
load(file.path(saveres, "gskcellines_frma.RData"))
## gene centric data
myx <- which(annot.ge.gsk[ ,"best"])
myx <- myx[!duplicated(annot.ge.gsk[myx,"jetset.EntrezID"])]
data.gsk <- data.ge.gsk[ ,myx,drop=FALSE]
annot.gsk <- annot.ge.gsk[myx, ,drop=FALSE]
colnames(data.gsk) <- rownames(annot.gsk) <- paste("geneid", annot.gsk[ ,"jetset.EntrezID"], sep="_")
## sample information
sampleinfo.gsk <- sampleinfo.ge.gsk
## drug sensitivity
drugpheno.gsk <- drugpheno.ge.gsk
## drug concentrations
drugconc.gsk <- drugconc.ge.gsk


## cell line annotations
## CGP cell line collection
# celline.cgp <- read.csv(file=file.path(saveres, "cell_line_collection_cgp.csv"), stringsAsFactors=FALSE)
# rownames(celline.cgp) <- as.character(celline.cgp[ , "cellid"])
# celline.cgp[ , "cellid"] <- as.character(celline.cgp[ , "cellid"])
if(any(!is.element(rownames(data.cgp), rownames(celline.cgp)))) { stop("some cell lines in CGP are not part of CGP cell line collection") }
## CCLE cell line collection
# celline.ccle <- read.csv(file=file.path(saveres, "cell_line_collection_ccle.csv"), stringsAsFactors=FALSE)
# rownames(celline.ccle) <- as.character(celline.ccle[ , "cellid"])
# celline.ccle[ , "cellid"] <- as.character(celline.ccle[ , "cellid"])
if(any(!is.element(rownames(data.ccle), rownames(celline.ccle)))) { stop("some cell lines in CCLE are not part of CCLE cell line collection") }

## read manual matching for cell lines in CCLE and CGP
match.ccle.cgp <- read.csv(file=file.path("matching_cell_line_CCLE_CGP.csv"), stringsAsFactors=FALSE)

## use CGP cell lines as reference 
## update ccle cell line collection
nn <- intersect(as.character(match.ccle.cgp[ , "CCLE.cell.line"]), as.character(celline.ccle[ , "cellid"]))
iix0 <- which(is.element(as.character(match.ccle.cgp[ , "CCLE.cell.line"]), nn))
iix <- match(as.character(match.ccle.cgp[iix0, "CCLE.cell.line"]), as.character(celline.ccle[ , "cellid"]))
celline.ccle[iix, "cellid"] <- as.character(match.ccle.cgp[iix0, "CGP.cell.line"])
rownames(celline.ccle) <- as.character(celline.ccle[ , "cellid"])
## update ccle gene expression
nn <- intersect(as.character(match.ccle.cgp[ , "CCLE.cell.line"]), rownames(data.ccle))
iix0 <- which(is.element(as.character(match.ccle.cgp[ , "CCLE.cell.line"]), nn))
iix <- match(as.character(match.ccle.cgp[iix0, "CCLE.cell.line"]), rownames(data.ccle))
rownames(data.ccle)[iix] <- as.character(match.ccle.cgp[iix0, "CGP.cell.line"])
## update ccle sample information
nn <- intersect(as.character(match.ccle.cgp[ , "CCLE.cell.line"]), as.character(sampleinfo.ccle[ , "cellid"]))
iix0 <- which(is.element(as.character(match.ccle.cgp[ , "CCLE.cell.line"]), nn))
iix <- match(as.character(match.ccle.cgp[iix0, "CCLE.cell.line"]), as.character(sampleinfo.ccle[ , "cellid"]))
sampleinfo.ccle[iix, "cellid"] <- as.character(match.ccle.cgp[iix0, "CGP.cell.line"])
rownames(sampleinfo.ccle) <- as.character(sampleinfo.ccle[ , "cellid"])
## update ccle drug phenotypes
nn <- intersect(as.character(match.ccle.cgp[ , "CCLE.cell.line"]), rownames(drugpheno.ccle))
iix0 <- which(is.element(as.character(match.ccle.cgp[ , "CCLE.cell.line"]), nn))
iix <- match(as.character(match.ccle.cgp[iix0, "CCLE.cell.line"]), rownames(drugpheno.ccle))
rownames(drugpheno.ccle)[iix] <- as.character(match.ccle.cgp[iix0, "CGP.cell.line"])
## update ccle drug concentrations
nn <- intersect(as.character(match.ccle.cgp[ , "CCLE.cell.line"]), drugconc.ccle[ , "cellid"])
iix0 <- which(is.element(as.character(match.ccle.cgp[ , "CCLE.cell.line"]), nn))
for(i in 1:length(iix0)) {
  myx <- !is.na(drugconc.ccle[ , "cellid"]) & drugconc.ccle[ , "cellid"] == match.ccle.cgp[iix0[i], "CCLE.cell.line"]
  drugconc.ccle[myx, "cellid"] <- as.character(match.ccle.cgp[iix0[i], "CGP.cell.line"])
}
rownames(drugconc.ccle) <- paste(as.character(drugconc.ccle[ , "cellid"]), as.character(drugconc.ccle[ , "drugid"]), sep="...")
## update ccle mutations
nn <- intersect(as.character(match.ccle.cgp[ , "CCLE.cell.line"]), rownames(mutation.ccle))
iix0 <- which(is.element(as.character(match.ccle.cgp[ , "CCLE.cell.line"]), nn))
iix <- match(as.character(match.ccle.cgp[iix0, "CCLE.cell.line"]), rownames(mutation.ccle))
rownames(mutation.ccle)[iix] <- as.character(match.ccle.cgp[iix0, "CGP.cell.line"])
## update ccle cell line collection
nn <- intersect(as.character(match.ccle.cgp[ , "CCLE.cell.line"]), as.character(drugconc.ccle[ , "cellid"]))
iix0 <- which(is.element(as.character(match.ccle.cgp[ , "CCLE.cell.line"]), nn))

## GSK
# ## matching with CGP and CCLE
# ## best matching based on cell line name
# tt1 <- lapply(rownames(data.gsk), agrep, x=rownames(data.cgp), value=TRUE)
# tt1 <- t(mapply(function(x, y) { 
#   if(length(y) >= 1) {
#     myx <- is.element(y, x)
#     if(any(myx)) { yy <- y[myx] } else { yy <- y[1] }
#     return(cbind(x, yy))
#   } else { return(cbind(x, NA)) }
# }, x=as.list(rownames(data.gsk)), y=tt1))
# tt2 <- lapply(rownames(data.gsk), agrep, x=rownames(data.ccle), value=TRUE)
# tt2 <- t(mapply(function(x, y) { 
#   if(length(y) >= 1) {
#     myx <- is.element(y, x)
#     if(any(myx)) { yy <- y[myx] } else { yy <- y[1] }
#     return(cbind(x, yy))
#   } else { return(cbind(x, NA)) }
# }, x=as.list(rownames(data.gsk)), y=tt2))
# match1 <- cbind(tt1, tt2[ , 2])
# colnames(match1) <- c("GSK", "CGP", "CCLE")
# ## best matching based on correlation of gene epxression profiles
# genen <- fold(intersect, colnames(data.gsk), colnames(data.cgp), colnames(data.ccle))
# genen <- genen[order(apply(data.gsk[ , genen, drop=FALSE], 2, var, na.rm=TRUE), decreasing=TRUE)[1:topvar]]
# rr <- cor(x=t(data.gsk[ , genen, drop=FALSE]), y=t(data.cgp[ , genen, drop=FALSE]))
# match2 <- cbind("GSK"=rownames(data.gsk), "CGP"=sapply(apply(rr, 1, which.max), function(x, y) { if(length(x) > 0) { return(y[x]) } else { return(NA)} }, y=rownames(data.cgp)))
# rr <- cor(x=t(data.gsk[ , genen, drop=FALSE]), y=t(data.ccle[ , genen, drop=FALSE]))
# match2 <- cbind(match2, "CCLE"=sapply(apply(rr, 1, which.max), function(x, y) { if(length(x) > 0) { return(y[x]) } else { return(NA)} }, y=rownames(data.ccle)))
# write.csv(cbind(match1, "", match2), "temp.csv", row.names=FALSE)
## read the (manually curated) matching file
mm <- read.csv("matching_cell_line_GSK_CCLE_CGP.csv", stringsAsFactors=FALSE)
mm[mm == ""] <- NA
rownames(mm) <- mm[ , 1]
nnn <- nnn0 <- rownames(sampleinfo.gsk)
names(nnn) <- nnn
nnn[rownames(mm)[!is.na(mm[ , "CGP"])]] <- mm[!is.na(mm[ , "CGP"]), "CGP"]
nnn[rownames(mm)[is.na(mm[ , "CGP"])]] <- mm[is.na(mm[ , "CGP"]), "CCLE"]
## update cell line names
rownames(data.gsk) <- rownames(sampleinfo.gsk) <- sampleinfo.gsk[ , "cellid"] <- rownames(drugpheno.gsk) <- nnn
## update ccle drug concentrations
nn <- intersect(as.character(match.ccle.cgp[ , "CCLE.cell.line"]), drugconc.ccle[ , "cellid"])
iix0 <- which(nnn != names(nnn))
for(i in 1:length(iix0)) {
  myx <- !is.na(drugconc.gsk[ , "cellid"]) & drugconc.gsk[ , "cellid"] == names(nnn)[iix0[i]]
  drugconc.gsk[myx, "cellid"] <- as.character(nnn[iix0[i]])
}
rownames(drugconc.gsk) <- paste(as.character(drugconc.gsk[ , "cellid"]), as.character(drugconc.gsk[ , "drugid"]), sep="...")

## curate tissue types
# cc <- intersect(rownames(sampleinfo.gsk), rownames(sampleinfo.cgp))
# tt <- as.matrix(table(as.character(sampleinfo.gsk[cc, "tissue.type"]), as.character(sampleinfo.cgp[cc, "tissue.type"])))
# match1 <- colnames(tt)[apply(tt, 1, which.max)]
# names(match1) <- rownames(tt)
# cc <- intersect(rownames(sampleinfo.gsk), rownames(sampleinfo.ccle))
# tt <- as.matrix(table(as.character(sampleinfo.gsk[cc, "tissue.type"]), as.character(sampleinfo.ccle[cc, "tissue.type"])))
# match2 <- colnames(tt)[apply(tt, 1, which.max)]
# names(match2) <- rownames(tt)
# nn <- union(names(match1), names(match2))
# nn <- matrix(NA, nrow=length(nn), ncol=3, dimnames=list(nn, c("GSK", "CGP", "CCLE")))
# nn[ , "GSK"] <- rownames(nn)
# nn[names(match1), "CGP"] <- match1
# nn[names(match2), "CCLE"] <- match2
# write.csv(nn, row.names=FALSE, file="temp2.csv")
mm <- read.csv("matching_tissue_type_GSK_CCLE_CGP.csv", stringsAsFactors=FALSE)
rownames(mm) <- mm[ , "GSK"]
tissue <- sampleinfo.gsk[ , "tissue.type"]
levels(tissue) <- mm[levels(tissue), "CGP"]
sampleinfo.gsk[ , "tissue.type"] <- as.character(tissue)

## intersection between GSK, CGP and CCLE
ll <- list("CGP"=sampleinfo.cgp[ , "cellid"], "CCLE"=sampleinfo.ccle[ , "cellid"], "GSK"=sampleinfo.gsk[ , "cellid"])
pdf(file.path(file.path(saveres, "intersection_cellines_cgp_ccle_gsk_paper.pdf")), width=5, height=5)
ww <- Vennerable::Venn(Sets=ll)
plot(ww, show=list(Universe=FALSE))
dev.off()
pdf(file.path(file.path(saveres, "intersection_cellines_cgp_ccle_paper.pdf")), width=5, height=5)
ww <- Vennerable::Venn(Sets=ll[c("CGP", "CCLE")])
plot(ww, show=list(Universe=FALSE))
dev.off()

## merge cell line annotations from CCLE and CGP for which genomic data are available
celln <- sort(unique(c(rownames(data.ccle), rownames(data.cgp))))
celline.collection <- data.frame(matrix(NA, nrow=length(celln), ncol=7, dimnames=list(celln, c("cellid", "CGP.cell.line", "CCLE.cell.line", "CGP.tissue.type", "CCLE.tissue.type", "CGP.link", "CCLE.link"))))
celline.collection[ , "cellid"] <- celln
celline.collection[rownames(data.cgp), c("CGP.cell.line", "CGP.tissue.type", "CGP.link")] <- celline.cgp[rownames(data.cgp), c("Sample.name", "Primary.site", "link")]
celline.collection[rownames(data.ccle), c("CCLE.cell.line", "CCLE.tissue.type", "CCLE.link")] <- celline.ccle[rownames(data.ccle), c("Cell.line.primary.name", "Site.Primary", "link")]

## merge tissue types
iix <- apply(celline.collection[ , c("CGP.tissue.type", "CCLE.tissue.type")], 1, function(x) {
  if(all(!is.na(x)) & (x[1] != x[2])) { 
    return(TRUE) 
  } else {
    return(FALSE)
  }
})
## in case of discrepancies, use CGP as reference
tt <- celline.collection[ , "CGP.tissue.type"]
names(tt) <- rownames(celline.collection)
## read manual curation of tissue types for cell lines with missing tissue type in CGP and CCLE
match.tissue.ccle.cgp <- read.csv(file=file.path("matching_tissue_type_CCLE_CGP.csv"), stringsAsFactors=FALSE)
tt[as.character(match.tissue.ccle.cgp[ , "cellid"])] <- as.character(match.tissue.ccle.cgp[ , "new.tissue.type"])
## use CCLE tissue type for missing entries
tt[is.na(tt)] <- celline.collection[is.na(tt), "CCLE.tissue.type"]
tissue.cgp <- tt[rownames(data.cgp)]
tissue.ccle <- tt[rownames(data.ccle)]
## save cell line collection
celline.collection <- data.frame(celline.collection, "tissue.type"=tt)
write.csv(celline.collection, file=file.path(saveres, "cell_line_collection_all.csv"), row.names=FALSE)
## update sample information with tissue types
sampleinfo.cgp <- data.frame(sampleinfo.cgp, "tissue.type"=tissue.cgp)
sampleinfo.ccle <- data.frame(sampleinfo.ccle, "tissue.type"=tissue.ccle)

## mutation data
message("\tFormat (missense) mutation data")
ccell <- intersect(rownames(mutation.cgp), rownames(mutation.ccle))
cgene <- intersect(colnames(mutation.cgp), colnames(mutation.ccle))
mutation.cgp <- mutation.cgp[ccell, cgene, drop=FALSE]
mutation.ccle <- mutation.ccle[ccell, cgene, drop=FALSE]
write.csv(mutation.cgp, file=file.path(saveres, "mutation_cgp_common.csv"))
write.csv(mutation.ccle, file=file.path(saveres, "mutation_ccle_common.csv"))

## drug sensitivity measures for CGP
message("\tFormat drug sensitivity measures")

## IC50 in micro molar
message("\tIC50 (CGP)")

myx <- sapply(strsplit(colnames(drugpheno.cgp), "_"), function(x) { return(all(x[c(length(x)-1, length(x))] == c("IC", "50"))) })
ic50.cgp <- drugpheno.cgp[ ,myx,drop=FALSE]
nn <- dimnames(ic50.cgp)
nn[[2]] <- gsub("_IC_50", "", nn[[2]])
ic50.cgp <- apply(ic50.cgp, 2, as.numeric)
dimnames(ic50.cgp) <- nn
ic50.cgp <- exp(ic50.cgp)
## boxplot
oo <- order(apply(ic50.cgp, 2, median, na.rm=TRUE))
xx <- apply(-log10(ic50.cgp[ , oo, drop=FALSE] / 10^6), 2, function(x) { return(as.numeric(x[!is.na(x)])) })
pdf(file.path(saveres, "cgp_ic50_boxplot.pdf"), width=23, height=10)
par(las=3, mar=c(10, 4, 4, 2) + 0.1)
boxplot(xx, outline=FALSE, ylab="-log10(IC50)", main="Drugs sensitivity\nCGP")
dev.off()

## sensitivity calling using waterfall plot
ic50.call.cgp <- NULL
pdf(file.path(saveres, "cgp_ic50_sensitivity_calling_drugs.pdf"), width=5, height=10)
for(i in 1:ncol(ic50.cgp)) {
  ic50.call.cgp <- cbind(ic50.call.cgp, sensitivity.calling.waterfall(x=ic50.cgp[ ,i], type="ic50", intermediate.fold=c(4, 1.2, 1.2), cor.min.linear=0.95, plot=TRUE, name=sprintf("%s (CGP)", colnames(ic50.cgp)[i])))
}
dev.off()
dimnames(ic50.call.cgp) <- dimnames(ic50.cgp)
## count how many drugs have enough resistant and sensitive cell lines
summcall <- t(apply(ic50.call.cgp, 2, function(x) {
  tt <- table(x)
  tab <- c("resistant"=0, "intermediate"=0, "sensitive"=0)
  tab[names(tt)] <- tt
  return(tab)
  # return(ifelse(tab["resistant"] >= minn && tab["sensitive"] >= minn, TRUE, FALSE))
}))
rr <- apply(summcall, 1, function(x, minn) { return(ifelse(x["resistant"] >= minn && x["sensitive"] >= minn, TRUE, FALSE)) }, minn=5)
message(sprintf("Sensitivity calling using waterfall method:\n %i/%i drugs have at least 5 cell lines either resistant or sensitive", sum(rr), length(rr)))

## filtering based on concentration range
## ic50 larger or equal to the maximum tested drug concentration are filtered out
ic50.filt.cgp <- matrix(FALSE, nrow=nrow(ic50.cgp), ncol=ncol(ic50.cgp), dimnames=dimnames(ic50.cgp))
drugc <- drugconc.cgp[is.element(drugconc.cgp[ , "cellid"], rownames(ic50.cgp)) & is.element(drugconc.cgp[ , "drugid"], colnames(ic50.cgp)), , drop=FALSE]
maxdose <- drugc[ , "max.Dose.uM"]
ffilt <- matrix(NA, nrow=nrow(ic50.cgp), ncol=ncol(ic50.cgp), dimnames=dimnames(ic50.cgp))
ffilt[as.matrix(drugc[ , c("cellid", "drugid")])] <- maxdose
ic50.filt.cgp[!is.na(ic50.cgp) & ic50.cgp < ffilt] <- TRUE
dimnames(ic50.filt.cgp) <- dimnames(ic50.cgp)
## ic50.filt is set to TRUE if the IC50 measurement passes the filter

## activity area
message("\tAUC (CGP)")

## continuous values
myx <- sapply(strsplit(colnames(drugpheno.cgp), "_"), function(x) { return(all(x[c(length(x))] == c("AUC"))) })
auc.cgp <- drugpheno.cgp[ , myx, drop=FALSE]
nn <- dimnames(auc.cgp)
nn[[2]] <- gsub("_AUC", "", nn[[2]])
auc.cgp <- apply(auc.cgp, 2, as.numeric)
## AUC for sensitivity
auc.cgp <- 1 - auc.cgp
dimnames(auc.cgp) <- nn
oo <- order(apply(auc.cgp, 2, median, na.rm=TRUE))
xx <- apply(auc.cgp[ , oo, drop=FALSE], 2, function(x) { return(as.numeric(x[!is.na(x)])) })
## diagnostic plots
pdf(file.path(saveres, "cgp_auc_boxplot.pdf"), width=23, height=10)
par(las=3, mar=c(10, 4, 4, 2) + 0.1)
boxplot(xx, outline=FALSE, ylab="AUC", main="Drugs sensitivity\nCGP", ylim=c(0, 1))
dev.off()

## sensitivity calling using waterfall plot
auc.call.cgp <- NULL
pdf(file.path(saveres, "cgp_auc_sensitivity_calling_drugs.pdf"), width=5, height=10)
for(i in 1:ncol(auc.cgp)) {
  auc.call.cgp <- cbind(auc.call.cgp, sensitivity.calling.waterfall(x=auc.cgp[ ,i], type="actarea", intermediate.fold=c(4, 1.2, 1.2), cor.min.linear=0.95, plot=TRUE, name=sprintf("%s (CGP)", colnames(auc.cgp)[i])))
}
dev.off()
dimnames(auc.call.cgp) <- dimnames(auc.cgp)
## count how many drugs have enough resistant and sensitive cell lines
summcall <- t(apply(auc.call.cgp, 2, function(x) {
  tt <- table(x)
  tab <- c("resistant"=0, "intermediate"=0, "sensitive"=0)
  tab[names(tt)] <- tt
  return(tab)
  # return(ifelse(tab["resistant"] >= minn && tab["sensitive"] >= minn, TRUE, FALSE))
}))
rr <- apply(summcall, 1, function(x, minn) { return(ifelse(x["resistant"] >= minn && x["sensitive"] >= minn, TRUE, FALSE)) }, minn=5)
message(sprintf("Sensitivity calling using waterfall method:\n %i/%i drugs have at least 5 cell lines either resistant or sensitive", sum(rr), length(rr)))

## IC50 in micro molar
message("\tIC50 (CCLE)")

myx <- sapply(strsplit(colnames(drugpheno.ccle), "_"), function(x) { return(x[length(x)] == "IC50") })
ic50.ccle <- drugpheno.ccle[ ,myx,drop=FALSE]
nn <- dimnames(ic50.ccle)
nn[[2]] <- gsub("_IC50", "", nn[[2]])
ic50.ccle <- apply(ic50.ccle, 2, as.numeric)
dimnames(ic50.ccle) <- nn
## boxplot
oo <- order(apply(ic50.ccle, 2, median, na.rm=TRUE))
xx <- apply(-log10(ic50.ccle[ , oo, drop=FALSE] / 10^6), 2, function(x) { return(as.numeric(x[!is.na(x)])) })
pdf(file.path(saveres, "ccle_ic50_boxplot.pdf"), width=10, height=10)
par(las=3, mar=c(10, 4, 4, 2) + 0.1)
boxplot(xx, outline=FALSE, ylab="-log10(IC50)", main="Drugs sensitivity\nCCLE")
dev.off()

## sensitivity calling using waterfall plot
## 8in CCLE IC50 = 8 is a placeholder as drug concentration was not enough to yield 50% of growth inhibition, so IC50 cannot be  estimated
## here we assume that these cells are resistant
resix <- !is.na(ic50.ccle) & ic50.ccle >= 8
ic50t.ccle <- ic50.ccle
ic50t.ccle[resix] <- NA
ic50.call.ccle <- NULL
pdf(file.path(saveres, "ccle_ic50_sensitivity_calling_drugs.pdf"), width=5, height=10)
for(i in 1:ncol(ic50.ccle)) {
  ic50.call.ccle <- cbind(ic50.call.ccle, sensitivity.calling.waterfall(x=ic50t.ccle[ ,i], type="ic50", intermediate.fold=c(4, 1.2, 1.2), cor.min.linear=0.95, plot=TRUE, name=sprintf("%s (CCLE)", gsub("drugid_", "", colnames(ic50.ccle)[i]))))
}
dev.off()
dimnames(ic50.call.ccle) <- dimnames(ic50t.ccle)
## set the initial IC50 = 8 to resistant
ic50.call.ccle[resix] <- "resistant"
## count how many drugs have enough resistant and sensitive cell lines
summcall <- t(apply(ic50.call.ccle, 2, function(x) {
  tt <- table(x)
  tab <- c("resistant"=0, "intermediate"=0, "sensitive"=0)
  tab[names(tt)] <- tt
  return(tab)
  # return(ifelse(tab["resistant"] >= minn && tab["sensitive"] >= minn, TRUE, FALSE))
}))
rr <- apply(summcall, 1, function(x, minn) {return(ifelse(x["resistant"] >= minn && x["sensitive"] >= minn, TRUE, FALSE))}, minn=5)
message(sprintf("Sensitivity calling using waterfall method:\n %i/%i drugs have at least 5 cell lines either resistant or sensitive", sum(rr), length(rr)))

## filtering based on concentration range
## ic50 larger or equal to the maximum tested drug concentration are filtered out
ic50.filt.ccle <- matrix(FALSE, nrow=nrow(ic50.ccle), ncol=ncol(ic50.ccle), dimnames=dimnames(ic50.ccle))
drugc <- drugconc.ccle[is.element(drugconc.ccle[ , "cellid"], rownames(ic50.ccle)) & is.element(drugconc.ccle[ , "drugid"], colnames(ic50.ccle)), , drop=FALSE]
maxdose <- apply(drugc[ , grep("^Dose", colnames(drugconc.ccle))], 1, max, na.rm=TRUE)
ffilt <- matrix(NA, nrow=nrow(ic50.ccle), ncol=ncol(ic50.ccle), dimnames=dimnames(ic50.ccle))
ffilt[as.matrix(drugc[ , c("cellid", "drugid")])] <- maxdose
ic50.filt.ccle[!is.na(ic50.ccle) & ic50.ccle < ffilt] <- TRUE
dimnames(ic50.filt.ccle) <- dimnames(ic50.ccle)
## ic50.filt is set to TRUE if the IC50 measurement passes the filter


## activity area
message("\tAUC (CCLE)")

## Note: to compute the activity area use the activity data (median)
## SUM( - Activity Data ) / 100
myx <- sapply(strsplit(colnames(drugpheno.ccle), "_"), function(x) { return(x[length(x)] == "ActivityArea") })
auc.ccle <- drugpheno.ccle[ ,myx,drop=FALSE]
nn <- dimnames(auc.ccle)
nn[[2]] <- gsub("_ActivityArea", "", nn[[2]])
auc.ccle <- apply(auc.ccle, 2, as.numeric)
dimnames(auc.ccle) <- nn
## division by the number of concentrations tested
myx <- sapply(strsplit(colnames(drugpheno.ccle), "_"), function(x) { return(x[length(x)] == "Doses") })
ndoses.ccle <- drugpheno.ccle[ ,myx,drop=FALSE]
ndoses <- apply(ndoses.ccle, 2, function(x) {
  return(sapply(x, function(x) {
    return(sapply(strsplit(as.character(x), ","), function(x) { if(is.na(x[1])) { return(NA) } else { return(length(x)) } }))
  }))
})
auc.ccle <- (auc.ccle / ndoses)
## boxplot
oo <- order(apply(auc.ccle, 2, median, na.rm=TRUE))
xx <- apply(auc.ccle[ , oo, drop=FALSE], 2, function(x) { return(as.numeric(x[!is.na(x)])) })
pdf(file.path(saveres, "ccle_auc_boxplot.pdf"), width=10, height=10)
par(las=3, mar=c(10, 4, 4, 2) + 0.1)
boxplot(xx, outline=FALSE, ylab="AUC", main="Drugs sensitivity\nCCLE")
dev.off()
## sensitivity calling using waterfall plot
auc.call.ccle <- NULL
pdf(file.path(saveres, "ccle_auc_sensitivity_calling_drugs.pdf"), width=5, height=10)
for(i in 1:ncol(auc.ccle)) {
  auc.call.ccle <- cbind(auc.call.ccle, sensitivity.calling.waterfall(x=auc.ccle[ ,i], type="actarea", intermediate.fold=c(4, 1.2, 1.2), cor.min.linear=0.95, plot=TRUE, name=sprintf("%s (CCLE)", gsub("drugid_", "", colnames(auc.ccle)[i]))))
}
dev.off()
dimnames(auc.call.ccle) <- dimnames(auc.ccle)
## count how many drugs have enough resistant and sensitive cell lines
summcall <- t(apply(auc.call.ccle, 2, function(x) {
  tt <- table(x)
  tab <- c("resistant"=0, "intermediate"=0, "sensitive"=0)
  tab[names(tt)] <- tt
  return(tab)
  # return(ifelse(tab["resistant"] >= minn && tab["sensitive"] >= minn, TRUE, FALSE))
}))
rr <- apply(summcall, 1, function(x, minn) { return(ifelse(x["resistant"] >= minn && x["sensitive"] >= minn, TRUE, FALSE)) }, minn=5)
message(sprintf("Sensitivity calling using waterfall method:\n %i/%i drugs have at least 5 cell lines either resistant or sensitive", sum(rr), length(rr)))


## IC50 in micro molar
message("\tIC50 (GSK)")
## transform drug concentration into Molar
ic50.gsk <- drugpheno.gsk
## boxplot
oo <- order(apply(ic50.gsk, 2, median, na.rm=TRUE))
xx <- apply(-log10(ic50.gsk[ , oo, drop=FALSE] / 10^6), 2, function(x) { return(as.numeric(x[!is.na(x)])) })
pdf(file.path(saveres, "gsk_ic50_boxplot.pdf"), width=23, height=10)
par(las=3, mar=c(10, 4, 4, 2) + 0.1)
boxplot(xx, outline=FALSE, ylab="-log10(IC50)", main="Drugs sensitivity\nGSK")
dev.off()

## filtering based on concentration range
## ic50 larger or equal to the maximum tested drug concentration are filtered out
ic50.filt.gsk <- matrix(FALSE, nrow=nrow(ic50.gsk), ncol=ncol(ic50.gsk), dimnames=dimnames(ic50.gsk))
drugc <- drugconc.gsk[is.element(drugconc.gsk[ , "cellid"], rownames(ic50.gsk)) & is.element(drugconc.gsk[ , "drugid"], colnames(ic50.gsk)), , drop=FALSE]
maxdose <- apply(drugc[ , grep("^Dose", colnames(drugconc.gsk))], 1, function(x) { return(max(as.numeric(x), na.rm=TRUE)) })
ffilt <- matrix(NA, nrow=nrow(ic50.gsk), ncol=ncol(ic50.gsk), dimnames=dimnames(ic50.gsk))
ffilt[as.matrix(drugc[ , c("cellid", "drugid")])] <- maxdose
ic50.filt.gsk[!is.na(ic50.gsk) & ic50.gsk < ffilt] <- TRUE
dimnames(ic50.filt.gsk) <- dimnames(ic50.gsk)
## ic50.filt is set to TRUE if the IC50 measurement passes the filter

## save full datasets
save(list=c("data.ccle", "data.cgp", "data.gsk", "mutation.cgp", "mutation.ccle", "annot.cgp", "annot.ccle", "annot.gsk", "druginfo.ccle", "druginfo.cgp", "ic50.ccle", "ic50.filt.ccle", "ic50.cgp", "ic50.filt.cgp", "ic50.call.cgp", "ic50.call.ccle", "ic50.gsk", "ic50.filt.gsk", "auc.cgp", "auc.ccle", "auc.call.cgp", "auc.call.ccle", "sampleinfo.cgp", "sampleinfo.ccle", "sampleinfo.gsk"), compress=TRUE, file=file.path(saveres, "CDRUG_cgp_ccle_gsk_full.RData"))

########################
## comparison between CCLE and CGP
## cell line id
cellid.common <- intersect(rownames(data.cgp), rownames(data.ccle))
## drug ids
## manual mapping CCLE vs CGP vs GSK
drug.map <- rbind(c("drugid_ERLOTINIB", "drugid_1", NA),
  c("drugid_LAPATINIB", "drugid_119", "drugid_LAPATINIB"),
  c("drugid_PHA665752", "drugid_6", NA),
  c("drugid_CRIZOTINIB", "drugid_37", NA),
  c("drugid_TAE684", "drugid_35", NA),
  c("drugid_VANDETANIB", NA, NA),
  c("drugid_NILOTINIB", "drugid_1013", NA),
  c("drugid_AZD0530", "drugid_38", NA),
  c("drugid_SORAFENIB", "drugid_30", NA),
  c("drugid_TKI258", NA, NA),
  c("drugid_PD0332991", "drugid_1054", NA),
  c("drugid_AEW541", NA, NA),
  c("drugid_RAF265", NA, NA),
  c("drugid_PLX4720", "drugid_1036", NA),
  c("drugid_PD0325901", "drugid_1060", NA),
  c("drugid_AZD6244", "drugid_1062", NA),
  c("drugid_NUTLIN3", "drugid_1047", NA),
  c("drugid_LBW242", NA, NA),
  c("drugid_17AAG", "drugid_1026", NA),
  c("drugid_L685458", NA, NA),
  c("drugid_PANOBINOSTAT", NA, NA),
  c("drugid_PACLITAXEL", "drugid_11", "drugid_PACLITAXEL"),
  # c("drugid_IRINOTECAN", "drugid_1003", NA),
  c("drugid_IRINOTECAN", NA, NA),
  c("drugid_TOPOTECAN", NA, NA)
)
colnames(drug.map) <- c("CCLE", "CGP", "GSK")
rownames(drug.map) <- gsub("drugid_", "", druginfo.ccle[drug.map[ ,"CCLE"], "drugid"])

## use CCLE drug names for consistency
iix <- drug.map[complete.cases(drug.map[ , c("CGP", "CCLE")]), "CGP"]
names(iix) <- drug.map[complete.cases(drug.map[ , c("CGP", "CCLE")]), "CCLE"]
## CGP
druginfo.cgp <- druginfo.cgp[iix, , drop=FALSE]
ic50.cgp <- ic50.cgp[ , iix, drop=FALSE]
ic50.filt.cgp <- ic50.filt.cgp[ , iix, drop=FALSE]
ic50.call.cgp <- ic50.call.cgp[ , iix, drop=FALSE]
auc.cgp <- auc.cgp[ , iix, drop=FALSE]
auc.call.cgp <- auc.call.cgp[ , iix, drop=FALSE]
rownames(druginfo.cgp) <- colnames(ic50.cgp) <- colnames(ic50.filt.cgp) <- colnames(ic50.call.cgp) <- colnames(auc.cgp) <- colnames(auc.call.cgp) <- names(iix)
## drug concentrations
drugconc.cgp <- drugconc.cgp[is.element(drugconc.cgp[ ,"drugid"], iix), , drop=FALSE]
drugconc.cgp[ , "drugid"] <- names(iix)[sapply(1:nrow(drugconc.cgp), function(x, y, z) { return(which(y[x] == z)) }, y=drugconc.cgp[ , "drugid"], z=iix)]
rownames(drugconc.cgp) <- paste(drugconc.cgp[ , "cellid"], drugconc.cgp[ , "drugid"], sep=".")
## CCLE
druginfo.ccle <- druginfo.ccle[names(iix), ,drop=FALSE]
ic50.ccle <- ic50.ccle[ , names(iix), drop=FALSE]
ic50.filt.ccle <- ic50.filt.ccle[ , names(iix), drop=FALSE]
ic50.call.ccle <- ic50.call.ccle[ , names(iix), drop=FALSE]
auc.ccle <- auc.ccle[ , names(iix), drop=FALSE]
auc.call.ccle <- auc.call.ccle[ , names(iix), drop=FALSE]
rownames(druginfo.ccle) <- colnames(ic50.ccle) <- colnames(ic50.filt.ccle) <- colnames(auc.ccle) <- colnames(auc.call.ccle) <- names(iix)
## drug concentrations
drugconc.ccle <- drugconc.ccle[is.element(drugconc.ccle[ ,"drugid"], names(iix)), , drop=FALSE]

## consider only common genes, drugs and cell lines between CCLE and CGP
cell.common <- intersect(rownames(data.ccle), rownames(data.cgp))
drug.common <- intersect(colnames(ic50.ccle), colnames(ic50.cgp))
gene.common <- intersect(colnames(data.ccle), colnames(data.cgp))

druginfo.ccle <- druginfo.ccle[drug.common, , drop=FALSE]
druginfo.cgp <- druginfo.cgp[drug.common, , drop=FALSE]
ic50.ccle <- ic50.ccle[cell.common, drug.common, drop=FALSE]
ic50.filt.ccle <- ic50.filt.ccle[cell.common, drug.common, drop=FALSE]
ic50.call.ccle <- ic50.call.ccle[cell.common, drug.common, drop=FALSE]
ic50.cgp <- ic50.cgp[cell.common, drug.common, drop=FALSE]
ic50.filt.cgp <- ic50.filt.cgp[cell.common, drug.common, drop=FALSE]
ic50.call.cgp <- ic50.call.cgp[cell.common, drug.common, drop=FALSE]
auc.ccle <- auc.ccle[cell.common, drug.common, drop=FALSE]
auc.call.ccle <- auc.call.ccle[cell.common, drug.common, drop=FALSE]
auc.cgp <- auc.cgp[cell.common, drug.common, drop=FALSE]
auc.call.cgp <- auc.call.cgp[cell.common, drug.common, drop=FALSE]
data.ccle <- data.ccle[cell.common, gene.common, drop=FALSE]
data.cgp <- data.cgp[cell.common, gene.common, drop=FALSE]
sampleinfo.ccle <- sampleinfo.ccle[cell.common, , drop=FALSE]
sampleinfo.cgp <- sampleinfo.cgp[cell.common, , drop=FALSE]
drugconc.cgp <- drugconc.cgp[is.element(drugconc.cgp[ , "cellid"], cell.common), , drop=FALSE]
drugconc.ccle <- drugconc.ccle[is.element(drugconc.ccle[ , "cellid"], cell.common), , drop=FALSE]
annot <- annot.ccle[gene.common, , drop=FALSE]

save(list=c("data.ccle", "data.cgp", "mutation.cgp", "mutation.ccle", "annot", "druginfo.ccle", "druginfo.cgp", "ic50.ccle", "ic50.filt.ccle", "ic50.cgp", "ic50.filt.cgp", "ic50.call.cgp", "ic50.call.ccle", "auc.cgp", "auc.ccle", "auc.call.cgp", "auc.call.ccle", "sampleinfo.cgp", "sampleinfo.ccle"), compress=TRUE, file=file.path(saveres, "CDRUG_cgp_ccle.RData"))



## end

