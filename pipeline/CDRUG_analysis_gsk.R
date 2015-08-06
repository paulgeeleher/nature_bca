########################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## February 11, 2013
########################


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

load(file.path(saveres, "CDRUG_cgp_ccle_gsk_full.RData"))

tissue.cgp <- as.character(sampleinfo.cgp[ , "tissue.type"])
tissue.gsk <- as.character(sampleinfo.gsk[ , "tissue.type"])
tissue.ccle <- as.character(sampleinfo.ccle[ , "tissue.type"])
tissue <- fold(union, tissue.cgp, tissue.ccle, tissue.gsk)
tissuen <- sort(unique(as.character(tissue)))
ic50f.ccle <- ic50.ccle <- -log10(ic50.ccle / 10^6)
ic50f.cgp <- ic50.cgp <- -log10(ic50.cgp / 10^6)
ic50f.gsk <- ic50.gsk <- -log10(ic50.gsk / 10^6)
## filter out low quality ic50 at high concentrations
ic50f.ccle[!ic50.filt.ccle] <- NA
ic50f.cgp[!ic50.filt.cgp] <- NA
ic50f.gsk[!ic50.filt.gsk] <- NA

drug.map <- rbind(c("drugid_LAPATINIB", "drugid_119", "drugid_LAPATINIB"),
  c("drugid_PACLITAXEL", "drugid_11", "drugid_PACLITAXEL")
)
colnames(drug.map) <- c("CCLE", "CGP", "GSK")
rownames(drug.map) <- gsub("drugid_", "", druginfo.ccle[drug.map[ ,"CCLE"], "drugid"])

## CGP vs GSK

###############################################################################
###############################################################################
## Correlation of sensitivity measurements
###############################################################################
###############################################################################

## all IC50s
for(i in 1:nrow(drug.map)) {
  pdf(file.path(saveres, sprintf("ic50_%s_ccle_cgp_gsk_paper.pdf", rownames(drug.map)[i])), height=6, width=12)
  par(mfrow=c(1, 2), cex=0.8, las=1)
  ## CGP
  # xxlim <- yylim <- round(range(c(ic50.cgp[, i], ic50.ccle[, i]), na.rm=TRUE) * 10) / 10
  cell.common <- intersect(rownames(data.cgp), rownames(data.gsk))
  ic50.1 <- ic50.cgp[cell.common, drug.map[i, "CGP"]]
  ic50.2 <- ic50.gsk[cell.common, drug.map[i, "GSK"]]
  yylim <- round(range(ic50.1, na.rm=TRUE) * 10) / 10
  xxlim <- round(range(ic50.2, na.rm=TRUE) * 10) / 10
  nnn <- sum(complete.cases(ic50.1, ic50.2))
  if(nnn >= minsample) {
    cc <- cor.test(ic50.1, ic50.2, method="spearman", use="complete.obs", alternative="greater")
  } else {
    cc <- list("estimate"=NA, "p.value"=NA)
  }
  myScatterPlot(x=ic50.2, y=ic50.1, xlab="-log10 IC50 (GSK)", ylab="-log10 IC50 (CGP)", main=sprintf("%s\nCGP vs. GSK", gsub("drugid_", "", rownames(drug.map)[i])), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75)
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
  legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, nnn), text.font=2)
  ## CCLE
  cell.common <- intersect(rownames(data.ccle), rownames(data.gsk))
  ic50.1 <- ic50.ccle[cell.common, drug.map[i, "CCLE"]]
  ic50.2 <- ic50.gsk[cell.common, drug.map[i, "GSK"]]
  yylim <- round(range(ic50.1, na.rm=TRUE) * 10) / 10
  xxlim <- round(range(ic50.2, na.rm=TRUE) * 10) / 10
  nnn <- sum(complete.cases(ic50.1, ic50.2))
  if(nnn >= minsample) {
    cc <- cor.test(ic50.1, ic50.2, method="spearman", use="complete.obs", alternative="greater")
  } else {
    cc <- list("estimate"=NA, "p.value"=NA)
  }
  myScatterPlot(x=ic50.2, y=ic50.1, xlab="-log10 IC50 (GSK)", ylab="-log10 IC50 (CCLE)", main=sprintf("%s\nCCLE vs. GSK", gsub("drugid_", "", rownames(drug.map)[i])), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75)
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
  legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, nnn), text.font=2)
  dev.off()
}

## filtered IC50s
for(i in 1:nrow(drug.map)) {
  pdf(file.path(saveres, sprintf("ic50_filt_%s_ccle_cgp_gsk.pdf", rownames(drug.map)[i])), height=6, width=12)
  par(mfrow=c(1, 2), cex=0.8, las=1)
  ## CGP
  # xxlim <- yylim <- round(range(c(ic50.cgp[, i], ic50.ccle[, i]), na.rm=TRUE) * 10) / 10
  cell.common <- intersect(rownames(data.cgp), rownames(data.gsk))
  ic50.1 <- ic50f.cgp[cell.common, drug.map[i, "CGP"]]
  ic50.2 <- ic50f.gsk[cell.common, drug.map[i, "GSK"]]
  yylim <- round(range(ic50.1, na.rm=TRUE) * 10) / 10
  xxlim <- round(range(ic50.2, na.rm=TRUE) * 10) / 10
  nnn <- sum(complete.cases(ic50.1, ic50.2))
  if(nnn >= minsample) {
    cc <- cor.test(ic50.1, ic50.2, method="spearman", use="complete.obs", alternative="greater")
  } else {
    cc <- list("estimate"=NA, "p.value"=NA)
  }
  myScatterPlot(x=ic50.2, y=ic50.1, xlab="-log10 IC50 (GSK)", ylab="-log10 IC50 (CGP)", main=sprintf("%s\nCGP vs. GSK", gsub("drugid_", "", rownames(drug.map)[i])), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75)
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
  legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, nnn), text.font=2)
  ## CCLE
  cell.common <- intersect(rownames(data.ccle), rownames(data.gsk))
  ic50.1 <- ic50f.ccle[cell.common, drug.map[i, "CCLE"]]
  ic50.2 <- ic50f.gsk[cell.common, drug.map[i, "GSK"]]
  yylim <- round(range(ic50.1, na.rm=TRUE) * 10) / 10
  xxlim <- round(range(ic50.2, na.rm=TRUE) * 10) / 10
  nnn <- sum(complete.cases(ic50.1, ic50.2))
  if(nnn >= minsample) {
    cc <- cor.test(ic50.1, ic50.2, method="spearman", use="complete.obs", alternative="greater")
  } else {
    cc <- list("estimate"=NA, "p.value"=NA)
  }
  myScatterPlot(x=ic50.2, y=ic50.1, xlab="-log10 IC50 (GSK)", ylab="-log10 IC50 (CCLE)", main=sprintf("%s\nCCLE vs. GSK", gsub("drugid_", "", rownames(drug.map)[i])), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75)
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
  legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, nnn), text.font=2)
  dev.off()
}



## end




