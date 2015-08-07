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

ccelln <- intersect(rownames(sampleinfo.cgp), rownames(sampleinfo.ccle))
tt <- table(sampleinfo.cgp[ccelln, "tissue.type"], sampleinfo.cgp[ccelln, "tissue.type"])
if(any(tt[upper.tri(tt)] > 0) || any(tt[lower.tri(tt)] > 0)) { warning("Discrepancies in tissue type classification between CGP and CCLE") }
mm <- cbind("Tissue type"=names(diag(tt)), "Number of cell lines"=diag(tt))
mm <- mm[mm[ , 2] != 0, , drop=FALSE]
xtable::print.xtable(xtable::xtable(mm), include.rownames=FALSE, floating=FALSE, file=file.path(saveres, "tissue_type_cgp_ccle_paper.tex"), append=FALSE)

utissue <- table(as.character(tissue))
coltissue <- coltissuet <- factor(names(utissue))
# levels(coltissue) <- gplots::rich.colors(length(utissue))
levels(coltissue) <- rainbow(n=length(utissue), s=0.5, v=1)
levels(coltissuet) <- rainbow(n=length(utissue), s=0.75, v=1)
coltissue <- as.character(coltissue)
coltissuet <- as.character(coltissuet)
names(coltissue) <- names(coltissuet) <- names(utissue)

###############################################################################
###############################################################################
## Correlation of sensitivity measurements
###############################################################################
###############################################################################

## create tables summarizing all the correlations
drugsn <- gsub("drugid_", "", colnames(ic50.ccle))
tt <- matrix(NA, nrow=length(drugsn), ncol=6, dimnames=list(drugsn, c("drug.sensitivity", "drug.sensitivity.filt", "gene.drug", "gene.drug.filt", "go.drug", "go.drug.filt")))
correlations <- list("ic50"=tt, "ic50.call"=tt, "auc"=tt, "auc.call"=tt)
## correlation statistics
tt <- matrix(NA, nrow=length(drugsn), ncol=5, dimnames=list(drugsn, c("rho", "lower", "upper", "p", "n")))
correlations.stats <- list("ic50"=tt, "ic50.filt"=tt, "ic50.call"=tt, "auc"=tt, "auc.call"=tt)

########################
## IC50

message("Correlation IC50 measures:")
## all IC50s
## consistency between IC50 with spearman correlation
pdf(file.path(saveres, "ic50_spearman_ccle_cgp_pres.pdf"), height=9, width=16)
par(mfrow=c(3, 5), cex=0.8, las=1)
for(i in 1:ncol(ic50.ccle)) {
  # xxlim <- yylim <- round(range(c(ic50.cgp[, i], ic50.ccle[, i]), na.rm=TRUE) * 10) / 10
  xxlim <- round(range(ic50.cgp[, i], na.rm=TRUE) * 10) / 10
  yylim <- round(range(ic50.ccle[, i], na.rm=TRUE) * 10) / 10
  nnn <- sum(complete.cases(ic50.ccle[, i], ic50.cgp[, i]))
  if(nnn >= minsample) {
    cc <- cor.test(ic50.ccle[, i], ic50.cgp[, i], method="spearman", use="complete.obs", alternative="greater")
  } else {
    cc <- list("estimate"=NA, "p.value"=NA)
  }
  # plot(x=ic50.cgp[ ,i], y=ic50.ccle[ ,i], xlab="-log10 IC50 (CGP)", ylab="-log10 IC50 (CCLE)", main=gsub("drugid_", "", colnames(ic50.ccle)[i]), pch=16, col=rgb(0, 200, 50, 100, maxColorValue=255))
  par(mar=c(4, 4, 3, 1) + 0.1)
  myScatterPlot(x=ic50.cgp[, i], y=ic50.ccle[, i], xlab=ifelse(i > 10, "-log10 IC50 (CGP)", ""), ylab=ifelse((i %% 5) == 1, "-log10 IC50 (CCLE)", ""), main=gsub("drugid_", "", colnames(ic50.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, nnn), text.font=2)
  correlations[["ic50"]][i, "drug.sensitivity"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
  correlations.stats[["ic50"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
}
dev.off()

## IC50s per tissue types
correlations.stats.ic50.tissue <- array(NA, dim=c(length(tissuen), length(drugn), 5), dimnames=list(tissuen, gsub("drugid_", "", drugn), c("rho", "upper", "lower", "p", "n")))
for(ii in 1:length(tissuen)) {
  for(i in 1:ncol(ic50.ccle)) {
    iix <- !is.na(tissue) & tissue == tissuen[ii]
    nnn <- sum(complete.cases(ic50.ccle[iix, i], ic50.cgp[iix, i]))
    if(nnn >= minsample) {
      cc <- cor.test(ic50.ccle[iix, i], ic50.cgp[iix, i], method="spearman", use="complete.obs", alternative="greater")
    } else {
      cc <- list("estimate"=NA, "p.value"=NA)
    }
    ## correlation statistics
    cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
    correlations.stats.ic50.tissue[ii, i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
  }
}

## boxplot for paper
pdf(file.path(saveres, sprintf("boxplot_ic50_ccle_cgp_tissue.pdf")), height=10, width=14)
## correlations for each drug per tissue type
tt <- correlations.stats.ic50.tissue[ , , "rho"]
## remove the non significant p-values
# tt[correlations.stats.ic50.tissue[ , , "p"] >= 0.05] <- NA
ll <- c(list("all_tissues"=correlations.stats[["ic50"]][ , "rho"]), lapply(apply(tt, 1, list), function(x) { return(x[[1]]) }))
ll <- lapply(ll, function(x) { return(x[!is.na(x)]) })
# ll <- ll[sapply(ll, length) > 0]
myx <- which(sapply(ll, length) < 3)
wt <- kruskal.test(x=ll[sapply(ll, length) > 0])
pp <- NULL
for(j in myx) { for(jj in ll[[j]]) { pp <- rbind(pp, c(j, jj)) } }
ll[sapply(ll, length) < 3] <- NA
par(las=2, mar=c(12, 4, 4, 2) + 0.1, xaxt="n")
mp <- boxplot(ll, outline=FALSE, ylab="Rs", main="Correlations of drug sensitivity (IC50) across tissue types\nCCLE vs. CGP", border="black", ylim=c(-1, 1), col=c("lightgrey", coltissue))
points(x=pp[ , 1], y=pp[ , 2], col=coltissue[pp[ , 1]], pch=20)
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.02), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", coltissuet))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
dev.off()

## boxplot for paper
pdf(file.path(saveres, sprintf("boxplot2_ic50_ccle_cgp_tissue_paper.pdf")), height=8, width=14)
## correlations for each drug per tissue type
tt <- correlations.stats.ic50.tissue[ , , "rho"]
## remove the non significant p-values
# tt[correlations.stats.ic50.tissue[ , , "p"] >= 0.05] <- NA
ll <- c(list("all_tissues"=correlations.stats[["ic50"]][ , "rho"]), lapply(apply(tt, 1, list), function(x) { return(x[[1]]) }))
ll <- lapply(ll, function(x) { return(x[!is.na(x)]) })
# ll <- ll[sapply(ll, length) > 0]
myx <- which(sapply(ll, length) < 3)
wt <- kruskal.test(x=ll[sapply(ll, length) > 0])
ll[sapply(ll, length) < 3] <- NA
par(las=1, mar=c(12, 4, 4, 2) + 0.1, xaxt="n")
mp <- boxplot(ll, outline=FALSE, ylab="Rs", main="Correlations of drug sensitivity (IC50) across tissue types\nCCLE vs. CGP", border="grey50", ylim=c(-1, 1), col="white", range=0)
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=c("lightgrey", coltissuet))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=1:length(ll))
points(x=jitter(unlist(xx), 1), y=unlist(ll), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", coltissuet))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
dev.off()

## test whether correlations are significantly higher in specific tissue types compared to all tissues combined
cor.diff.test <- NULL
ll <- c(list("all_tissues"=correlations.stats[["ic50"]][ , "rho"]), lapply(apply(correlations.stats.ic50.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ddnn <- names(ll[[1]])
for(i in 2:length(ll)) {
  if(sum(complete.cases(ll[[1]], y=ll[[i]])) < 3) {
   wt <- list("estimate"=NA, "p.value"=NA) 
  } else {
    wt <- wilcox.test(x=ll[[1]], y=ll[[i]], alternative="less", conf.int=TRUE)
  }
  cor.diff.test <- rbind(cor.diff.test, c(wt$estimate, wt$p.value))
}
dimnames(cor.diff.test) <- list(names(ll)[-1], c("difference", "p"))
# cor.diff.test <- cbind(cor.diff.test, "lfdr"=fdrtool::fdrtool(cor.diff.test[ , 2], statistic="pvalue", verbose=FALSE)$lfdr)
cor.diff.test <- cbind(cor.diff.test, "fdr"=p.adjust(cor.diff.test[ , "p"], method="fdr"))
myx <- complete.cases(cor.diff.test) & (cor.diff.test[ , "p"] < 0.05 | cor.diff.test[ , "fdr"] < myfdr)
tt <- cor.diff.test[myx, , drop=FALSE]
tt <- tt[order(tt[ , "p"], decreasing=FALSE), , drop=FALSE]
message("Are correlations between IC50s in some tissue types significantly higher than in all tissues combined?")
print(tt)

## discretized IC50s (sensitivity calling)

## we used Cohen Kappa coefficient (k) (67), as implemented in the R package epibasix ; k ranges from 0 to 1, with 0 indicating no rela- tion and 1 indicating a perfect concordance. Typically qualitative descriptions are associated with intervals [k ≤ 0.20, slight agreement; 0.20 < k ≤ 0.40, fair agreement; 0.40 < k ≤ 0.60, moderate agreement, 0.60 < k ≤ 0.80, substantial agreement; and 0.80 < k ≤ 0.99, almost perfect agreement, as described in Weigelt et al. (22)].

message("Correlation IC50 sensitivity calls:")

pdf(file.path(saveres, "ic50_sensitivity_calling_ccle_cgp.pdf"), width=9, height=11)
par(mfrow=c(5, 3), mar=c(8, 4, 4, 2) + 0.1)
mykap <- NULL
for (i in 1:ncol(ic50.call.ccle)) {
  ff1 <- factor(x=ic50.call.ccle[ , i], levels=c("resistant", "intermediate", "sensitive"))
  ff2 <- factor(x=ic50.call.cgp[ , i], levels=c("resistant", "intermediate", "sensitive"))
  levels(ff1)[2] <- NA
  levels(ff2)[2] <- NA
  tt <- table("CCLE"=ff1, "CGP"=ff2)
  tts <- vcd::assocstats(tt)
  ## check if all elements are on the diagonal
  ttt <- tt
  diag(ttt) <- 0
  if(sum(ttt) == 0) {
    rr <- rr.l <- rr.u <- 1
  } else {
    err <- try(rr <- epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if(class(err) == "try-error") { rr <- rr.l <- rr.u <- NA }
    rr.l <- rr$CIL
    rr.u <- rr$CIU
    rr <- rr$kappa
  }
  mykap <- c(mykap, rr)
  plot.new()
  plotrix::addtable2plot(x=par("usr")[1] + 0.15, y=par("usr")[4] - 0.5, xjust=0, yjust=0, table=as.matrix(tt), bty="o", display.rownames=TRUE, display.colnames=TRUE, hlines=FALSE, vlines=FALSE, title="CCLE   vs  CGP", cex=1.25)
  title(main=sprintf("IC50 sensitivity calling\n%s", gsub("drugid_", "", colnames(auc.call.ccle)[i])), sub=sprintf("Kappa=%.2g, 95%%CI [%.2g,%.2g], p=%.1E", rr, rr.l, rr.u, tts$chisq_tests[1,3]), cex=1)
}
names(mykap) <- colnames(auc.call.ccle)
dev.off()
## with intermediate
pdf(file.path(saveres, "ic50_sensitivity_calling_intermediate_ccle_cgp_paper.pdf"), width=8, height=9)
par(mfrow=c(5, 3), mar=c(5, 4, 4, 2) + 0.1)
mykap <- NULL
for (i in 1:ncol(ic50.call.ccle)) {
  ff1 <- factor(x=ic50.call.ccle[ , i], levels=c("resistant", "intermediate", "sensitive"))
  ff2 <- factor(x=ic50.call.cgp[ , i], levels=c("resistant", "intermediate", "sensitive"))
  levels(ff1) <- c("res", "inter", "sens")
  levels(ff2) <- c("res", "inter", "sens")
  nnn <- sum(complete.cases(ff1, ff2))
  tt <- table("CCLE"=ff1, "CGP"=ff2)
  tts <- vcd::assocstats(tt)
  ## check if all elements are on the diagonal
  ttt <- tt
  diag(ttt) <- 0
  if(sum(ttt) == 0) {
    rr <- rr.l <- rr.u <- 1
  } else {
    err <- try(rr <- epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if(class(err) == "try-error") { rr <- rr.l <- rr.u <- NA }
    rr.l <- rr$CIL
    rr.u <- rr$CIU
    rr <- rr$kappa
  }
  mykap <- c(mykap, rr)
  plot.new()
  plotrix::addtable2plot(x=par("usr")[1] + 0.1, y=par("usr")[4] - 0.4, xjust=0, yjust=0, table=as.matrix(tt), bty="o", display.rownames=TRUE, display.colnames=TRUE, hlines=FALSE, vlines=FALSE, title="CCLE   vs   CGP", cex=1.25)
  title(main=sprintf("IC50 sensitivity calling\n%s", gsub("drugid_", "", colnames(ic50.call.ccle)[i])), sub=sprintf("Kappa=%.2g, 95%%CI [%.2g,%.2g], p=%.1E", rr, rr.l, rr.u, tts$chisq_tests[1,3]))
  correlations[["ic50.call"]][i, "drug.sensitivity"] <- rr
  ## statistics
  correlations.stats[["ic50.call"]][i, ] <- c(rr, rr.l, rr.u, tts$chisq_tests[1,3], nnn)
}
names(mykap) <- colnames(ic50.call.ccle)
dev.off()
## histogram of Kappa coefficients
pdf(file.path(saveres, "ic50_sensitivity_calling_kappa_ccle_cgp.pdf"), width=7, height=7)
hist(mykap, main=sprintf("Sensitivity calling based on IC50\nKappa coefficients (CCLE vs CGP)", xlab="Kappa"))
dev.off()


########################
## AUC

message("Correlation AUC measures")
## consistency between AUC with spearman correlation
pdf(file.path(saveres, "auc_spearman_ccle_cgp_pres.pdf"), height=9, width=16)
par(mfrow=c(3, 5), cex=0.8, las=1)
for(i in 1:ncol(auc.ccle)) {
  xxlim <- yylim <- round(range(c(auc.cgp[, i], auc.ccle[, i]), na.rm=TRUE) * 10) / 10
  nnn <- sum(complete.cases(auc.ccle[, i], auc.cgp[, i]))
  if(nnn >= minsample) {
    cc <- cor.test(auc.ccle[, i], auc.cgp[, i], method="spearman", use="complete.obs", alternative="greater")
  } else {
    cc <- list("estimate"=NA, "p.value"=NA)
  }
  par(mar=c(4, 4, 3, 1) + 0.1)
  myScatterPlot(x=auc.cgp[, i], y=auc.ccle[, i], xlab=ifelse(i > 10, "AUC (CGP)", ""), ylab=ifelse((i %% 5) == 1, "AUC (CCLE)", ""), main=gsub("drugid_", "", colnames(auc.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, sum(complete.cases(auc.ccle[ ,i], auc.cgp[ ,i]))), text.font=2)
  correlations[["auc"]][i, "drug.sensitivity"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
  correlations.stats[["auc"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
}
dev.off()


######################################################################
############## Begin Geeleher et. al Code ############################
######################################################################

## Calculate the AUC correlations *between* cell lines (rather than across cell lines).
corOutSpearman <- numeric()
corOutPearson <- numeric()
for(j in 1:nrow(auc.ccle))
{
  if((sum(is.na(auc.ccle[j,])) < 12) & (sum(is.na(auc.cgp[j,])) < 12)) #ensure that we reach the minimum number of samples to run a Spearman's correlation test
  {
    corOutSpearman[j] <- cor.test(auc.ccle[j,], auc.cgp[j,], method="spearman")$estimate
    corOutPearson[j] <- cor.test(auc.ccle[j,], auc.cgp[j,], method="pearson")$estimate
  }
}
 # Result for point #1, i.e. different in Spearman's correlation when reported "across" or "within" samples for drug sensitivity data:
median(na.omit(corOutSpearman))
# [1] 0.6243905
median(na.omit(corOutPearson))
# [1] 0.8033401

## Nilotinib rank of BCR-ABL1 positive cell lines in CCLE and CGP overlapping cell lines:
rank(na.omit(auc.ccle[!is.na(auc.ccle[, "drugid_NILOTINIB"]) & !is.na(auc.cgp[, "drugid_NILOTINIB"]), "drugid_NILOTINIB"]))[c("MEG-01", "EM-2")]
# MEG-01   EM-2 
#    187    188

rank(na.omit(auc.cgp[!is.na(auc.ccle[, "drugid_NILOTINIB"]) & !is.na(auc.cgp[, "drugid_NILOTINIB"]), "drugid_NILOTINIB"]))[c("MEG-01", "EM-2")]
# MEG-01   EM-2 
#    188    187

## Scatterplot nilotinib AUC in CGP (NB: Fig 1(a))
cgpNilotinibDat <- na.omit(auc.cgp[!is.na(auc.ccle[, "drugid_NILOTINIB"]) & !is.na(auc.cgp[, "drugid_NILOTINIB"]), "drugid_NILOTINIB"])
cgpNilotinibDat <- c(cgpNilotinibDat, 0.64) ## During manuscript review, Haibe-Kains et al. kindly informed us that, having re-curated the data, there was a 3rd BCR-ABL positive cell line (KU812) among those overlapping CCLE and CGP. We include the AUC value for this cell line here.
dataSort <- 1 - sort(cgpNilotinibDat, decreasing=TRUE)
seqData <- 1:length(dataSort)
df <- data.frame(dataSort, seqData)

pdf("fig1a_new_1.pdf", width=7, height=5)
ggplot(data=df, aes(x=seqData, y=dataSort)) + theme_bw() + geom_point(position = position_jitter(height = 0), size=2, alpha = I(0.5), color="#003333") + xlab("Cell lines ordered by Nilotinib AUC") + ylab("Nilotinib AUC")
dev.off()

## Calculate the probability of the three BCR-ABL1 positive cell lines being ranked as the 3 most sensitive cell lines in CCLE (P-value is reported in paper).
1 / choose(189, 3)
# [1] 9.030047e-07

## Calculate median AUC values for CGP and CCLE samples for each drug. Result reported in main text and Supplementary Table.
## In the cell lines compared, 9 drugs have median AUC of > 0.95 in CGP:
medianAucsCGP <- apply(auc.cgp, 2, function(vec)return(1 - median(na.omit(vec)))) # CGP
#  drugid_ERLOTINIB  drugid_LAPATINIB  drugid_PHA665752 drugid_CRIZOTINIB 
#          0.993255          0.991030          0.995560          0.988340 
#     drugid_TAE684  drugid_NILOTINIB    drugid_AZD0530  drugid_SORAFENIB 
#          0.890670          0.991560          0.987060          0.972170 
#  drugid_PD0332991    drugid_PLX4720  drugid_PD0325901    drugid_AZD6244 
#          0.894780          0.978875          0.899840          0.925260 
#    drugid_NUTLIN3      drugid_17AAG drugid_PACLITAXEL 
#          0.984035          0.774210          0.796730

## In the cell lines compared, 9 drugs have AUC > 0.9 in CCLE.
medianAucsCcle <- apply(auc.ccle, 2, function(vec)return(1 - median(na.omit(vec)))) # CCLE
#  drugid_ERLOTINIB  drugid_LAPATINIB  drugid_PHA665752 drugid_CRIZOTINIB 
#         0.9320375         0.9253750         0.9410250         0.9107500 
#     drugid_TAE684  drugid_NILOTINIB    drugid_AZD0530  drugid_SORAFENIB 
#         0.8326250         0.9315000         0.8863750         0.9392500 
#  drugid_PD0332991    drugid_PLX4720  drugid_PD0325901    drugid_AZD6244 
#         0.9287750         0.9457125         0.7733750         0.8598037 
#    drugid_NUTLIN3      drugid_17AAG drugid_PACLITAXEL 
#         0.9411250         0.5585000         0.3206250 


## For each drug, calculate the correlation between variance in CCLE/CGP and the correlation between the two studies.
## Unsurprisingly, In both cases there is a significant correlation between variability and Spearman correlation between the two studies (i.e. Haibe-Kains measure of "concordance").
## These results were used for figure 1(b)
drugVarCgp <- apply(auc.cgp, 2, function(vec)return(var(na.omit(vec))))

## Variance in CGP vs Concordance between CCLE and CGP
cor.test(drugVarCgp, correlations.stats$auc[, "rho"], method="spearman")
# 	Spearman's rank correlation rho
# 
# data:  drugVarCgp and correlations.stats$auc[, "rho"]
# S = 240, p-value = 0.02862
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5714286

## Variance in CCLE vs Concordance between CCLE and CGP
drugVarCcle <- apply(auc.ccle, 2, function(vec)return(var(na.omit(vec))))
cor.test(drugVarCcle, correlations.stats$auc[, "rho"], method="spearman")
# 	Spearman's rank correlation rho
# 
# data:  drugVarCcle and correlations.stats$auc[, "rho"]
# S = 96, p-value = 0.0001883
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.8285714 

## The same drugs are also tend to induce high levels of variability in both studies (result not reported in paper).
cor.test(drugVarCgp, drugVarCcle, method="spearman")
# 	Spearman's rank correlation rho
# 
# data:  drugVarCgp and drugVarCcle
# S = 178, p-value = 0.006525
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.6821429 

## Fig 1(b): a scatterplot of Variability in either dataset vs Concordance between the datasets.
drugNames <- rownames(correlations.stats$auc)

# make the drug names look better when they appear on the plot.
drugNames[1] <- "Erlotinib"
drugNames[2] <- "Lapatinib"
drugNames[4] <- "Crizotinib"
drugNames[6] <- "Nilotinib"
drugNames[8] <- "Sorafenib"
drugNames[13] <- "Nutlin-3"
drugNames[15] <- "Paclitaxel"

library(ggplot2)
dat <- data.frame(x=correlations.stats$auc[, "rho"], y=drugVarCcle, Drug=drugNames, varCpg=drugVarCgp)
pdf("fig1b_bw_1.pdf", width=7, height=5)
ggplot(data=dat, aes(x=y, y=x)) + theme_bw() + geom_point(aes(color=varCpg), size=I(3)) + geom_text(aes(label=Drug), vjust=-.5, hjust=-.24, size=2.5, angle=15) + ylab("Correlation of AUC between CCLE and CGP") + xlab("Variance of AUC in CCLE") + scale_color_continuous(low="steelblue4",high="tomato2", name=expression(sigma^2 ~ "CPG")) + theme(legend.position=c(.9,.2))
dev.off()

## The Spearman correlation tests below show the relationship between lack of drug response (i.e. high median AUC) and lack of variability.
## (Result not included in paper)
cor.test(drugVarCgp, medianAucsCGP, method="spearman") # CGP
#         Spearman's rank correlation rho
# 
# data:  drugVarCgp and medianAucsCGP
# S = 1060, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.8928571

cor.test(drugVarCcle, medianAucsCcle, method="spearman") # CCLE
#         Spearman's rank correlation rho
# 
# data:  drugVarCcle and medianAucsCcle
# S = 1026, p-value = 0.0001578
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.8321429


## Here we show that the drug response data used in the Haibe-Kains et al Figure 3 break the 
## assumptions of a linear model fit, thus rendering the results presented in Figure 3 highly questionable.
## This point was removed from the final version of the paper during review due to word count constraints.
## 1. Do the shapiro-wilk tests on the distribution of the AUC data. 
## 2. Create QQplots.
# first prepare the data
gdscDrugName <- c("Erlotinib", "PHA-665752", "PF-02341066", "NVP-TAE684", "Nilotinib", "AZD-0530", "Sorafenib", "PD-0332991", "PLX4720", "PD-0325901", "WO2009093972", "Nutlin-3a", "17-AAG", "Paclitaxel", "Lapatinib")
drugAucFile <- read.csv("gdsc_manova_input_w2.csv", as.is=TRUE) # Drug data from CGP (ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/current_release/gdsc_manova_input_w2.csv)
sensList <- list()
sensListNas <- list()
sensBigVec <- numeric()
colnames(drugAucFile)[colnames(drugAucFile) == "X17.AAG_AUC"] <- "17.AAG_AUC"
drugVecBig <- character()
for(i in 1:15)
{
  sensDrug <- as.numeric(drugAucFile[, paste(gsub("-", ".", gdscDrugName[i]), "_AUC", sep="")])
  names(sensDrug) <- drugAucFile[, 1]
  drugVecBig <- c(drugVecBig, rep(gdscDrugName[i], length(sensDrug)))
  sensBigVec <- c(sensBigVec, ((sensDrug - mean(na.omit(sensDrug))) / sd(na.omit(sensDrug))) )
  sensList[[i]] <- na.omit(sensDrug)
  sensListNas[[i]] <- sensDrug
  print(i)
}

# do a shapiro test on all auc data
pOutShapirosCGP <- numeric()
for(i in 1:length(sensList))
{
  pOutShapirosCGP[i] <- shapiro.test(sensList[[i]] / sd(sensList[[i]]))$p.value
}

# All shaprio-wilk results are highly significant.
print(pOutShapirosCGP)
#  [1] 3.048542e-27 2.620912e-33 6.318809e-31 4.571068e-15 3.162312e-40
#  [6] 2.850502e-26 1.653762e-23 1.813438e-13 1.583423e-31 3.162116e-20
# [11] 3.598103e-13 4.543510e-29 3.886809e-11 3.928300e-10 1.402829e-29

# calculate the maximum p-value from there (this is reported in point 2 of the paper...)
print(max(pOutShapirosCGP))
#[1] 3.9283e-10

# QQplots for these AUCs.
par(mfrow=c(3,5))
for(i in 1:15)
{
  qqnorm(sensList[[i]], main=gdscDrugName[i])
  qqline(sensList[[i]], col="red")
}

####################################################################################
##################### End Geeleher et. al Code #####################################
####################################################################################




## AUCs per tissue types
correlations.stats.auc.tissue <- array(NA, dim=c(length(tissuen), length(drugn), 5), dimnames=list(tissuen, gsub("drugid_", "", drugn), c("rho", "upper", "lower", "p", "n")))
for(ii in 1:length(tissuen)) {
  for(i in 1:ncol(auc.ccle)) {
    iix <- !is.na(tissue) & tissue == tissuen[ii]
    nnn <- sum(complete.cases(auc.ccle[iix, i], auc.cgp[iix, i]))
    if(nnn >= minsample) {
      cc <- cor.test(auc.ccle[iix, i], auc.cgp[iix, i], method="spearman", use="complete.obs", alternative="greater")
    } else {
      cc <- list("estimate"=NA, "p.value"=NA)
    }
    ## correlation statistics
    cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
    correlations.stats.auc.tissue[ii, i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
  }
}

## boxplot for auc
pdf(file.path(saveres, sprintf("boxplot_auc_ccle_cgp_tissue.pdf")), height=10, width=14)
## correlations for each drug per tissue type
tt <- correlations.stats.auc.tissue[ , , "rho"]
## remove the non significant p-values
# tt[correlations.stats.auc.tissue[ , , "p"] >= 0.05] <- NA
ll <- c(list("all_tissues"=correlations.stats[["auc"]][ , "rho"]), lapply(apply(tt, 1, list), function(x) { return(x[[1]]) }))
ll <- lapply(ll, function(x) { return(x[!is.na(x)]) })
# ll <- ll[sapply(ll, length) > 0]
myx <- which(sapply(ll, length) < 3)
wt <- kruskal.test(x=ll[sapply(ll, length) > 0])
pp <- NULL
for(j in myx) { for(jj in ll[[j]]) { pp <- rbind(pp, c(j, jj)) } }
ll[sapply(ll, length) < 3] <- NA
par(las=2, mar=c(12, 4, 4, 2) + 0.1, xaxt="n")
mp <- boxplot(ll, outline=FALSE, ylab="Rs", main="Correlations of drug sensitivity (AUC) across tissue types\nCCLE vs. CGP", border="black", ylim=c(-1, 1), col=c("lightgrey", coltissue))
points(x=pp[ , 1], y=pp[ , 2], col=coltissue[pp[ , 1]], pch=20)
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.02), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", coltissuet))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
dev.off()

## boxplot for paper
pdf(file.path(saveres, sprintf("boxplot2_auc_ccle_cgp_tissue_paper.pdf")), height=8, width=14)
## correlations for each drug per tissue type
tt <- correlations.stats.auc.tissue[ , , "rho"]
## remove the non significant p-values
# tt[correlations.stats.auc.tissue[ , , "p"] >= 0.05] <- NA
ll <- c(list("all_tissues"=correlations.stats[["auc"]][ , "rho"]), lapply(apply(tt, 1, list), function(x) { return(x[[1]]) }))
ll <- lapply(ll, function(x) { return(x[!is.na(x)]) })
# ll <- ll[sapply(ll, length) > 0]
myx <- which(sapply(ll, length) < 3)
wt <- kruskal.test(x=ll[sapply(ll, length) > 0])
ll[sapply(ll, length) < 3] <- NA
par(las=1, mar=c(12, 4, 4, 2) + 0.1, xaxt="n")
mp <- boxplot(ll, outline=FALSE, ylab="Rs", main="Correlations of drug sensitivity (AUC) across tissue types\nCCLE vs. CGP", border="grey50", ylim=c(-1, 1), col="white", range=0)
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=c("lightgrey", coltissuet))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=1:length(ll))
points(x=jitter(unlist(xx), 1), y=unlist(ll), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", coltissuet))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
dev.off()

## test whether correlations are significantly higher in specific tissue types compared to all tissues combined
cor.diff.test <- NULL
ll <- c(list("all_tissues"=correlations.stats[["auc"]][ , "rho"]), lapply(apply(correlations.stats.auc.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ddnn <- names(ll[[1]])
for(i in 2:length(ll)) {
  if(sum(complete.cases(ll[[1]], y=ll[[i]])) < 3) {
   wt <- list("estimate"=NA, "p.value"=NA) 
  } else {
    wt <- wilcox.test(x=ll[[1]], y=ll[[i]], alternative="less", conf.int=TRUE)
  }
  cor.diff.test <- rbind(cor.diff.test, c(wt$estimate, wt$p.value))
}
dimnames(cor.diff.test) <- list(names(ll)[-1], c("difference", "p"))
# cor.diff.test <- cbind(cor.diff.test, "lfdr"=fdrtool::fdrtool(cor.diff.test[ , 2], statistic="pvalue", verbose=FALSE)$lfdr)
cor.diff.test <- cbind(cor.diff.test, "fdr"=p.adjust(cor.diff.test[ , "p"], method="fdr"))
myx <- complete.cases(cor.diff.test) & (cor.diff.test[ , "p"] < 0.05 | cor.diff.test[ , "fdr"] < myfdr)
tt <- cor.diff.test[myx, , drop=FALSE]
tt <- tt[order(tt[ , "p"], decreasing=FALSE), , drop=FALSE]
message("Are correlations between AUCs in some tissue types significantly higher than in all tissues combined?")
print(tt)

## discretized AUCs (sensitivity calling)

## To quantitatively assess concordance of drug sensitvity calling between CGP and CCLE, we used Cohen Kappa coefficient (k) (67), as implemented in the R package epibasix ; k ranges from 0 to 1, with 0 indicating no relation and 1 indicating a perfect concordance. Typically qualitative descriptions are associated with intervals [k ≤ 0.20, slight agreement; 0.20 < k ≤ 0.40, fair agreement; 0.40 < k ≤ 0.60, moderate agreement, 0.60 < k ≤ 0.80, substantial agreement; and 0.80 < k ≤ 0.99, almost perfect agreement, as described in Weigelt et al. (22)].

message("Correlation AUC sensitivity calls:")

pdf(file.path(saveres, "auc_sensitivity_calling_ccle_cgp.pdf"), width=9, height=11)
par(mfrow=c(5, 3), mar=c(8, 4, 4, 2) + 0.1)
mykap <- NULL
for (i in 1:ncol(auc.call.ccle)) {
  ff1 <- factor(x=auc.call.ccle[ , i], levels=c("resistant", "intermediate", "sensitive"))
  ff2 <- factor(x=auc.call.cgp[ , i], levels=c("resistant", "intermediate", "sensitive"))
  levels(ff1)[2] <- NA
  levels(ff2)[2] <- NA
  tt <- table("CCLE"=ff1, "CGP"=ff2)
  tts <- vcd::assocstats(tt)
  ## check if all elements are on the diagonal
  ttt <- tt
  diag(ttt) <- 0
  if(sum(ttt) == 0) {
    rr <- rr.l <- rr.u <- 1
  } else {
    err <- try(rr <- epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if(class(err) == "try-error") { rr <- rr.l <- rr.u <- NA }
    rr.l <- rr$CIL
    rr.u <- rr$CIU
    rr <- rr$kappa
  }
  mykap <- c(mykap, rr)
  plot.new()
  plotrix::addtable2plot(x=par("usr")[1] + 0.15, y=par("usr")[4] - 0.5, xjust=0, yjust=0, table=as.matrix(tt), bty="o", display.rownames=TRUE, display.colnames=TRUE, hlines=FALSE, vlines=FALSE, title="CCLE   vs  CGP", cex=1.25)
  title(main=sprintf("AUC sensitivity calling\n%s", gsub("drugid_", "", colnames(auc.call.ccle)[i])), sub=sprintf("Kappa=%.2g, 95%%CI [%.2g,%.2g], p=%.1E", rr, rr.l, rr.u, tts$chisq_tests[1,3]), cex=1)
}
names(mykap) <- colnames(auc.call.ccle)
dev.off()
## with intermediate
pdf(file.path(saveres, "auc_sensitivity_calling_intermediate_ccle_cgp_paper.pdf"), width=8, height=9)
par(mfrow=c(5, 3), mar=c(5, 4, 4, 2) + 0.1)
mykap <- NULL
for (i in 1:ncol(auc.call.ccle)) {
  ff1 <- factor(x=auc.call.ccle[ , i], levels=c("resistant", "intermediate", "sensitive"))
  ff2 <- factor(x=auc.call.cgp[ , i], levels=c("resistant", "intermediate", "sensitive"))
  levels(ff1) <- c("res", "inter", "sens")
  levels(ff2) <- c("res", "inter", "sens")
  nnn <- sum(complete.cases(ff1, ff2))
  tt <- table("CCLE"=ff1, "CGP"=ff2)
  tts <- vcd::assocstats(tt)
  ## check if all elements are on the diagonal
  ttt <- tt
  diag(ttt) <- 0
  if(sum(ttt) == 0) {
    rr <- rr.l <- rr.u <- 1
  } else {
    err <- try(rr <- epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if(class(err) == "try-error") { rr <- rr.l <- rr.u <- NA }
    rr.l <- rr$CIL
    rr.u <- rr$CIU
    rr <- rr$kappa
  }
  mykap <- c(mykap, rr)
  plot.new()
  plotrix::addtable2plot(x=par("usr")[1] + 0.1, y=par("usr")[4] - 0.4, xjust=0, yjust=0, table=as.matrix(tt), bty="o", display.rownames=TRUE, display.colnames=TRUE, hlines=FALSE, vlines=FALSE, title="CCLE   vs   CGP", cex=1.25)
  title(main=sprintf("AUC sensitivity calling\n%s", gsub("drugid_", "", colnames(ic50.call.ccle)[i])), sub=sprintf("Kappa=%.2g, 95%%CI [%.2g,%.2g], p=%.1E", rr, rr.l, rr.u, tts$chisq_tests[1,3]))
  correlations[["auc.call"]][i, "drug.sensitivity"] <- rr
  ## statistics
  correlations.stats[["auc.call"]][i, ] <- c(rr, rr.l, rr.u, tts$chisq_tests[1,3], nnn)
}
names(mykap) <- colnames(ic50.call.ccle)
dev.off()
## histogram of Kappa coefficients
pdf(file.path(saveres, "auc_sensitivity_calling_kappa_ccle_cgp.pdf"), width=7, height=7)
hist(mykap, main=sprintf("Sensitivity calling based on AUC\nKappa coefficients (CCLE vs CGP)", xlab="Kappa"))
dev.off()


## save all correlations statistics
savecor <- NULL
for(i in 1:length(correlations.stats)) {
  tt <- correlations.stats[[i]][ , apply(correlations.stats[[i]], 2, function(x) { return(!all(is.na(x))) }), drop=FALSE]
  savecor <- c(savecor, list(data.frame(tt)))
}
names(savecor) <- names(correlations.stats)
WriteXLS::WriteXLS("savecor", ExcelFileName=file.path(saveres, "correlations_stats.xls"), row.names=TRUE)

save(list=c("correlations", "correlations.stats"), compress=TRUE, file=file.path(saveres, "correlations_stats.RData"))
save(list=c("correlations.stats.ic50.tissue", "correlations.stats.auc.tissue"), compress=TRUE, file=file.path(saveres, "correlations_stats_tissue.RData"))

## barplots for paper
pdf(file.path(saveres, "barplot_ic50_auc_cgp_ccle_paper.pdf"), height=6, width=6)
par(las=3, mar=c(6, 4, 1, 0), xaxt="n")
## IC50
xx <- correlations.stats[["ic50"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.stats[["ic50"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.stats[["ic50"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.stats[["ic50"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.stats[["ic50"]])
cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
## AUC
xx <- correlations.stats[["auc"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.stats[["auc"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.stats[["auc"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.stats[["auc"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.stats[["auc"]])
cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main="Drug sensitivity measures")
legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
# plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
text(x=mp[1, ] + 1.55, y=rbind(cors.ic50$rho, cors.auc$rho)[1, ], pos=2, labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
text(x=mp[2, ] + 1.55, y=rbind(cors.ic50$rho, cors.auc$rho)[2, ], pos=2, labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()

## barplot for paper
pdf(file.path(saveres, "barplot_ic50_auc_tissue_cgp_ccle.pdf"), height=6, width=6)
for(ii in 1:length(tissuen)) {
  ncelline <- sum(tissue == tissuen[ii])
  if(max(correlations.stats.ic50.tissue[ii, , "n"], na.rm=TRUE) >= minsample) {
    ## IC50
    xx <- correlations.stats.ic50.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.stats.ic50.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.stats.ic50.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.stats.ic50.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.stats.ic50.tissue)[[2]]
    cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    ## AUC
    xx <- correlations.stats.auc.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.stats.auc.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.stats.auc.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.stats.auc.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.stats.auc.tissue)[[2]]
    cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    par(las=1, mar=c(6, 4, 6, 0), xaxt="n")
    # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
    yylim <- c(0, 1)
    mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main=sprintf("Drug sensitivity measures\n%s (%i cell lines)", gsub("_", " ", tissuen[ii]), ncelline))
    # legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
    # plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
    axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
    text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
    text(x=mp[1, ], y=cors.ic50$rho, pos=ifelse(is.na(cors.ic50$rho) | cors.ic50$rho < 0, 1, 3), labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
    text(x=mp[2, ], y=cors.auc$rho, pos=ifelse(is.na(cors.auc$rho) | cors.auc$rho < 0, 1, 3), labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
  }
}
dev.off()

## multiple barplots per page
count <- 0
nf <- 6
for(ii in 1:length(tissuen)) {
  ncelline <- sum(tissue == tissuen[ii])
  if(max(correlations.stats.ic50.tissue[ii, , "n"], na.rm=TRUE) >= minsample) {
    count <- count + 1
    ## IC50
    xx <- correlations.stats.ic50.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.stats.ic50.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.stats.ic50.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.stats.ic50.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.stats.ic50.tissue)[[2]]
    cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    ## AUC
    xx <- correlations.stats.auc.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.stats.auc.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.stats.auc.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.stats.auc.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.stats.auc.tissue)[[2]]
    cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    if((count %% nf) == 1) {
      pdf(file.path(saveres, sprintf("barplot_ic50_auc_tissue%i_cgp_ccle_paper.pdf", ceiling(count / nf))), height=14, width=10)
      par(las=1, mar=c(6, 4, 6, 0), xaxt="n", mfrow=c(3, 2))
    }
    # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
    yylim <- c(0, 1)
    mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main=sprintf("Drug sensitivity measures\n%s (%i cell lines)", gsub("_", " ", tissuen[ii]), ncelline))
    legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
    # plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
    axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
    text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
    text(x=mp[1, ], y=cors.ic50$rho, pos=ifelse(is.na(cors.ic50$rho) | cors.ic50$rho < 0, 1, 3), labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
    text(x=mp[2, ], y=cors.auc$rho, pos=ifelse(is.na(cors.auc$rho) | cors.auc$rho < 0, 1, 3), labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
    if((count %% nf) == 0) { dev.off() }
  }
}
if((count %% nf) != 0) { dev.off() }

## IC50
pdf(file.path(saveres, "barplot_ic50_cgp_ccle.pdf"), height=6, width=6)
par(las=3, mar=c(6, 4, 3, 0), xaxt="n")
xx <- correlations.stats[["ic50"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.stats[["ic50"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.stats[["ic50"]][ , "upper"]
ll[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.stats[["ic50"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.stats[["ic50"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
title(main="IC50\nCGP vs. CCLE")
# plotrix::plotCI(x=mp, y=xx, li=ll, ui=uu, err="y", pch=".", add=TRUE)
dev.off()
## with scatter plot
pdf(file.path(saveres, "scatterbarplot_ic50_cgp_ccle_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
for(i in 1:ncol(ic50.ccle)) {
  xxlim <- round(range(ic50.cgp[, i], na.rm=TRUE) * 10) / 10
  yylim <- round(range(ic50.ccle[, i], na.rm=TRUE) * 10) / 10
  par(mar=c(4, 4, 3, 1) + 0.1)
  myScatterPlot(x=ic50.cgp[, i], y=ic50.ccle[, i], xlab=ifelse(i > 11, "-log10 IC50 (CGP)", ""), ylab=ifelse((i %% 4) == 1, "-log10 IC50 (CCLE)", ""), main=gsub("drugid_", "", colnames(ic50.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75)
}
par(mar=c(6, 4, 3, 1) + 0.1, xaxt="n", las=1)
xx <- correlations.stats[["ic50"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.stats[["ic50"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.stats[["ic50"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.stats[["ic50"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.stats[["ic50"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
# plotrix::plotCI(x=mp, y=xx, li=ll, ui=uu, err="y", pch=".", add=TRUE)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()

## IC50 within tested drug concentrations
ic50.ccle.t <- ic50.ccle
ic50.ccle.t[!ic50.filt.ccle] <- NA
ic50.cgp.t <- ic50.cgp
ic50.cgp.t[!ic50.filt.cgp] <- NA
pdf(file.path(saveres, "barplot_ic50_cgp_ccle_filtconc.pdf"), height=6, width=6)
ccall <- NULL
for(i in 1:ncol(ic50.ccle.t)) {
  xxlim <- round(range(ic50.cgp.t[, i], na.rm=TRUE) * 10) / 10
  yylim <- round(range(ic50.ccle.t[, i], na.rm=TRUE) * 10) / 10
  nnn <- sum(complete.cases(ic50.ccle.t[, i], ic50.cgp.t[, i]))
  if(nnn >= 3) {
    cc <- cor.test(ic50.ccle.t[, i], ic50.cgp.t[, i], method="spearman", use="complete.obs", alternative="greater")
    cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
    ccall <- rbind(ccall, c("rho"=cc$estimate, "lower"=cci[1], "upper"=cci[2], "p"=cc$p.value, "n"=nnn))
    correlations[["ic50"]][i, "drug.sensitivity.filt"] <- cc$estimate
    ## correlation statistics
    cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
    correlations.stats[["ic50.filt"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
  } else {
    ccall <- rbind(ccall, c("rho"=NA, "lower"=NA, "upper"=NA, "p"=NA, "n"=nnn))
  }
}
dimnames(ccall) <- list(gsub("drugid_", "", colnames(ic50.ccle.t)), c("rho", "lower", "upper", "p", "n"))
par(mar=c(6, 4, 3, 1) + 0.1, xaxt="n", las=1)
xx <- ccall[ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- ccall[ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- ccall[ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- ccall[ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(ccall)
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
# plotrix::plotCI(x=mp, y=xx, li=ll, ui=uu, err="y", pch=".", add=TRUE)
dev.off()
## with scatterplot
pdf(file.path(saveres, "scatterbarplot_ic50_cgp_ccle_filtconc_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
ccall <- NULL
for(i in 1:ncol(ic50.ccle)) { 
  xxlim <- round(range(ic50.cgp.t[, i], na.rm=TRUE) * 10) / 10
  yylim <- round(range(ic50.ccle.t[, i], na.rm=TRUE) * 10) / 10
  nnn <- sum(complete.cases(ic50.ccle.t[, i], ic50.cgp.t[, i]))
  if(nnn >= minsample) {
    cc <- cor.test(ic50.ccle.t[, i], ic50.cgp.t[, i], method="spearman", use="complete.obs", alternative="greater")
    cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
    ccall <- rbind(ccall, c("rho"=cc$estimate, "lower"=cci[1], "upper"=cci[2], "p"=cc$p.value, "n"=nnn))
  } else {
    ccall <- rbind(ccall, c("rho"=NA, "lower"=NA, "upper"=NA, "p"=NA, "n"=nnn))
  }
  par(mar=c(4, 4, 3, 1) + 0.1)
  myScatterPlot(x=ic50.cgp.t[, i], y=ic50.ccle.t[, i], xlab=ifelse(i > 11, "-log10 IC50 (CGP)", ""), ylab=ifelse((i %% 4) == 1, "-log10 IC50 (CCLE)", ""), main=gsub("drugid_", "", colnames(ic50.ccle.t)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75)
}
dimnames(ccall) <- list(gsub("drugid_", "", colnames(ic50.ccle.t)), c("rho", "lower", "upper", "p", "n"))
par(mar=c(6, 4, 3, 1) + 0.1, xaxt="n", las=1)
xx <- ccall[ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- ccall[ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- ccall[ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- ccall[ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(ccall)
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()
# print(cbind(correlations[["ic50"]][ ,1], ccall[ ,1]))
# print(wilcox.test(correlations[["ic50"]][ ,1], ccall[ ,1], paired=TRUE))

## AUC
pdf(file.path(saveres, "barplot_auc_cgp_ccle.pdf"), height=6, width=6)
par(las=3, mar=c(6, 4, 3, 0), xaxt="n")
xx <- correlations.stats[["auc"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.stats[["auc"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.stats[["auc"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.stats[["auc"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.stats[["auc"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
title(main="AUC\nCGP vs. CCLE")
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()
## with scatter plot
pdf(file.path(saveres, "scatterbarplot_auc_cgp_ccle_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
for(i in 1:ncol(auc.ccle)) {
  xxlim <- round(range(auc.cgp[, i], na.rm=TRUE) * 10) / 10
  yylim <- round(range(auc.ccle[, i], na.rm=TRUE) * 10) / 10
  ccix <- complete.cases(auc.ccle[ ,i], auc.cgp[ ,i])
  nnn <- sum(ccix)
  par(mar=c(6, 4, 3, 1) + 0.1)
  myScatterPlot(x=auc.cgp[ ,i], y=auc.ccle[ ,i], xlab=ifelse(i > 11, "AUC (CGP)", ""), ylab=ifelse((i %% 4) == 1, "AUC (CCLE)", ""), main=gsub("drugid_", "", colnames(auc.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75)
}
par(mar=c(6, 4, 3, 1) + 0.1, xaxt="n", las=1)
xx <- correlations.stats[["auc"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.stats[["auc"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.stats[["auc"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.stats[["auc"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.stats[["auc"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()

## Kappa
## IC50
pdf(file.path(saveres, "barplot_ic50_call_cgp_ccle.pdf"), height=6, width=6)
par(las=3, mar=c(7, 4, 3, 0), xt="n")
xx <- correlations.stats[["ic50.call"]][ , "rho"]
ll <- correlations.stats[["ic50.call"]][ , "lower"]
ll[!is.na(ll) & ll < -1] <- -1
uu <- correlations.stats[["ic50.call"]][ , "upper"]
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.stats[["ic50.call"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.stats[["ic50.call"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.3)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Kappa", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp+1.25, y=par("usr")[3] - 0.03, pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
title(main="IC50 calls\nCGP vs. CCLE")
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()

## AUC
pdf(file.path(saveres, "barplot_auc_call_cgp_ccle.pdf"), height=6, width=6)
par(las=3, mar=c(6, 4, 3, 0), xaxt="n")
xx <- correlations.stats[["auc.call"]][ , "rho"]
ll <- correlations.stats[["auc.call"]][ , "lower"]
ll[!is.na(ll) & ll < -1] <- -1
uu <- correlations.stats[["auc.call"]][ , "upper"]
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.stats[["auc.call"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.stats[["auc"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.5)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Kappa", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp+1.25, y=par("usr")[3] - 0.01, pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
title(main="AUC calls\nCGP vs. CCLE")
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()

pdf(file.path(saveres, "barplot_ic50_auc_call_cgp_ccle.pdf"), height=12, width=6)
par(mfrow=c(2, 1))
## IC50
par(las=3, mar=c(6, 4, 1, 0), xaxt="n")
xx <- correlations.stats[["ic50.call"]][ , "rho"]
ll <- correlations.stats[["ic50.call"]][ , "lower"]
ll[!is.na(ll) & ll < -1] <- -1
uu <- correlations.stats[["ic50.call"]][ , "upper"]
pp <- correlations.stats[["ic50.call"]][ , "p"]
uu[!is.na(uu) & uu > 1] <- 1
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.stats[["ic50.call"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(-0.05, 0.5)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Kappa (IC50 calls)", ylim=yylim)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
## AUC
par(las=3, mar=c(6, 4, 1, 0), xaxt="n")
xx <- correlations.stats[["auc.call"]][ , "rho"]
ll <- correlations.stats[["auc.call"]][ , "lower"]
ll[!is.na(ll) & ll < -1] <- -1
uu <- correlations.stats[["auc.call"]][ , "upper"]
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.stats[["auc.call"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.stats[["auc.call"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(-0.05, 0.5)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Kappa (AUC calls)", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp+1.25, y=par("usr")[3], pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()


###############################################################################
###############################################################################
## compute associations between genes and drugs
###############################################################################
###############################################################################

## correlation statistics
tt <- matrix(NA, nrow=length(drugsn), ncol=5, dimnames=list(drugsn, c("rho", "lower", "upper", "p", "n")))
correlations.assoc.stats <- list("ic50"=tt, "ic50.filt"=tt, "ic50.call"=tt, "ic50.call.filt"=tt, "auc"=tt, "auc.filt"=tt, "auc.call"=tt, "auc.call.filt"=tt)

########################
## IC50 per tissue type

for(ii in 1:length(tissuen)) {
  myfn <- file.path(saveres, sprintf("assoc_ic50_cgp_ccle_%s.RData", toupper(gsub(badchars, "", tissuen[ii]))))
  if(!file.exists(myfn)) {
    ## CCLE
    message(sprintf("Gene-drug association based on IC50 for tissue %s (CCLE)", tissuen[ii]))
    assoc.ic50.ccle <- NULL
    for(i in 1:ncol(ic50.ccle)) {
      message("Computation for drug ", gsub("drugid_", "", colnames(ic50.ccle)[i]))
      myx <- !is.na(tissue.ccle) & tissue.ccle == tissuen[ii]
      splitix <- parallel::splitIndices(nx=ncol(data.ccle), ncl=nbcore)
      mcres <- parallel::mclapply(splitix, function(x, data, ic50, tissue, tissuex) {
        res <- apply(X=data[ , x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=ic50, z=tissue, method=genedrugm)
        return(res)      
      }, data=data.ccle[myx, , drop=FALSE], ic50=ic50.ccle[myx, i], tissue=tissue.ccle[myx])
      mcres <- t(do.call(cbind, mcres))
      mcres <- mcres[colnames(data.ccle), , drop=FALSE]
      mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
      assoc.ic50.ccle <- c(assoc.ic50.ccle, list(mcres))
    }
    message("")
    names(assoc.ic50.ccle) <- colnames(ic50.ccle)
    ## CGP
    message(sprintf("Gene-drug association based on IC50 for tissue %s (CGP)", tissuen[ii]))
    assoc.ic50.cgp <- NULL
    for(i in 1:ncol(ic50.cgp)) {
      message("Computation for drug ", gsub("drugid_", "", colnames(ic50.cgp)[i]))
      myx <- !is.na(tissue.cgp) & tissue.cgp == tissuen[ii]
      splitix <- parallel::splitIndices(nx=ncol(data.cgp), ncl=nbcore)
      mcres <- parallel::mclapply(splitix, function(x, data, ic50, tissue) {
        res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=ic50, z=tissue, method=genedrugm)
        return(res)      
      }, data=data.cgp[myx, , drop=FALSE], ic50=ic50.cgp[myx, i], tissue=tissue.cgp[myx])
      mcres <- t(do.call(cbind, mcres))
      mcres <- mcres[colnames(data.cgp), , drop=FALSE]
      mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
      assoc.ic50.cgp <- c(assoc.ic50.cgp, list(mcres))
    }
    message("")
    names(assoc.ic50.cgp) <- colnames(ic50.cgp)
    ## save all associations
    ## CCLE
    rr <- NULL
    for(i in 1:length(assoc.ic50.ccle)) {
      tt <- cbind(assoc.ic50.ccle[[i]], "entrez_gene_id"=annot[rownames(assoc.ic50.ccle[[i]]), "jetset.EntrezID"], "gene_symbol"=annot[rownames(assoc.ic50.ccle[[i]]), "jetset.symbol"])
      rr <- c(rr, list(data.frame(tt)))
    }
    names(rr) <- names(assoc.ic50.ccle)
    WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, sprintf("ic50_ccle_results_gene_drug_%s.xls", toupper(gsub(badchars, "", tissuen[ii])))), row.names=TRUE)
    ## CGP
    rr <- NULL
    for(i in 1:length(assoc.ic50.cgp)) {
      tt <- cbind(assoc.ic50.cgp[[i]], "entrez_gene_id"=annot[rownames(assoc.ic50.cgp[[i]]), "jetset.EntrezID"], "gene_symbol"=annot[rownames(assoc.ic50.cgp[[i]]), "jetset.symbol"])
      rr <- c(rr, list(data.frame(tt)))
    }
    names(rr) <- names(assoc.ic50.cgp)
    WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, sprintf("ic50_cgp_results_gene_drug_%s.xls", toupper(gsub(badchars, "", tissuen[ii])))), row.names=TRUE)
    save(list=c("assoc.ic50.cgp", "assoc.ic50.ccle"), compress=TRUE, file=myfn)
  }
}

correlations.assoc.stats.ic50.tissue <- array(NA, dim=c(length(tissuen), length(drugn), 5), dimnames=list(tissuen, gsub("drugid_", "", drugn), c("rho", "upper", "lower", "p", "n")))
for(ii in 1:length(tissuen)) {
  myfn <- file.path(saveres, sprintf("assoc_ic50_cgp_ccle_%s.RData", toupper(gsub(badchars, "", tissuen[ii]))))
  if(file.exists(myfn)) {
    load(myfn)
    for(i in 1:length(assoc.ic50.ccle)) {
      nnn <- sum(complete.cases(assoc.ic50.ccle[[i]][ , "estimate"], assoc.ic50.cgp[[i]][ , "estimate"]))
      if(nnn >= minsample) {
        cc <- cor.test(assoc.ic50.ccle[[i]][ , "estimate"], assoc.ic50.cgp[[i]][ , "estimate"], method="spearman", use="complete.obs", alternative="greater")
      } else {
        cc <- list("estimate"=NA, "p.value"=NA)
      }
      ## correlation statistics
      cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
      correlations.assoc.stats.ic50.tissue[ii, i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
    }
  }
}

correlations.assoc.stats.ic50.signif.tissue <- array(NA, dim=c(length(tissuen), length(drugn), 5), dimnames=list(tissuen, gsub("drugid_", "", drugn), c("rho", "upper", "lower", "p", "n")))
for(ii in 1:length(tissuen)) {
  myfn <- file.path(saveres, sprintf("assoc_ic50_cgp_ccle_%s.RData", toupper(gsub(badchars, "", tissuen[ii]))))
  if(file.exists(myfn)) {
    load(myfn)
    for(i in 1:length(assoc.ic50.ccle)) {
      myx <- assoc.ic50.ccle[[i]][ , "fdr"] < myfdr | assoc.ic50.cgp[[i]][ , "fdr"] < myfdr
      nnn <- sum(complete.cases(assoc.ic50.ccle[[i]][myx, "estimate"], assoc.ic50.cgp[[i]][myx, "estimate"]))
      if(nnn >= minsample) {
        cc <- cor.test(assoc.ic50.ccle[[i]][myx, "estimate"], assoc.ic50.cgp[[i]][myx, "estimate"], method="spearman", use="complete.obs", alternative="greater")
      } else {
        cc <- list("estimate"=NA, "p.value"=NA)
      }
      ## correlation statistics
      cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
      correlations.assoc.stats.ic50.signif.tissue[ii, i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
    }
  }
}

## IC50 all tissues
myfn <- file.path(saveres, "assoc_ic50_cgp_ccle.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on IC50 (CCLE)")
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
  message("Gene-drug association based on IC50 (CGP)")
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
  ## save all associations
  ## CCLE
  rr <- NULL
  for(i in 1:length(assoc.ic50.ccle)) {
    tt <- cbind(assoc.ic50.ccle[[i]], "entrez_gene_id"=annot[rownames(assoc.ic50.ccle[[i]]), "jetset.EntrezID"], "gene_symbol"=annot[rownames(assoc.ic50.ccle[[i]]), "jetset.symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.ic50.ccle)
  WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "ic50_ccle_results_gene_drug_paper.xls"), row.names=TRUE)
  ## CGP
  rr <- NULL
  for(i in 1:length(assoc.ic50.cgp)) {
    tt <- cbind(assoc.ic50.cgp[[i]], "entrez_gene_id"=annot[rownames(assoc.ic50.cgp[[i]]), "jetset.EntrezID"], "gene_symbol"=annot[rownames(assoc.ic50.cgp[[i]]), "jetset.symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.ic50.cgp)
  WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "ic50_cgp_results_gene_drug_paper.xls"), row.names=TRUE)
  save(list=c("assoc.ic50.cgp", "assoc.ic50.ccle"), compress=TRUE, file=myfn)
} else { load(myfn) }

## consistency between associations
## all genes
## IC50
message("Correlation gene-drug associations based on IC50, all genes:")
pdf(file.path(saveres, "assoc_ic50_ccle_cgp_all.pdf"), height=10, width=15)
par(mfrow=c(3, 5))
for(i in 1:length(assoc.ic50.ccle)) {
  tt <- min(round(quantile(abs(c(assoc.ic50.cgp[[i]][ , "estimate"], assoc.ic50.ccle[[i]][ , "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 5)
  xxlim <- yylim <- c(-tt, tt)
  myx <- complete.cases(assoc.ic50.ccle[[i]][ ,"estimate"], assoc.ic50.cgp[[i]][ ,"estimate"])
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(assoc.ic50.cgp[[i]][ ,"estimate"], assoc.ic50.ccle[[i]][ ,"estimate"], method="spearman", use="complete.obs", alternative="greater")
  }
  # message("SCC for ", names(assoc.ic50.cgp)[i], ": ", cc$estimate)
  # plot(x=assoc.ic50.cgp[[i]][ ,"estimate"], y=assoc.ic50.ccle[[i]][ ,"estimate"], xlab="Gene-drug association IC50 (CGP)", ylab="Gene-drug association IC50 (CCLE)", main=gsub("drugid_", "", names(assoc.ic50.ccle)[i]), pch=16, col=rgb(0, 200, 50, 10, maxColorValue=255), xlim=c(0.2, 0.8), ylim=c(0.2, 0.8))
  ## smoothScatter
  # ramp <- colorRamp(c("white", "yellow", "orange", "red"))
  # ramp <- colorRamp(c("white", "green", "turquoise", "blue"))
  # myramp <- rgb(ramp(seq(0, 1, length=10)), max=255)
  # mypalette <- colorRampPalette(myramp)
  myScatterPlot(x=assoc.ic50.cgp[[i]][ ,"estimate"], y=assoc.ic50.ccle[[i]][ ,"estimate"], xlab="Gene-drug association IC50 (CGP)", ylab="Gene-drug association IC50 (CCLE)", main=gsub("drugid_", "", names(assoc.ic50.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.1, yjust=0.8, bty="n", legend=sprintf("R=%.3g, p=%.1E", cc$estimate, cc$p.value), text.font=2)
  correlations[["ic50"]][i, "gene.drug"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=sum(myx), alpha=0.05)
  correlations.assoc.stats[["ic50"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, sum(myx))
}
message("")
dev.off()
## genes with significant gene-drug association in at least one dataset (FDR)
message("Correlation gene-drug associations based on IC50, significant genes:")
pdf(file.path(saveres, "assoc_ic50_ccle_cgp_fdr.pdf"), height=10, width=15)
par(mfrow=c(3, 5))
for(i in 1:length(assoc.ic50.ccle)) {
  tt <- min(round(quantile(abs(c(assoc.ic50.cgp[[i]][ , "estimate"], assoc.ic50.ccle[[i]][ , "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 5)
  xxlim <- yylim <- c(-tt, tt)
  myx <- complete.cases(assoc.ic50.ccle[[i]][ ,"fdr"], assoc.ic50.cgp[[i]][ ,"fdr"]) & (assoc.ic50.ccle[[i]][ ,"fdr"] < myfdr | assoc.ic50.cgp[[i]][ ,"fdr"] < myfdr)
  if(sum(myx) < 3) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(assoc.ic50.cgp[[i]][myx, "estimate"], assoc.ic50.ccle[[i]][myx, "estimate"], method="spearman", use="complete.obs", alternative="greater")
    tt <- min(round(quantile(abs(c(assoc.ic50.cgp[[i]][myx, "estimate"], assoc.ic50.ccle[[i]][myx, "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 3)
    xxlim <- yylim <- c(-tt, tt)  
  }
  # message("SCC for ", names(assoc.ic50.cgp)[i], ": ", cc$estimate)
  myScatterPlot(assoc.ic50.cgp[[i]][myx, "estimate"], y=assoc.ic50.ccle[[i]][myx, "estimate"], xlab="Gene-drug association IC50 (CGP)", ylab="Gene-drug association IC50 (CCLE)", main=gsub("drugid_", "", names(assoc.ic50.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.1, yjust=0.8, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, sum(myx)), text.font=2)
  correlations[["ic50"]][i, "gene.drug.filt"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=sum(myx), alpha=0.05)
  correlations.assoc.stats[["ic50.filt"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, sum(myx))
}
message("")
dev.off()

## boxplot for gene-drug associations based on ic50
pdf(file.path(saveres, sprintf("boxplot2_assoc_ic50_ccle_cgp_tissue_paper.pdf")), height=8, width=14)
## correlations for each drug per tissue type
ll <- c(list("all_tissues"=correlations.assoc.stats[["ic50"]][ , "rho"]), lapply(apply(correlations.assoc.stats.ic50.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ll <- lapply(ll, function(x) { return(x[!is.na(x)]) })
# ll <- ll[sapply(ll, length) > 0]
myx <- which(sapply(ll, length) < 3)
wt <- kruskal.test(x=ll[sapply(ll, length) > 0])
ll[sapply(ll, length) < 3] <- NA
par(las=1, mar=c(12, 4, 4, 2) + 0.1, xaxt="n")
mp <- boxplot(ll, outline=FALSE, ylab="Rs", main="Correlations of gene-drug associations (IC50) across tissue types\nCCLE vs. CGP", border="grey50", ylim=c(-1, 1), col="white", range=0)
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=c("lightgrey", coltissuet))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=1:length(ll))
points(x=jitter(unlist(xx), 1), y=unlist(ll), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", coltissuet))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
dev.off()

## test whether correlations are significantly higher in specific tissue types compared to all tissues combined
cor.diff.test <- NULL
ll <- c(list("all_tissues"=correlations.assoc.stats[["ic50"]][ , "rho"]), lapply(apply(correlations.assoc.stats.ic50.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ddnn <- names(ll[[1]])
for(i in 2:length(ll)) {
  if(sum(complete.cases(ll[[1]], y=ll[[i]])) < 3) {
   wt <- list("estimate"=NA, "p.value"=NA) 
  } else {
    wt <- wilcox.test(x=ll[[1]], y=ll[[i]], alternative="less", conf.int=TRUE)
  }
  cor.diff.test <- rbind(cor.diff.test, c(wt$estimate, wt$p.value))
}
dimnames(cor.diff.test) <- list(names(ll)[-1], c("difference", "p"))
# cor.diff.test <- cbind(cor.diff.test, "lfdr"=fdrtool::fdrtool(cor.diff.test[ , 2], statistic="pvalue", verbose=FALSE)$lfdr)
cor.diff.test <- cbind(cor.diff.test, "fdr"=p.adjust(cor.diff.test[ , "p"], method="fdr"))
myx <- complete.cases(cor.diff.test) & (cor.diff.test[ , "p"] < 0.05 | cor.diff.test[ , "fdr"] < myfdr)
tt <- cor.diff.test[myx, , drop=FALSE]
tt <- tt[order(tt[ , "p"], decreasing=FALSE), , drop=FALSE]
message("Are correlations between gene-drug associations (IC50) in some tissue types\n significantly higher than in all tissues combined?")
if(nrow(tt) == 0) { message("\t-> none") } else { print(tt) }

## boxplot for significant gene-drug associations based on ic50
pdf(file.path(saveres, sprintf("boxplot2_assoc_ic50_ccle_cgp_signif_tissue_paper.pdf")), height=8, width=14)
## correlations for each drug per tissue type
ll <- c(list("all_tissues"=correlations.assoc.stats[["ic50.filt"]][ , "rho"]), lapply(apply(correlations.assoc.stats.ic50.signif.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ll <- lapply(ll, function(x) { return(x[!is.na(x)]) })
# ll <- ll[sapply(ll, length) > 0]
myx <- which(sapply(ll, length) < 3)
wt <- kruskal.test(x=ll[sapply(ll, length) > 0])
ll[sapply(ll, length) < 3] <- NA
par(las=1, mar=c(12, 4, 4, 2) + 0.1, xaxt="n")
mp <- boxplot(ll, outline=FALSE, ylab="Rs", main=sprintf("Correlations of gene-drug associations (IC50) across tissue types (FDR < %i%%)\nCCLE vs. CGP", round(myfdr * 100)), border="grey50", ylim=c(-1, 1), col="white", range=0)
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=c("lightgrey", coltissuet))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=1:length(ll))
points(x=jitter(unlist(xx), 1), y=unlist(ll), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", coltissuet))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
dev.off()

## test whether correlations are significantly higher in specific tissue types compared to all tissues combined
cor.diff.test <- NULL
ll <- c(list("all_tissues"=correlations.assoc.stats[["ic50.filt"]][ , "rho"]), lapply(apply(correlations.assoc.stats.ic50.signif.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ddnn <- names(ll[[1]])
for(i in 2:length(ll)) {
  if(sum(complete.cases(ll[[1]], y=ll[[i]])) < 3) {
   wt <- list("estimate"=NA, "p.value"=NA) 
  } else {
    wt <- wilcox.test(x=ll[[1]], y=ll[[i]], alternative="less", conf.int=TRUE)
  }
  cor.diff.test <- rbind(cor.diff.test, c(wt$estimate, wt$p.value))
}
dimnames(cor.diff.test) <- list(names(ll)[-1], c("difference", "p"))
# cor.diff.test <- cbind(cor.diff.test, "lfdr"=fdrtool::fdrtool(cor.diff.test[ , 2], statistic="pvalue", verbose=FALSE)$lfdr)
cor.diff.test <- cbind(cor.diff.test, "fdr"=p.adjust(cor.diff.test[ , "p"], method="fdr"))
myx <- complete.cases(cor.diff.test) & (cor.diff.test[ , "p"] < 0.05 | cor.diff.test[ , "fdr"] < myfdr)
tt <- cor.diff.test[myx, , drop=FALSE]
tt <- tt[order(tt[ , "p"], decreasing=FALSE), , drop=FALSE]
message("Are correlations between significant gene-drug associations (IC50) in some tissue types\n significantly higher than in all tissues combined?")
if(nrow(tt) == 0) { message("\t-> none") } else { print(tt) }


########################
## IC50 sensitivity calling

myfn <- file.path(saveres, "assoc_ic50_call_cgp_ccle.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on IC50 sensitivity calling (CCLE)")
  assoc.ic50.call.ccle <- NULL
  for(i in 1:ncol(ic50.call.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(ic50.call.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, ic50, tissue) {
      ic50 <- factor(ic50, levels=c("resistant", "intermediate", "sensitive"))
      ## remove the intermediate
      levels(ic50)[2] <- NA
      tissue <- factor(tissue)
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=ic50, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.ccle, ic50=ic50.call.ccle[ ,i], tissue=tissue.ccle)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ccle), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.ic50.call.ccle <- c(assoc.ic50.call.ccle, list(mcres))
  }
  message("")
  names(assoc.ic50.call.ccle) <- colnames(ic50.call.ccle)
  ## CGP
  message("Gene-drug association based on IC50 sensitivity calling (CGP)")
  assoc.ic50.call.cgp <- NULL
  for(i in 1:ncol(ic50.call.cgp)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(ic50.call.cgp)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.cgp), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, ic50, tissue) {
      ic50 <- factor(ic50, levels=c("resistant", "intermediate", "sensitive"))
      ## remove the intermediate
      levels(ic50)[2] <- NA
      tissue <- factor(tissue)
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=ic50, z=tissue, method=genedrugm)
      return(res) 
    }, data=data.cgp, ic50=ic50.call.cgp[ ,i], tissue=tissue.cgp)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.cgp), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.ic50.call.cgp <- c(assoc.ic50.call.cgp, list(mcres))
  }
  message("")
  names(assoc.ic50.call.cgp) <- colnames(ic50.call.cgp)
  ## save all associations
  ## CCLE
  rr <- NULL
  for(i in 1:length(assoc.ic50.call.ccle)) {
    tt <- cbind(assoc.ic50.call.ccle[[i]], "entrez_gene_id"=annot[rownames(assoc.ic50.call.ccle[[i]]), "jetset.EntrezID"], "gene_symbol"=annot[rownames(assoc.ic50.call.ccle[[i]]), "jetset.symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.ic50.call.ccle)
  WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "ic50_call_ccle_results_gene_drug.xls"), row.names=TRUE)
  ## CGP
  rr <- NULL
  for(i in 1:length(assoc.ic50.call.cgp)) {
    tt <- cbind(assoc.ic50.call.cgp[[i]], "entrez_gene_id"=annot[rownames(assoc.ic50.call.cgp[[i]]), "jetset.EntrezID"], "gene_symbol"=annot[rownames(assoc.ic50.call.cgp[[i]]), "jetset.symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.ic50.call.cgp)
  WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "ic50_call_cgp_results_gene_drug.xls"), row.names=TRUE)
  save(list=c("assoc.ic50.call.cgp", "assoc.ic50.call.ccle"), compress=TRUE, file=myfn)
} else { load(myfn) }

## consistency between associations
## all genes
message("Correlation gene-drug associations based on IC50 sensitivity calling, all genes:")
pdf(file.path(saveres, "assoc_ic50_call_ccle_cgp_all.pdf"), height=10, width=15)
par(mfrow=c(3, 5))
for(i in 1:length(assoc.ic50.call.ccle)) {
  tt <- min(round(quantile(abs(c(assoc.ic50.call.cgp[[i]][ , "estimate"], assoc.ic50.call.ccle[[i]][ , "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 10)
  xxlim <- yylim <- c(-tt, tt)
  myx <- complete.cases(assoc.ic50.call.ccle[[i]][ ,"estimate"], assoc.ic50.call.cgp[[i]][ ,"estimate"])
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(assoc.ic50.call.cgp[[i]][ ,"estimate"], assoc.ic50.call.ccle[[i]][ ,"estimate"], method="spearman", use="complete.obs", alternative="greater")
  }
  # message("SCC for ", names(assoc.ic50.call.cgp)[i], ": ", cc$estimate)
  myScatterPlot(x=assoc.ic50.call.cgp[[i]][ ,"estimate"], y=assoc.ic50.call.ccle[[i]][ ,"estimate"], xlab="Gene-drug association IC50 calling (CGP)", ylab="Gene-drug association IC50 calling (CCLE)", main=gsub("drugid_", "", names(assoc.ic50.call.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.1, yjust=0.8, bty="n", legend=sprintf("R=%.3g, p=%.1E", cc$estimate, cc$p.value), text.font=2)
  correlations[["ic50.call"]][i, "gene.drug"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=sum(myx), alpha=0.05)
  correlations.assoc.stats[["ic50.call"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, sum(myx))
}
message("")
dev.off()
## genes with significant gene-drug association in at least one dataset (FDR)
message("Correlation gene-drug associations based on IC50 sensitivity calling, significant genes:")
pdf(file.path(saveres, "assoc_ic50_call_ccle_cgp_fdr.pdf"), height=10, width=15)
par(mfrow=c(3, 5))
for(i in 1:length(assoc.ic50.call.ccle)) {
  tt <- min(round(quantile(abs(c(assoc.ic50.call.cgp[[i]][ , "estimate"], assoc.ic50.call.ccle[[i]][ , "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 10)
  xxlim <- yylim <- c(-tt, tt)
  myx <- complete.cases(assoc.ic50.call.ccle[[i]][ ,"fdr"], assoc.ic50.call.cgp[[i]][ ,"fdr"]) & (assoc.ic50.call.ccle[[i]][ ,"fdr"] < myfdr | assoc.ic50.call.cgp[[i]][ ,"fdr"] < myfdr)
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(assoc.ic50.call.cgp[[i]][myx, "estimate"], assoc.ic50.call.ccle[[i]][myx, "estimate"], method="spearman", use="complete.obs", alternative="greater")
    tt <- min(round(quantile(abs(c(assoc.ic50.call.cgp[[i]][myx, "estimate"], assoc.ic50.call.ccle[[i]][myx, "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 10)
    xxlim <- yylim <- c(-tt, tt)
  }
  # message("SCC for ", names(assoc.ic50.call.cgp)[i], ": ", cc$estimate)
  myScatterPlot(x=assoc.ic50.call.cgp[[i]][myx, "estimate"], y=assoc.ic50.call.ccle[[i]][myx, "estimate"], xlab="Gene-drug association IC50 calling (CGP)", ylab="Gene-drug association IC50 calling (CCLE)", main=gsub("drugid_", "", names(assoc.ic50.call.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.1, yjust=0.8, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, sum(myx)), text.font=2)
  correlations[["ic50.call"]][i, "gene.drug.filt"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=sum(myx), alpha=0.05)
  correlations.assoc.stats[["ic50.call.filt"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, sum(myx))
}
message("")
dev.off()

########################
## AUC per tissue types

for(ii in 1:length(tissuen)) {
  myfn <- file.path(saveres, sprintf("assoc_auc_cgp_ccle_%s.RData", toupper(gsub(badchars, "", tissuen[ii]))))
  if(!file.exists(myfn)) {
    ## CCLE
    message(sprintf("Gene-drug association based on AUC for tissue %s (CCLE)", tissuen[ii]))
    assoc.auc.ccle <- NULL
    for(i in 1:ncol(auc.ccle)) {
      message("Computation for drug ", gsub("drugid_", "", colnames(auc.ccle)[i]))
      myx <- !is.na(tissue.ccle) & tissue.ccle == tissuen[ii]
      splitix <- parallel::splitIndices(nx=ncol(data.ccle), ncl=nbcore)
      mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue, tissuex) {
        res <- apply(X=data[ , x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm)
        return(res)      
      }, data=data.ccle[myx, , drop=FALSE], auc=auc.ccle[myx, i], tissue=tissue.ccle[myx])
      mcres <- t(do.call(cbind, mcres))
      mcres <- mcres[colnames(data.ccle), , drop=FALSE]
      mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
      assoc.auc.ccle <- c(assoc.auc.ccle, list(mcres))
    }
    message("")
    names(assoc.auc.ccle) <- colnames(auc.ccle)
    ## CGP
    message(sprintf("Gene-drug association based on AUC for tissue %s (CGP)", tissuen[ii]))
    assoc.auc.cgp <- NULL
    for(i in 1:ncol(auc.cgp)) {
      message("Computation for drug ", gsub("drugid_", "", colnames(auc.cgp)[i]))
      myx <- !is.na(tissue.cgp) & tissue.cgp == tissuen[ii]
      splitix <- parallel::splitIndices(nx=ncol(data.cgp), ncl=nbcore)
      mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
        res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm)
        return(res)      
      }, data=data.cgp[myx, , drop=FALSE], auc=auc.cgp[myx, i], tissue=tissue.cgp[myx])
      mcres <- t(do.call(cbind, mcres))
      mcres <- mcres[colnames(data.cgp), , drop=FALSE]
      mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
      assoc.auc.cgp <- c(assoc.auc.cgp, list(mcres))
    }
    message("")
    names(assoc.auc.cgp) <- colnames(auc.cgp)
    ## save all associations
    ## CCLE
    rr <- NULL
    for(i in 1:length(assoc.auc.ccle)) {
      tt <- cbind(assoc.auc.ccle[[i]], "entrez_gene_id"=annot[rownames(assoc.auc.ccle[[i]]), "jetset.EntrezID"], "gene_symbol"=annot[rownames(assoc.auc.ccle[[i]]), "jetset.symbol"])
      rr <- c(rr, list(data.frame(tt)))
    }
    names(rr) <- names(assoc.auc.ccle)
    WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, sprintf("auc_ccle_results_gene_drug_%s.xls", toupper(gsub(badchars, "", tissuen[ii])))), row.names=TRUE)
    ## CGP
    rr <- NULL
    for(i in 1:length(assoc.auc.cgp)) {
      tt <- cbind(assoc.auc.cgp[[i]], "entrez_gene_id"=annot[rownames(assoc.auc.cgp[[i]]), "jetset.EntrezID"], "gene_symbol"=annot[rownames(assoc.auc.cgp[[i]]), "jetset.symbol"])
      rr <- c(rr, list(data.frame(tt)))
    }
    names(rr) <- names(assoc.auc.cgp)
    WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, sprintf("auc_cgp_results_gene_drug_%s.xls", toupper(gsub(badchars, "", tissuen[ii])))), row.names=TRUE)
    save(list=c("assoc.auc.cgp", "assoc.auc.ccle"), compress=TRUE, file=myfn)
  }
}

correlations.assoc.stats.auc.tissue <- array(NA, dim=c(length(tissuen), length(drugn), 5), dimnames=list(tissuen, gsub("drugid_", "", drugn), c("rho", "upper", "lower", "p", "n")))
for(ii in 1:length(tissuen)) {
  myfn <- file.path(saveres, sprintf("assoc_auc_cgp_ccle_%s.RData", toupper(gsub(badchars, "", tissuen[ii]))))
  if(file.exists(myfn)) {
    load(myfn)
    for(i in 1:length(assoc.auc.ccle)) {
      nnn <- sum(complete.cases(assoc.auc.ccle[[i]][ , "estimate"], assoc.auc.cgp[[i]][ , "estimate"]))
      if(nnn >= minsample) {
        cc <- cor.test(assoc.auc.ccle[[i]][ , "estimate"], assoc.auc.cgp[[i]][ , "estimate"], method="spearman", use="complete.obs", alternative="greater")
      } else {
        cc <- list("estimate"=NA, "p.value"=NA)
      }
      ## correlation statistics
      cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
      correlations.assoc.stats.auc.tissue[ii, i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
    }
  }
}

correlations.assoc.stats.auc.signif.tissue <- array(NA, dim=c(length(tissuen), length(drugn), 5), dimnames=list(tissuen, gsub("drugid_", "", drugn), c("rho", "upper", "lower", "p", "n")))
for(ii in 1:length(tissuen)) {
  myfn <- file.path(saveres, sprintf("assoc_auc_cgp_ccle_%s.RData", toupper(gsub(badchars, "", tissuen[ii]))))
  if(file.exists(myfn)) {
    load(myfn)
    for(i in 1:length(assoc.auc.ccle)) {
      myx <- assoc.auc.ccle[[i]][ , "fdr"] < myfdr | assoc.auc.cgp[[i]][ , "fdr"] < myfdr
      nnn <- sum(complete.cases(assoc.auc.ccle[[i]][myx, "estimate"], assoc.auc.cgp[[i]][myx, "estimate"]))
      if(nnn >= minsample) {
        cc <- cor.test(assoc.auc.ccle[[i]][myx, "estimate"], assoc.auc.cgp[[i]][myx, "estimate"], method="spearman", use="complete.obs", alternative="greater")
      } else {
        cc <- list("estimate"=NA, "p.value"=NA)
      }
      ## correlation statistics
      cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
      correlations.assoc.stats.auc.signif.tissue[ii, i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
    }
  }
}

## AUC for all tissues
myfn <- file.path(saveres, "assoc_auc_cgp_ccle.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on AUC (CCLE)")
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
  message("Gene-drug association based on AUC (CGP)")
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
  ## save all associations
  ## CCLE
  rr <- NULL
  for(i in 1:length(assoc.auc.ccle)) {
    tt <- cbind(assoc.auc.ccle[[i]], "entrez_gene_id"=annot[rownames(assoc.auc.ccle[[i]]), "jetset.EntrezID"], "gene_symbol"=annot[rownames(assoc.auc.ccle[[i]]), "jetset.symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.auc.ccle)
  WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "auc_ccle_results_gene_drug_paper.xls"), row.names=TRUE)
  ## CGP
  rr <- NULL
  for(i in 1:length(assoc.auc.cgp)) {
    tt <- cbind(assoc.auc.cgp[[i]], "entrez_gene_id"=annot[rownames(assoc.auc.cgp[[i]]), "jetset.EntrezID"], "gene_symbol"=annot[rownames(assoc.auc.cgp[[i]]), "jetset.symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.auc.cgp)
  WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "auc_cgp_results_gene_drug_paper.xls"), row.names=TRUE)
  save(list=c("assoc.auc.cgp", "assoc.auc.ccle"), compress=TRUE, file=myfn)
} else { load(myfn) }

## consistency between associations
## all genes
## AUC
message("Correlation gene-drug associations based on AUC, all genes:")
pdf(file.path(saveres, "assoc_auc_ccle_cgp_all.pdf"), height=10, width=15)
par(mfrow=c(3, 5))
for(i in 1:length(assoc.auc.ccle)) {
  tt <- min(round(quantile(abs(c(assoc.auc.cgp[[i]][ , "estimate"], assoc.auc.ccle[[i]][ , "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 5)
  xxlim <- yylim <- c(-tt, tt)
  myx <- complete.cases(assoc.auc.ccle[[i]][ ,"estimate"], assoc.auc.cgp[[i]][ ,"estimate"])
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(assoc.auc.cgp[[i]][ ,"estimate"], assoc.auc.ccle[[i]][ ,"estimate"], method="spearman", use="complete.obs", alternative="greater")
  }
  # message("SCC for ", names(assoc.auc.cgp)[i], ": ", cc$estimate)
  # plot(x=assoc.auc.cgp[[i]][ ,"estimate"], y=assoc.auc.ccle[[i]][ ,"estimate"], xlab="Gene-drug association AUC (CGP)", ylab="Gene-drug association AUC (CCLE)", main=gsub("drugid_", "", names(assoc.auc.ccle)[i]), pch=16, col=rgb(0, 200, 50, 10, maxColorValue=255), xlim=c(0.2, 0.8), ylim=c(0.2, 0.8))
  ## smoothScatter
  # ramp <- colorRamp(c("white", "yellow", "orange", "red"))
  # ramp <- colorRamp(c("white", "green", "turquoise", "blue"))
  # myramp <- rgb(ramp(seq(0, 1, length=10)), max=255)
  # mypalette <- colorRampPalette(myramp)
  myScatterPlot(x=assoc.auc.cgp[[i]][ ,"estimate"], y=assoc.auc.ccle[[i]][ ,"estimate"], xlab="Gene-drug association AUC (CGP)", ylab="Gene-drug association AUC (CCLE)", main=gsub("drugid_", "", names(assoc.auc.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.1, yjust=0.8, bty="n", legend=sprintf("R=%.3g, p=%.1E", cc$estimate, cc$p.value), text.font=2)
  correlations[["auc"]][i, "gene.drug"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=sum(myx), alpha=0.05)
  correlations.assoc.stats[["auc"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, sum(myx))
}
message("")
dev.off()
## genes with significant gene-drug association in at least one dataset (FDR)
message("Correlation gene-drug associations based on AUC, significant genes:")
pdf(file.path(saveres, "assoc_auc_ccle_cgp_fdr.pdf"), height=10, width=15)
par(mfrow=c(3, 5))
for(i in 1:length(assoc.auc.ccle)) {
  tt <- min(round(quantile(abs(c(assoc.auc.cgp[[i]][ , "estimate"], assoc.auc.ccle[[i]][ , "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 5)
  xxlim <- yylim <- c(-tt, tt)
  myx <- complete.cases(assoc.auc.ccle[[i]][ ,"fdr"], assoc.auc.cgp[[i]][ ,"fdr"]) & (assoc.auc.ccle[[i]][ ,"fdr"] < myfdr | assoc.auc.cgp[[i]][ ,"fdr"] < myfdr)
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(assoc.auc.cgp[[i]][myx, "estimate"], assoc.auc.ccle[[i]][myx, "estimate"], method="spearman", use="complete.obs", alternative="greater")
    tt <- min(round(quantile(abs(c(assoc.auc.cgp[[i]][myx, "estimate"], assoc.auc.ccle[[i]][myx, "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 5)
    xxlim <- yylim <- c(-tt, tt)  
  }
  # message("SCC for ", names(assoc.auc.cgp)[i], ": ", cc$estimate)
  myScatterPlot(assoc.auc.cgp[[i]][myx, "estimate"], y=assoc.auc.ccle[[i]][myx, "estimate"], xlab="Gene-drug association AUC (CGP)", ylab="Gene-drug association AUC (CCLE)", main=gsub("drugid_", "", names(assoc.auc.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.1, yjust=0.8, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, sum(myx)), text.font=2)
  correlations[["auc"]][i, "gene.drug.filt"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=sum(myx), alpha=0.05)
  correlations.assoc.stats[["auc.filt"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, sum(myx))
}
message("")
dev.off()

save(list=c("correlations", "correlations.assoc.stats"), compress=TRUE, file=file.path(saveres, "correlations_assoc_stats.RData"))
save(list=c("correlations.assoc.stats.ic50.tissue", "correlations.assoc.stats.auc.tissue", "correlations.assoc.stats.ic50.signif.tissue", "correlations.assoc.stats.auc.signif.tissue"), compress=TRUE, file=file.path(saveres, "correlations_assoc_stats_tissue.RData"))

## boxplot for gene-drug associations based on auc
pdf(file.path(saveres, sprintf("boxplot2_assoc_auc_ccle_cgp_tissue_paper.pdf")), height=8, width=14)
## correlations for each drug per tissue type
ll <- c(list("all_tissues"=correlations.assoc.stats[["auc"]][ , "rho"]), lapply(apply(correlations.assoc.stats.auc.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ll <- lapply(ll, function(x) { return(x[!is.na(x)]) })
# ll <- ll[sapply(ll, length) > 0]
myx <- which(sapply(ll, length) < 3)
wt <- kruskal.test(x=ll[sapply(ll, length) > 0])
ll[sapply(ll, length) < 3] <- NA
par(las=1, mar=c(12, 4, 4, 2) + 0.1, xaxt="n")
mp <- boxplot(ll, outline=FALSE, ylab="Rs", main="Correlations of gene-drug associations (AUC) across tissue types\nCCLE vs. CGP", border="grey50", ylim=c(-1, 1), col="white", range=0)
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=c("lightgrey", coltissuet))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=1:length(ll))
points(x=jitter(unlist(xx), 1), y=unlist(ll), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", coltissuet))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
dev.off()

## test whether correlations are significantly higher in specific tissue types compared to all tissues combined
cor.diff.test <- NULL
ll <- c(list("all_tissues"=correlations.assoc.stats[["auc"]][ , "rho"]), lapply(apply(correlations.assoc.stats.auc.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ddnn <- names(ll[[1]])
for(i in 2:length(ll)) {
  if(sum(complete.cases(ll[[1]], y=ll[[i]])) < 3) {
   wt <- list("estimate"=NA, "p.value"=NA) 
  } else {
    wt <- wilcox.test(x=ll[[1]], y=ll[[i]], alternative="less", conf.int=TRUE)
  }
  cor.diff.test <- rbind(cor.diff.test, c(wt$estimate, wt$p.value))
}
dimnames(cor.diff.test) <- list(names(ll)[-1], c("difference", "p"))
# cor.diff.test <- cbind(cor.diff.test, "lfdr"=fdrtool::fdrtool(cor.diff.test[ , 2], statistic="pvalue", verbose=FALSE)$lfdr)
cor.diff.test <- cbind(cor.diff.test, "fdr"=p.adjust(cor.diff.test[ , "p"], method="fdr"))
myx <- complete.cases(cor.diff.test) & (cor.diff.test[ , "p"] < 0.05 | cor.diff.test[ , "fdr"] < myfdr)
tt <- cor.diff.test[myx, , drop=FALSE]
tt <- tt[order(tt[ , "p"], decreasing=FALSE), , drop=FALSE]
message("Are correlations between gene-drug associations (AUC) in some tissue types\n significantly higher than in all tissues combined?")
if(nrow(tt) == 0) { message("\t-> none") } else { print(tt) }

## boxplot for significant gene-drug associations based on auc
pdf(file.path(saveres, sprintf("boxplot2_assoc_auc_ccle_cgp_signif_tissue_paper.pdf")), height=8, width=14)
## correlations for each drug per tissue type
ll <- c(list("all_tissues"=correlations.assoc.stats[["auc.filt"]][ , "rho"]), lapply(apply(correlations.assoc.stats.auc.signif.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ll <- lapply(ll, function(x) { return(x[!is.na(x)]) })
# ll <- ll[sapply(ll, length) > 0]
myx <- which(sapply(ll, length) < 3)
wt <- kruskal.test(x=ll[sapply(ll, length) > 0])
ll[sapply(ll, length) < 3] <- NA
par(las=1, mar=c(12, 4, 4, 2) + 0.1, xaxt="n")
mp <- boxplot(ll, outline=FALSE, ylab="Rs", main=sprintf("Correlations of gene-drug associations (AUC) across tissue types (FDR < %i%%)\nCCLE vs. CGP", round(myfdr * 100)), border="grey50", ylim=c(-1, 1), col="white", range=0)
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=c("lightgrey", coltissuet))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=1:length(ll))
points(x=jitter(unlist(xx), 1), y=unlist(ll), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", coltissuet))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
dev.off()

## test whether correlations are significantly higher in specific tissue types compared to all tissues combined
cor.diff.test <- NULL
ll <- c(list("all_tissues"=correlations.assoc.stats[["auc.filt"]][ , "rho"]), lapply(apply(correlations.assoc.stats.auc.signif.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ddnn <- names(ll[[1]])
for(i in 2:length(ll)) {
  if(sum(complete.cases(ll[[1]], y=ll[[i]])) < 3) {
   wt <- list("estimate"=NA, "p.value"=NA) 
  } else {
    wt <- wilcox.test(x=ll[[1]], y=ll[[i]], alternative="less", conf.int=TRUE)
  }
  cor.diff.test <- rbind(cor.diff.test, c(wt$estimate, wt$p.value))
}
dimnames(cor.diff.test) <- list(names(ll)[-1], c("difference", "p"))
# cor.diff.test <- cbind(cor.diff.test, "lfdr"=fdrtool::fdrtool(cor.diff.test[ , 2], statistic="pvalue", verbose=FALSE)$lfdr)
cor.diff.test <- cbind(cor.diff.test, "fdr"=p.adjust(cor.diff.test[ , "p"], method="fdr"))
myx <- complete.cases(cor.diff.test) & (cor.diff.test[ , "p"] < 0.05 | cor.diff.test[ , "fdr"] < myfdr)
tt <- cor.diff.test[myx, , drop=FALSE]
tt <- tt[order(tt[ , "p"], decreasing=FALSE), , drop=FALSE]
message("Are correlations between significant gene-drug associations (AUC) in some tissue types\n significantly higher than in all tissues combined?")
if(nrow(tt) == 0) { message("\t-> none") } else { print(tt) }

########################
## AUC sensitivity calling

myfn <- file.path(saveres, "assoc_auc_call_cgp_ccle.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on AUC sensitivity calling (CCLE)")
  assoc.auc.call.ccle <- NULL
  for(i in 1:ncol(auc.call.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.call.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      auc <- factor(auc, levels=c("resistant", "intermediate", "sensitive"))
      ## remove the intermediate
      levels(auc)[2] <- NA
      tissue <- factor(tissue)
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm)
      return(res)      
    }, data=data.ccle, auc=auc.call.ccle[ ,i], tissue=tissue.ccle)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ccle), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.call.ccle <- c(assoc.auc.call.ccle, list(mcres))
  }
  message("")
  names(assoc.auc.call.ccle) <- colnames(auc.call.ccle)
  ## CGP
  message("Gene-drug association based on AUC sensitivity calling (CGP)")
  assoc.auc.call.cgp <- NULL
  for(i in 1:ncol(auc.call.cgp)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.call.cgp)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.cgp), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      auc <- factor(auc, levels=c("resistant", "intermediate", "sensitive"))
      ## remove the intermediate
      levels(auc)[2] <- NA
      tissue <- factor(tissue)
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm)
      return(res) 
    }, data=data.cgp, auc=auc.call.cgp[ ,i], tissue=tissue.cgp)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.cgp), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.call.cgp <- c(assoc.auc.call.cgp, list(mcres))
  }
  message("")
  names(assoc.auc.call.cgp) <- colnames(auc.call.cgp)
  ## save all associations
  ## CCLE
  rr <- NULL
  for(i in 1:length(assoc.auc.call.ccle)) {
    tt <- cbind(assoc.auc.call.ccle[[i]], "entrez_gene_id"=annot[rownames(assoc.auc.call.ccle[[i]]), "jetset.EntrezID"], "gene_symbol"=annot[rownames(assoc.auc.call.ccle[[i]]), "jetset.symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.auc.call.ccle)
  WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "auc_call_ccle_results_gene_drug.xls"), row.names=TRUE)
  ## CGP
  rr <- NULL
  for(i in 1:length(assoc.auc.call.cgp)) {
    tt <- cbind(assoc.auc.call.cgp[[i]], "entrez_gene_id"=annot[rownames(assoc.auc.call.cgp[[i]]), "jetset.EntrezID"], "gene_symbol"=annot[rownames(assoc.auc.call.cgp[[i]]), "jetset.symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.auc.call.cgp)
  WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "auc_call_cgp_results_gene_drug.xls"), row.names=TRUE)
  save(list=c("assoc.auc.call.cgp", "assoc.auc.call.ccle"), compress=TRUE, file=myfn)
} else { load(myfn) }

## consistency between associations
## all genes
message("Correlation gene-drug associations based on AUC sensitivity calling, all genes:")
pdf(file.path(saveres, "assoc_auc_call_ccle_cgp_all.pdf"), height=10, width=15)
par(mfrow=c(3, 5))
for(i in 1:length(assoc.auc.call.ccle)) {
  tt <- min(round(quantile(abs(c(assoc.auc.call.cgp[[i]][ , "estimate"], assoc.auc.call.ccle[[i]][ , "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 10)
  xxlim <- yylim <- c(-tt, tt)
  myx <- complete.cases(assoc.auc.call.ccle[[i]][ ,"estimate"], assoc.auc.call.cgp[[i]][ ,"estimate"])
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(assoc.auc.call.cgp[[i]][ ,"estimate"], assoc.auc.call.ccle[[i]][ ,"estimate"], method="spearman", use="complete.obs", alternative="greater")
  }
  # message("SCC for ", names(assoc.auc.call.cgp)[i], ": ", cc$estimate)
  myScatterPlot(x=assoc.auc.call.cgp[[i]][ ,"estimate"], y=assoc.auc.call.ccle[[i]][ ,"estimate"], xlab="Gene-drug association AUC calling (CGP)", ylab="Gene-drug association AUC calling (CCLE)", main=gsub("drugid_", "", names(assoc.auc.call.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.1, yjust=0.8, bty="n", legend=sprintf("R=%.3g, p=%.1E", cc$estimate, cc$p.value), text.font=2)
  correlations[["auc.call"]][i, "gene.drug"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=sum(myx), alpha=0.05)
  correlations.assoc.stats[["auc.call"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, sum(myx))
}
message("")
dev.off()
## genes with significant gene-drug association in at least one dataset (FDR)
message("Correlation gene-drug associations based on AUC sensitivity calling, significant genes:")
pdf(file.path(saveres, "assoc_auc_call_ccle_cgp_fdr.pdf"), height=10, width=15)
par(mfrow=c(3, 5))
for(i in 1:length(assoc.auc.call.ccle)) {
  tt <- min(round(quantile(abs(c(assoc.auc.call.cgp[[i]][ , "estimate"], assoc.auc.call.ccle[[i]][ , "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 10)
  xxlim <- yylim <- c(-tt, tt)
  myx <- complete.cases(assoc.auc.call.ccle[[i]][ ,"fdr"], assoc.auc.call.cgp[[i]][ ,"fdr"]) & (assoc.auc.call.ccle[[i]][ ,"fdr"] < myfdr | assoc.auc.call.cgp[[i]][ ,"fdr"] < myfdr)
  if(sum(myx) < minsample) {
     cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(assoc.auc.call.cgp[[i]][myx, "estimate"], assoc.auc.call.ccle[[i]][myx, "estimate"], method="spearman", use="complete.obs", alternative="greater")
    tt <- min(round(quantile(abs(c(assoc.auc.call.cgp[[i]][myx, "estimate"], assoc.auc.call.ccle[[i]][myx, "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 10)
    xxlim <- yylim <- c(-tt, tt)
  }
  # message("SCC for ", names(assoc.auc.call.cgp)[i], ": ", cc$estimate)
  myScatterPlot(x=assoc.auc.call.cgp[[i]][myx, "estimate"], y=assoc.auc.call.ccle[[i]][myx, "estimate"], xlab="Gene-drug association AUC calling (CGP)", ylab="Gene-drug association AUC calling (CCLE)", main=gsub("drugid_", "", names(assoc.auc.call.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.1, yjust=0.8, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, sum(myx)), text.font=2)
  correlations[["auc.call"]][i, "gene.drug.filt"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=sum(myx), alpha=0.05)
  correlations.assoc.stats[["auc.call.filt"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, sum(myx))
}
message("")
dev.off()


## save all correlations statistics
savecor <- NULL
for(i in 1:length(correlations.assoc.stats)) {
  tt <- correlations.assoc.stats[[i]][ , apply(correlations.assoc.stats[[i]], 2, function(x) { return(!all(is.na(x))) }), drop=FALSE]
  savecor <- c(savecor, list(data.frame(tt)))
}
names(savecor) <- names(correlations.assoc.stats)
WriteXLS::WriteXLS("savecor", ExcelFileName=file.path(saveres, "correlations_assoc_stats.xls"), row.names=TRUE)


## barpots for correlations
## Spearman
for (i in 1:length(correlations.assoc.stats)) {
  fnnn <- sprintf("barplot_assoc_%s_cgp_ccle.pdf", gsub(badchars, "_", names(correlations.assoc.stats)[i]))
  pdf(file.path(saveres, fnnn), height=6, width=6)
  par(las=3, mar=c(6, 4, 3, 0), xaxt="n")
  xx <- correlations.assoc.stats[[i]][ , "rho"]
  xx[!is.na(xx) & xx < 0] <- 0
  ll <- correlations.assoc.stats[[i]][ , "lower"]
  ll[!is.na(ll) & ll < 0] <- 0
  uu <- correlations.assoc.stats[[i]][ , "upper"]
  uu[!is.na(uu) & uu < 0] <- 0
  uu[!is.na(uu) & uu > 1] <- 1
  pp <- correlations.assoc.stats[[i]][ , "p"]
  names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.assoc.stats[[i]])
  # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
  yylim <- c(0, 0.71)
  mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
  axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
  text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
  title(main=sprintf("Gene-drug associations %s\nCGP vs. CCLE", toupper(gsub(badchars, " ", names(correlations.assoc.stats)[i]))))
  text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
  dev.off()
}

## gene-drug associations IC50 with scatter plot
pdf(file.path(saveres, "scatterbarplot_assoc_ic50_cgp_ccle_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
for(i in 1:length(assoc.ic50.ccle)) {
  tt <- min(round(quantile(abs(c(assoc.ic50.cgp[[i]][ , "estimate"], assoc.ic50.ccle[[i]][ , "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 5)
  xxlim <- yylim <- c(-tt, tt)
  myx <- complete.cases(assoc.ic50.ccle[[i]][ ,"estimate"], assoc.ic50.cgp[[i]][ ,"estimate"])
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(assoc.ic50.cgp[[i]][ ,"estimate"], assoc.ic50.ccle[[i]][ ,"estimate"], method="spearman", use="complete.obs", alternative="greater")
  }
  par(mar=c(6, 4, 3, 1) + 0.1, las=1)
  myScatterPlot(x=assoc.ic50.cgp[[i]][ ,"estimate"], y=assoc.ic50.ccle[[i]][ ,"estimate"], xlab=ifelse(i > 11, "Gene-drug association IC50 (CGP)", ""), ylab=ifelse((i %% 4) == 1, "Gene-drug association IC50 (CCLE)", ""), main=gsub("drugid_", "", names(assoc.ic50.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
}
par(mar=c(6, 4, 3, 1) + 0.1, xaxt="n", las=1)
xx <- correlations.assoc.stats[["ic50"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.assoc.stats[["ic50"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.assoc.stats[["ic50"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.assoc.stats[["ic50"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.assoc.stats[["ic50"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()

## gene-drug associations IC50 fdr with scatter plot
pdf(file.path(saveres, "scatterbarplot_assoc_ic50_filt_cgp_ccle_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
for(i in 1:length(assoc.ic50.ccle)) {
  tt <- min(round(quantile(abs(c(assoc.ic50.cgp[[i]][ , "estimate"], assoc.ic50.ccle[[i]][ , "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 5)
  xxlim <- yylim <- c(-tt, tt)
  myx <- complete.cases(assoc.ic50.ccle[[i]][ ,"fdr"], assoc.ic50.cgp[[i]][ ,"fdr"]) & (assoc.ic50.ccle[[i]][ ,"fdr"] < myfdr | assoc.ic50.cgp[[i]][ ,"fdr"] < myfdr)
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(assoc.ic50.cgp[[i]][myx, "estimate"], assoc.ic50.ccle[[i]][myx, "estimate"], method="spearman", use="complete.obs", alternative="greater")
    tt <- min(round(quantile(abs(c(assoc.ic50.cgp[[i]][myx, "estimate"], assoc.ic50.ccle[[i]][myx, "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 3)
    xxlim <- yylim <- c(-tt, tt)
  }
  par(mar=c(6, 4, 3, 1) + 0.1, las=1)
  myScatterPlot(x=assoc.ic50.cgp[[i]][myx, "estimate"], y=assoc.ic50.ccle[[i]][myx, "estimate"], xlab=ifelse(i > 11, "Gene-drug association IC50 (CGP)", ""), ylab=ifelse((i %% 4) == 1, "Gene-drug association IC50 (CCLE)", ""), main=gsub("drugid_", "", names(assoc.ic50.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
}
par(mar=c(6, 4, 3, 1) + 0.1, xaxt="n", las=1)
xx <- correlations.assoc.stats[["ic50.filt"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.assoc.stats[["ic50.filt"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.assoc.stats[["ic50.filt"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.assoc.stats[["ic50.filt"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.assoc.stats[["ic50.filt"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()

## gene-drug associations AUC with scatter plot
pdf(file.path(saveres, "scatterbarplot_assoc_auc_cgp_ccle_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
for(i in 1:length(assoc.auc.ccle)) {
  tt <- min(round(quantile(abs(c(assoc.auc.cgp[[i]][ , "estimate"], assoc.auc.ccle[[i]][ , "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 5)
  xxlim <- yylim <- c(-tt, tt)
  myx <- complete.cases(assoc.auc.ccle[[i]][ ,"estimate"], assoc.auc.cgp[[i]][ ,"estimate"])
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(assoc.auc.cgp[[i]][ ,"estimate"], assoc.auc.ccle[[i]][ ,"estimate"], method="spearman", use="complete.obs", alternative="greater")
  }
  par(mar=c(6, 4, 3, 1) + 0.1, las=1)
  myScatterPlot(x=assoc.auc.cgp[[i]][ ,"estimate"], y=assoc.auc.ccle[[i]][ ,"estimate"], xlab=ifelse(i > 11, "Gene-drug association AUC (CGP)", ""), ylab=ifelse((i %% 4) == 1, "Gene-drug association AUC (CCLE)", ""), main=gsub("drugid_", "", names(assoc.auc.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
}
par(mar=c(6, 4, 3, 1) + 0.1, xaxt="n", las=1)
xx <- correlations.assoc.stats[["auc"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.assoc.stats[["auc"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.assoc.stats[["auc"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.assoc.stats[["auc"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.assoc.stats[["auc"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()

## gene-drug associations AUC fdr with scatter plot
pdf(file.path(saveres, "scatterbarplot_assoc_auc_filt_cgp_ccle_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
for(i in 1:length(assoc.auc.ccle)) {
  tt <- min(round(quantile(abs(c(assoc.auc.cgp[[i]][ , "estimate"], assoc.auc.ccle[[i]][ , "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 5)
  xxlim <- yylim <- c(-tt, tt)
  myx <- complete.cases(assoc.auc.ccle[[i]][ ,"fdr"], assoc.auc.cgp[[i]][ ,"fdr"]) & (assoc.auc.ccle[[i]][ ,"fdr"] < myfdr | assoc.auc.cgp[[i]][ ,"fdr"] < myfdr)
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(assoc.auc.cgp[[i]][myx, "estimate"], assoc.auc.ccle[[i]][myx, "estimate"], method="spearman", use="complete.obs", alternative="greater")
    tt <- min(round(quantile(abs(c(assoc.auc.cgp[[i]][myx, "estimate"], assoc.auc.ccle[[i]][myx, "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 5)
    xxlim <- yylim <- c(-tt, tt)
  }
  par(mar=c(6, 4, 3, 1) + 0.1, las=1)
  myScatterPlot(x=assoc.auc.cgp[[i]][myx, "estimate"], y=assoc.auc.ccle[[i]][myx, "estimate"], xlab=ifelse(i > 11, "Gene-drug association AUC (CGP)", ""), ylab=ifelse((i %% 4) == 1, "Gene-drug association AUC (CCLE)", ""), main=gsub("drugid_", "", names(assoc.auc.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
}
par(mar=c(6, 4, 3, 1) + 0.1, xaxt="n", las=1)
xx <- correlations.assoc.stats[["auc.filt"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.assoc.stats[["auc.filt"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.assoc.stats[["auc.filt"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.assoc.stats[["auc.filt"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.assoc.stats[["auc.filt"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()

## barplot for paper
pdf(file.path(saveres, "barplot_assoc_ic50_auc_tissue_cgp_ccle.pdf"), height=6, width=6)
for(ii in 1:length(tissuen)) {
  ncelline <- sum(tissue == tissuen[ii])
  if(max(correlations.assoc.stats.ic50.tissue[ii, , "n"], na.rm=TRUE) >= minsample && max(correlations.assoc.stats.ic50.tissue[ii, , "rho"], na.rm=TRUE) > 0) {
    ## IC50
    xx <- correlations.assoc.stats.ic50.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.assoc.stats.ic50.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.assoc.stats.ic50.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.assoc.stats.ic50.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.assoc.stats.ic50.tissue)[[2]]
    cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    ## AUC
    xx <- correlations.assoc.stats.auc.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.assoc.stats.auc.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.assoc.stats.auc.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.assoc.stats.auc.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.assoc.stats.auc.tissue)[[2]]
    cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    par(las=1, mar=c(6, 4, 6, 0), xaxt="n")
    # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
    yylim <- c(0, 1)
    mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main=sprintf("Drug sensitivity measures (all)\n%s (%i cell lines)", gsub("_", " ", tissuen[ii]), ncelline))
    # legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
    # plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
    axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
    text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
    text(x=mp[1, ], y=cors.ic50$rho, pos=ifelse(is.na(cors.ic50$rho) | cors.ic50$rho < 0, 1, 3), labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
    text(x=mp[2, ], y=cors.auc$rho, pos=ifelse(is.na(cors.auc$rho) | cors.auc$rho < 0, 1, 3), labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
  }
}
dev.off()

## multiple barplots per page
count <- 0
nf <- 6
for(ii in 1:length(tissuen)) {
  ncelline <- sum(tissue == tissuen[ii])
  if(max(correlations.assoc.stats.ic50.tissue[ii, , "n"], na.rm=TRUE) >= minsample && max(correlations.assoc.stats.ic50.tissue[ii, , "rho"], na.rm=TRUE) > 0) {
    count <- count + 1
    ## IC50
    xx <- correlations.assoc.stats.ic50.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.assoc.stats.ic50.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.assoc.stats.ic50.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.assoc.stats.ic50.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.assoc.stats.ic50.tissue)[[2]]
    cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    ## AUC
    xx <- correlations.assoc.stats.auc.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.assoc.stats.auc.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.assoc.stats.auc.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.assoc.stats.auc.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.assoc.stats.auc.tissue)[[2]]
    cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    if((count %% nf) == 1) {
      pdf(file.path(saveres, sprintf("barplot_assoc_ic50_auc_tissue%i_cgp_ccle_paper.pdf", ceiling(count / nf))), height=14, width=10)
      par(las=1, mar=c(6, 4, 6, 0), xaxt="n", mfrow=c(3, 2))
    }
    # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
    yylim <- c(0, 1)
    mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main=sprintf("Gene-drug associations (all)\n%s (%i cell lines)", gsub("_", " ", tissuen[ii]), ncelline))
    legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
    # plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
    axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
    text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
    text(x=mp[1, ], y=cors.ic50$rho, pos=ifelse(is.na(cors.ic50$rho) | cors.ic50$rho < 0, 1, 3), labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
    text(x=mp[2, ], y=cors.auc$rho, pos=ifelse(is.na(cors.auc$rho) | cors.auc$rho < 0, 1, 3), labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
    if((count %% nf) == 0) { dev.off() }
  }
}
if((count %% nf) != 0) { dev.off() }

## barplot for paper
pdf(file.path(saveres, "barplot_assoc_ic50_auc_signif_tissue_cgp_ccle.pdf"), height=6, width=6)
for(ii in 1:length(tissuen)) {
  ncelline <- sum(tissue == tissuen[ii])
  if(max(correlations.assoc.stats.ic50.signif.tissue[ii, , "n"], na.rm=TRUE) >= minsample) {
    ## IC50
    xx <- correlations.assoc.stats.ic50.signif.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.assoc.stats.ic50.signif.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.assoc.stats.ic50.signif.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.assoc.stats.ic50.signif.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.assoc.stats.ic50.signif.tissue)[[2]]
    cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    ## AUC
    xx <- correlations.assoc.stats.auc.signif.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.assoc.stats.auc.signif.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.assoc.stats.auc.signif.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.assoc.stats.auc.signif.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.assoc.stats.auc.signif.tissue)[[2]]
    cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    par(las=1, mar=c(6, 4, 6, 0), xaxt="n")
    # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
    yylim <- c(0, 1)
    mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main=sprintf("Gene-drug association (FDR < %i%%)\n%s (%i cell lines)", round(myfdr * 100), gsub("_", " ", tissuen[ii]), ncelline))
    # legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
    # plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
    axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
    text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
    text(x=mp[1, ], y=cors.ic50$rho, pos=ifelse(is.na(cors.ic50$rho) | cors.ic50$rho < 0, 1, 3), labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
    text(x=mp[2, ], y=cors.auc$rho, pos=ifelse(is.na(cors.auc$rho) | cors.auc$rho < 0, 1, 3), labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
  }
}
dev.off()

## multiple barplots per page
count <- 0
nf <- 6
for(ii in 1:length(tissuen)) {
  ncelline <- sum(tissue == tissuen[ii])
  if(max(correlations.assoc.stats.ic50.signif.tissue[ii, , "n"], na.rm=TRUE) >= minsample && max(correlations.assoc.stats.ic50.signif.tissue[ii, , "rho"], na.rm=TRUE) > 0) {
    count <- count + 1
    ## IC50
    xx <- correlations.assoc.stats.ic50.signif.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.assoc.stats.ic50.signif.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.assoc.stats.ic50.signif.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.assoc.stats.ic50.signif.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.assoc.stats.ic50.signif.tissue)[[2]]
    cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    ## AUC
    xx <- correlations.assoc.stats.auc.signif.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.assoc.stats.auc.signif.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.assoc.stats.auc.signif.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.assoc.stats.auc.signif.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.assoc.stats.auc.signif.tissue)[[2]]
    cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    if((count %% nf) == 1) {
      pdf(file.path(saveres, sprintf("barplot_assoc_ic50_auc_signif_tissue%i_cgp_ccle_paper.pdf", ceiling(count / nf))), height=14, width=10)
      par(las=1, mar=c(6, 4, 6, 0), xaxt="n", mfrow=c(3, 2))
    }
    # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
    yylim <- c(0, 1)
    mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main=sprintf("Gene-drug associations (FDR < %i%%)\n%s (%i cell lines)", round(myfdr * 100), gsub("_", " ", tissuen[ii]), ncelline))
    legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
    # plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
    axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
    text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
    text(x=mp[1, ], y=cors.ic50$rho, pos=ifelse(is.na(cors.ic50$rho) | cors.ic50$rho < 0, 1, 3), labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
    text(x=mp[2, ], y=cors.auc$rho, pos=ifelse(is.na(cors.auc$rho) | cors.auc$rho < 0, 1, 3), labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
    if((count %% nf) == 0) { dev.off() }
  }
}
if((count %% nf) != 0) { dev.off() }


###############################################################################
###############################################################################
## compute geneset enrichment score using prerank GSEA
## based on associations between genes and drugs
###############################################################################
###############################################################################

## correlation statistics
tt <- matrix(NA, nrow=length(drugsn), ncol=5, dimnames=list(drugsn, c("rho", "lower", "upper", "p", "n")))
correlations.gsea.stats <- list("ic50"=tt, "ic50.filt"=tt, "ic50.call"=tt, "ic50.call.filt"=tt, "auc"=tt, "auc.filt"=tt, "auc.call"=tt, "auc.call.filt"=tt)

########################
## gene sets

myfn <- file.path(saveres, "genesets.RData")
genesets.filen2 <- unlist(strsplit(x=genesets.filen, split="[.]"))
genesets.filen2 <- paste(c(genesets.filen2[1:(length(genesets.filen2) - 1)], "format", genesets.filen2[length(genesets.filen2)]), collapse=".")
if(!file.exists(myfn)) {
  genesets <- formatGMT(infile=genesets.filen, outfile=file.path(saveres, genesets.filen2), replace=TRUE)
  save(list="genesets", compress=TRUE, file=myfn)
} else { load(myfn) }

########################
## IC50 per tissue types

for(ii in 1:length(tissuen)) {
  ttt <- toupper(gsub(badchars, "", tissuen[ii]))
  myfn <- file.path(saveres, sprintf("gsea_ic50_cgp_ccle_%s.RData", ttt))
  if(!file.exists(myfn)) {
    myfn2 <- file.path(saveres, sprintf("assoc_ic50_cgp_ccle_%s.RData", ttt))
    if(file.exists(myfn2)) {
      ## load gene-drug associations for a specific tissue type
      load(myfn2)
      ## CCLE
      message(sprintf("GSEA based on IC50 in %s (CCLE)", ttt))
      ## create rankings for each drug
      dir.create(file.path(saveres, "GSEA", "rankings", "CCLE", ttt), recursive=TRUE, showWarnings=FALSE)
      rank.files <- NULL
      for(i in 1:length(assoc.ic50.ccle)) {
        ss <- assoc.ic50.ccle[[i]][ , "stat"]
        ss <- sort(ss, decreasing=TRUE, na.last=NA)
        ss[ss == Inf] <- .Machine$double.xmax
        ss[ss == -Inf] <- -.Machine$double.xmax
        rankg <- cbind(annot[names(ss), "jetset.EntrezID"], ss)
        fff <- file.path(saveres, "GSEA", "rankings", "CCLE", ttt, sprintf("ic50_%s.rnk", gsub("drugid_", "", names(assoc.ic50.ccle)[i])))
        write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
        rank.files <- c(rank.files, fff)  
      }
      names(rank.files) <- gsub("drugid_", "", names(assoc.ic50.ccle))
      ## GSEA
      dir.create(file.path(saveres, "GSEA", "reports", "CCLE", ttt), recursive=TRUE, showWarnings=FALSE)
      gsea.res <- mclapply(1:length(rank.files), function(x, y, z, ...) { 
      	tt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=y[x], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=file.path(saveres, "GSEA", "reports", "CCLE", ttt), replace.res=FALSE, gsea.seed=987654321)
        tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)], "description"=z[rownames(tt)])
        tt[!is.na(tt[ , "NOM.p.val"]) & tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (gsea.nperm + 1)
        return(tt)
      }, y=rank.files, z=genesets$description, gsea.exec=gsea.exec, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, saveres=saveres)
      names(gsea.res) <- names(rank.files)
      ## save all class comparisons into a multisheets xls file
      WriteXLS::WriteXLS("gsea.res", ExcelFileName=file.path(saveres, sprintf("ic50_ccle_results_gsea_gene_drug_%s.xls", ttt)))
      gsea.ic50.ccle <- lapply(gsea.res, function(x, y) {
        rr <- matrix(NA, nrow=length(y), ncol=3, dimnames=list(y, c("nes", "pvalue", "fdr")))
        rr[as.character(x[ , "geneset"]), ] <- data.matrix(x[ , c("NES", "NOM.p.val", "FDR.q.val")])
        return(rr)
      }, y=names(genesets$geneset))
      message("")
      ## CGP
      message(sprintf("GSEA based on IC50 in %s (CGP)", ttt))
      ## create rankings for each drug
      dir.create(file.path(saveres, "GSEA", "rankings", "CGP", ttt), recursive=TRUE, showWarnings=FALSE)
      rank.files <- NULL
      for(i in 1:length(assoc.ic50.cgp)) {
        ss <- assoc.ic50.cgp[[i]][ , "stat"]
        ss <- sort(ss, decreasing=TRUE, na.last=NA)
        ss[ss == Inf] <- .Machine$double.xmax
        ss[ss == -Inf] <- -.Machine$double.xmax
        rankg <- cbind(annot[names(ss), "jetset.EntrezID"], ss)
        fff <- file.path(saveres, "GSEA", "rankings", "CGP", ttt, sprintf("ic50_%s.rnk", gsub("drugid_", "", names(assoc.ic50.cgp)[i])))
        write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
        rank.files <- c(rank.files, fff)  
      }
      names(rank.files) <- gsub("drugid_", "", names(assoc.ic50.cgp))
      ## GSEA
      dir.create(file.path(saveres, "GSEA", "reports", "CGP", ttt), recursive=TRUE, showWarnings=FALSE)
      gsea.res <- mclapply(1:length(rank.files), function(x, y, z, ...) { 
      	tt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=y[x], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=file.path(saveres, "GSEA", "reports", "CGP", ttt), replace.res=FALSE, gsea.seed=987654321)
        tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)], "description"=z[rownames(tt)])
         tt[!is.na(tt[ , "NOM.p.val"]) & tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (gsea.nperm + 1)
        return(tt)
      }, y=rank.files, z=genesets$description, gsea.exec=gsea.exec, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, saveres=saveres)
      names(gsea.res) <- names(rank.files)
      ## save all class comparisons into a multisheets xls file
      WriteXLS::WriteXLS("gsea.res", ExcelFileName=file.path(saveres, sprintf("ic50_cgp_results_gsea_gene_drug_%s.xls", ttt)))
      gsea.ic50.cgp <- lapply(gsea.res, function(x, y) {
        rr <- matrix(NA, nrow=length(y), ncol=3, dimnames=list(y, c("nes", "pvalue", "fdr")))
        rr[as.character(x[ , "geneset"]), ] <- data.matrix(x[ , c("NES", "NOM.p.val", "FDR.q.val")])
        return(rr)
      }, y=names(genesets$geneset))
      message("")
      save(list=c("gsea.ic50.ccle", "gsea.ic50.cgp"), compress=TRUE, file=myfn)
      rm(list=c("assoc.ic50.ccle", "assoc.ic50.cgp"))
    }
  }
}

correlations.gsea.stats.ic50.tissue <- array(NA, dim=c(length(tissuen), length(drugn), 5), dimnames=list(tissuen, gsub("drugid_", "", drugn), c("rho", "upper", "lower", "p", "n")))
for(ii in 1:length(tissuen)) {
  myfn <- file.path(saveres, sprintf("gsea_ic50_cgp_ccle_%s.RData", toupper(gsub(badchars, "", tissuen[ii]))))
  if(file.exists(myfn)) {
    load(myfn)
    for(i in 1:length(gsea.ic50.ccle)) {
      nnn <- sum(complete.cases(gsea.ic50.ccle[[i]][ , "nes"], gsea.ic50.cgp[[i]][ , "nes"]))
      if(nnn >= minsample) {
        cc <- cor.test(gsea.ic50.ccle[[i]][ , "nes"], gsea.ic50.cgp[[i]][ , "nes"], method="spearman", use="complete.obs", alternative="greater")
      } else {
        cc <- list("estimate"=NA, "p.value"=NA)
      }
      ## correlation statistics
      cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
      correlations.gsea.stats.ic50.tissue[ii, i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
    }
  }
}

correlations.gsea.stats.ic50.signif.tissue <- array(NA, dim=c(length(tissuen), length(drugn), 5), dimnames=list(tissuen, gsub("drugid_", "", drugn), c("rho", "upper", "lower", "p", "n")))
for(ii in 1:length(tissuen)) {
  myfn <- file.path(saveres, sprintf("gsea_ic50_cgp_ccle_%s.RData", toupper(gsub(badchars, "", tissuen[ii]))))
  if(file.exists(myfn)) {
    load(myfn)
    for(i in 1:length(gsea.ic50.ccle)) {
      myx <- gsea.ic50.ccle[[i]][ , "fdr"] < myfdr | gsea.ic50.cgp[[i]][ , "fdr"] < myfdr
      nnn <- sum(complete.cases(gsea.ic50.ccle[[i]][myx, "nes"], gsea.ic50.cgp[[i]][myx, "nes"]))
      if(nnn >= minsample) {
        cc <- cor.test(gsea.ic50.ccle[[i]][myx, "nes"], gsea.ic50.cgp[[i]][myx, "nes"], method="spearman", use="complete.obs", alternative="greater")
      } else {
        cc <- list("estimate"=NA, "p.value"=NA)
      }
      ## correlation statistics
      cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
      correlations.gsea.stats.ic50.signif.tissue[ii, i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
    }
  }
}

## IC50 for all tissues
## CCLE
message("GSEA based on IC50 (CCLE)")
myfn2 <- file.path(saveres, "gsea_ic50_cgp_ccle.RData")
if(!file.exists(myfn2)) {
  myfn <- file.path(saveres, "assoc_ic50_cgp_ccle.RData")
  load(myfn)
  ## create rankings for each drug
  dir.create(file.path(saveres, "GSEA", "rankings", "CCLE", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
  rank.files <- NULL
  for(i in 1:length(assoc.ic50.ccle)) {
    ss <- assoc.ic50.ccle[[i]][ , "stat"]
    ss <- sort(ss, decreasing=TRUE, na.last=NA)
    ss[ss == Inf] <- .Machine$double.xmax
    ss[ss == -Inf] <- -.Machine$double.xmax
    rankg <- cbind(annot[names(ss), "jetset.EntrezID"], ss)
    fff <- file.path(saveres, "GSEA", "rankings", "CCLE", "ALLTISSUES", sprintf("ic50_%s.rnk", gsub("drugid_", "", names(assoc.ic50.ccle)[i])))
    write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    rank.files <- c(rank.files, fff)  
  }
  names(rank.files) <- gsub("drugid_", "", names(assoc.ic50.ccle))
  ## GSEA
  dir.create(file.path(saveres, "GSEA", "reports", "CCLE", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
  gsea.res <- mclapply(1:length(rank.files), function(x, y, z, ...) { 
  	tt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=y[x], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=file.path(saveres, "GSEA", "reports", "CCLE", "ALLTISSUES"), replace.res=FALSE, gsea.seed=987654321)
    tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)], "description"=z[rownames(tt)])
    tt[tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (gsea.nperm + 1)
    return(tt)
  }, y=rank.files, z=genesets$description, gsea.exec=gsea.exec, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, saveres=saveres)
  names(gsea.res) <- names(rank.files)
  ## save all class comparisons into a multisheets xls file
  WriteXLS::WriteXLS("gsea.res", ExcelFileName=file.path(saveres, "ic50_ccle_results_gsea_gene_drug_paper.xls"))
  gsea.ic50.ccle <- lapply(gsea.res, function(x, y) {
    rr <- matrix(NA, nrow=length(y), ncol=3, dimnames=list(y, c("nes", "pvalue", "fdr")))
    rr[as.character(x[ , "geneset"]), ] <- data.matrix(x[ , c("NES", "NOM.p.val", "FDR.q.val")])
    return(rr)
  }, y=names(genesets$geneset))
  message("")
  ## CGP
  message("GSEA based on IC50 (CGP)")
  ## create rankings for each drug
  dir.create(file.path(saveres, "GSEA", "rankings", "CGP", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
  rank.files <- NULL
  for(i in 1:length(assoc.ic50.cgp)) {
    ss <- assoc.ic50.cgp[[i]][ , "stat"]
    ss <- sort(ss, decreasing=TRUE, na.last=NA)
    ss[ss == Inf] <- .Machine$double.xmax
    ss[ss == -Inf] <- -.Machine$double.xmax
    rankg <- cbind(annot[names(ss), "jetset.EntrezID"], ss)
    fff <- file.path(saveres, "GSEA", "rankings", "CGP", "ALLTISSUES", sprintf("ic50_%s.rnk", gsub("drugid_", "", names(assoc.ic50.cgp)[i])))
    write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    rank.files <- c(rank.files, fff)  
  }
  names(rank.files) <- gsub("drugid_", "", names(assoc.ic50.cgp))
  ## GSEA
  dir.create(file.path(saveres, "GSEA", "reports", "CGP", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
  gsea.res <- mclapply(1:length(rank.files), function(x, y, z, ...) { 
  	tt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=y[x], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=file.path(saveres, "GSEA", "reports", "CGP", "ALLTISSUES"), replace.res=FALSE, gsea.seed=987654321)
    tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)], "description"=z[rownames(tt)])
    tt[tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (gsea.nperm + 1)
    return(tt)
  }, y=rank.files, z=genesets$description, gsea.exec=gsea.exec, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, saveres=saveres)
  names(gsea.res) <- names(rank.files)
  ## save all class comparisons into a multisheets xls file
  WriteXLS::WriteXLS("gsea.res", ExcelFileName=file.path(saveres, "ic50_cgp_results_gsea_gene_drug_paper.xls"))
  gsea.ic50.cgp <- lapply(gsea.res, function(x, y) {
    rr <- matrix(NA, nrow=length(y), ncol=3, dimnames=list(y, c("nes", "pvalue", "fdr")))
    rr[as.character(x[ , "geneset"]), ] <- data.matrix(x[ , c("NES", "NOM.p.val", "FDR.q.val")])
    return(rr)
  }, y=names(genesets$geneset))
  message("")
  save(list=c("gsea.ic50.ccle", "gsea.ic50.cgp"), compress=TRUE, file=myfn2)
} else { load(myfn2) }


## consistency between associations
## all GO terms
message("Correlation enrichment scores based on IC50, all genesets:")
pdf(file.path(saveres, "gsea_ic50_ccle_cgp_all.pdf"), height=10, width=15)
par(mfrow=c(3, 5))
for(i in 1:length(gsea.ic50.ccle)) {
  myx <- complete.cases(gsea.ic50.ccle[[i]][ ,"nes"], gsea.ic50.cgp[[i]][ ,"nes"])
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA)
  } else {
    cc <- cor.test(gsea.ic50.cgp[[i]][ ,"nes"], gsea.ic50.ccle[[i]][ ,"nes"], method="spearman", use="complete.obs", alternative="greater")
  }
  myScatterPlot(x=gsea.ic50.cgp[[i]][ ,"nes"], y=gsea.ic50.ccle[[i]][ ,"nes"], xlim=c(-3, 3), ylim=c(-3, 3), xlab="Enrichment score IC50 (CGP)", ylab="Enrichment score IC50 (CCLE)", main=gsub("drugid_", "", names(gsea.ic50.ccle)[i]), pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.1, yjust=0.8, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, sum(myx)), text.font=2)
  # message("SCC for ", names(gsea.ic50.cgp)[i], ": ", cc$estimate)
  correlations[["ic50"]][i, "go.drug"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=sum(myx), alpha=0.05)
  correlations.gsea.stats[["ic50"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, sum(myx))
}
message("")
dev.off()
## genes with significant gene-drug association in at least one dataset (FDR)
message("Correlation enrichment scores based on IC50, significant genesets:")
pdf(file.path(saveres, "gsea_ic50_ccle_cgp_fdr.pdf"), height=10, width=15)
par(mfrow=c(3, 5))
for(i in 1:length(gsea.ic50.ccle)) {
  myx <- complete.cases(gsea.ic50.ccle[[i]][ ,"fdr"], gsea.ic50.cgp[[i]][ ,"fdr"]) & (gsea.ic50.ccle[[i]][ ,"fdr"] < myfdr | gsea.ic50.cgp[[i]][ ,"fdr"] < myfdr)
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA)
  } else {
    cc <- cor.test(gsea.ic50.cgp[[i]][myx, "nes"], gsea.ic50.ccle[[i]][myx, "nes"], method="spearman", use="complete.obs", alternative="greater")
  }
  myScatterPlot(x=gsea.ic50.cgp[[i]][myx, "nes"], y=gsea.ic50.ccle[[i]][myx, "nes"], xlim=c(-3, 3), ylim=c(-3, 3), xlab="Enrichment score IC50 (CGP)", ylab="Enrichment score IC50 (CCLE)", main=gsub("drugid_", "", names(gsea.ic50.ccle)[i]), pch=16, method="transparent", transparency=0.5)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.1, yjust=0.8, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, sum(myx)), text.font=2)
  # message("SCC for ", names(gsea.ic50.cgp)[i], ": ", cc$estimate)
  correlations[["ic50"]][i, "go.drug.filt"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=sum(myx), alpha=0.05)
  correlations.gsea.stats[["ic50.filt"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, sum(myx))
}
message("")
dev.off()

## boxplot for geneset-drug associations based on ic50
pdf(file.path(saveres, sprintf("boxplot2_gsea_ic50_ccle_cgp_tissue_paper.pdf")), height=8, width=14)
## correlations for each drug per tissue type
ll <- c(list("all_tissues"=correlations.gsea.stats[["ic50"]][ , "rho"]), lapply(apply(correlations.gsea.stats.ic50.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ll <- lapply(ll, function(x) { return(x[!is.na(x)]) })
# ll <- ll[sapply(ll, length) > 0]
myx <- which(sapply(ll, length) < 3)
wt <- kruskal.test(x=ll[sapply(ll, length) > 0])
ll[sapply(ll, length) < 3] <- NA
par(las=1, mar=c(12, 4, 4, 2) + 0.1, xaxt="n")
mp <- boxplot(ll, outline=FALSE, ylab="Rs", main="Correlations of geneset-drug associations (IC50) across tissue types\nCCLE vs. CGP", border="grey50", ylim=c(-1, 1), col="white", range=0)
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=c("lightgrey", coltissuet))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=1:length(ll))
points(x=jitter(unlist(xx), 1), y=unlist(ll), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", coltissuet))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
dev.off()

## test whether correlations are significantly higher in specific tissue types compared to all tissues combined
cor.diff.test <- NULL
ll <- c(list("all_tissues"=correlations.gsea.stats[["ic50"]][ , "rho"]), lapply(apply(correlations.gsea.stats.ic50.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ddnn <- names(ll[[1]])
for(i in 2:length(ll)) {
  if(sum(complete.cases(ll[[1]], y=ll[[i]])) < 3) {
   wt <- list("estimate"=NA, "p.value"=NA) 
  } else {
    wt <- wilcox.test(x=ll[[1]], y=ll[[i]], alternative="less", conf.int=TRUE)
  }
  cor.diff.test <- rbind(cor.diff.test, c(wt$estimate, wt$p.value))
}
dimnames(cor.diff.test) <- list(names(ll)[-1], c("difference", "p"))
# cor.diff.test <- cbind(cor.diff.test, "lfdr"=fdrtool::fdrtool(cor.diff.test[ , 2], statistic="pvalue", verbose=FALSE)$lfdr)
cor.diff.test <- cbind(cor.diff.test, "fdr"=p.adjust(cor.diff.test[ , "p"], method="fdr"))
myx <- complete.cases(cor.diff.test) & (cor.diff.test[ , "p"] < 0.05 | cor.diff.test[ , "fdr"] < myfdr)
tt <- cor.diff.test[myx, , drop=FALSE]
tt <- tt[order(tt[ , "p"], decreasing=FALSE), , drop=FALSE]
message("Are correlations between pathway-drug associations (IC50) in some tissue types\n significantly higher than in all tissues combined?")
if(nrow(tt) == 0) { message("\t-> none") } else { print(tt) }

## boxplot for significant geneset-drug associations based on ic50
pdf(file.path(saveres, sprintf("boxplot2_gsea_ic50_ccle_cgp_signif_tissue_paper.pdf")), height=8, width=14)
## correlations for each drug per tissue type
ll <- c(list("all_tissues"=correlations.gsea.stats[["ic50.filt"]][ , "rho"]), lapply(apply(correlations.gsea.stats.ic50.signif.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ll <- lapply(ll, function(x) { return(x[!is.na(x)]) })
# ll <- ll[sapply(ll, length) > 0]
myx <- which(sapply(ll, length) < 3)
wt <- kruskal.test(x=ll[sapply(ll, length) > 0])
ll[sapply(ll, length) < 3] <- NA
par(las=1, mar=c(12, 4, 4, 2) + 0.1, xaxt="n")
mp <- boxplot(ll, outline=FALSE, ylab="Rs", main=sprintf("Correlations of geneset-drug associations (IC50) across tissue types (FDR < %i%%)\nCCLE vs. CGP", round(myfdr * 100)), border="grey50", ylim=c(-1, 1), col="white", range=0)
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=c("lightgrey", coltissuet))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=1:length(ll))
points(x=jitter(unlist(xx), 1), y=unlist(ll), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", coltissuet))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
dev.off()

## test whether correlations are significantly higher in specific tissue types compared to all tissues combined
cor.diff.test <- NULL
ll <- c(list("all_tissues"=correlations.gsea.stats[["ic50.filt"]][ , "rho"]), lapply(apply(correlations.gsea.stats.ic50.signif.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ddnn <- names(ll[[1]])
for(i in 2:length(ll)) {
  if(sum(complete.cases(ll[[1]], y=ll[[i]])) < 3) {
   wt <- list("estimate"=NA, "p.value"=NA) 
  } else {
    wt <- wilcox.test(x=ll[[1]], y=ll[[i]], alternative="less", conf.int=TRUE)
  }
  cor.diff.test <- rbind(cor.diff.test, c(wt$estimate, wt$p.value))
}
dimnames(cor.diff.test) <- list(names(ll)[-1], c("difference", "p"))
# cor.diff.test <- cbind(cor.diff.test, "lfdr"=fdrtool::fdrtool(cor.diff.test[ , 2], statistic="pvalue", verbose=FALSE)$lfdr)
cor.diff.test <- cbind(cor.diff.test, "fdr"=p.adjust(cor.diff.test[ , "p"], method="fdr"))
myx <- complete.cases(cor.diff.test) & (cor.diff.test[ , "p"] < 0.05 | cor.diff.test[ , "fdr"] < myfdr)
tt <- cor.diff.test[myx, , drop=FALSE]
tt <- tt[order(tt[ , "p"], decreasing=FALSE), , drop=FALSE]
message("Are correlations between significant pathway-drug associations (IC50) in some tissue types\n significantly higher than in all tissues combined?")
if(nrow(tt) == 0) { message("\t-> none") } else { print(tt) }

########################
## IC50 sensitivity calling

## CCLE
message("GSEA based on IC50 sensitivity calling (CCLE)")
## create rankings for each drug
dir.create(file.path(saveres, "GSEA", "rankings", "CCLE", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
rank.files <- NULL
for(i in 1:length(assoc.ic50.call.ccle)) {
  ss <- assoc.ic50.call.ccle[[i]][ , "stat"]
  ss <- sort(ss, decreasing=TRUE, na.last=NA)
  ss[ss == Inf] <- .Machine$double.xmax
  ss[ss == -Inf] <- -.Machine$double.xmax
  rankg <- cbind(annot[names(ss), "jetset.EntrezID"], ss)
  fff <- file.path(saveres, "GSEA", "rankings", "CCLE", "ALLTISSUES", sprintf("ic50_call_%s.rnk", gsub("drugid_", "", names(assoc.ic50.call.ccle)[i])))
  write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
  rank.files <- c(rank.files, fff)  
}
names(rank.files) <- gsub("drugid_", "", names(assoc.ic50.call.ccle))
## GSEA
dir.create(file.path(saveres, "GSEA", "reports", "CCLE", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
gsea.res <- mclapply(1:length(rank.files), function(x, y, z, ...) { 
	tt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=y[x], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=file.path(saveres, "GSEA", "reports", "CCLE", "ALLTISSUES"), replace.res=FALSE, gsea.seed=987654321)
  tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)], "description"=z[rownames(tt)])
  tt[tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (gsea.nperm + 1)
  return(tt)
}, y=rank.files, z=genesets$description, gsea.exec=gsea.exec, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, saveres=saveres)
names(gsea.res) <- names(rank.files)
## save all class comparisons into a multisheets xls file
WriteXLS::WriteXLS("gsea.res", ExcelFileName=file.path(saveres, "ic50_call_ccle_results_gsea_gene_drug.xls"))
gsea.ic50.call.ccle <- lapply(gsea.res, function(x, y) {
  rr <- matrix(NA, nrow=length(y), ncol=3, dimnames=list(y, c("nes", "pvalue", "fdr")))
  rr[as.character(x[ , "geneset"]), ] <- data.matrix(x[ , c("NES", "NOM.p.val", "FDR.q.val")])
  return(rr)
}, y=names(genesets$geneset))
message("")
## CGP
message("GSEA based on IC50 sensitivity calling (CGP)")
## create rankings for each drug
dir.create(file.path(saveres, "GSEA", "rankings", "CGP", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
rank.files <- NULL
for(i in 1:length(assoc.ic50.call.cgp)) {
  ss <- assoc.ic50.call.cgp[[i]][ , "stat"]
  ss <- sort(ss, decreasing=TRUE, na.last=NA)
  ss[ss == Inf] <- .Machine$double.xmax
  ss[ss == -Inf] <- -.Machine$double.xmax
  rankg <- cbind(annot[names(ss), "jetset.EntrezID"], ss)
  fff <- file.path(saveres, "GSEA", "rankings", "CGP", "ALLTISSUES", sprintf("ic50_call_%s.rnk", gsub("drugid_", "", names(assoc.ic50.call.cgp)[i])))
  write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
  rank.files <- c(rank.files, fff)  
}
names(rank.files) <- gsub("drugid_", "", names(assoc.ic50.call.cgp))
## GSEA
dir.create(file.path(saveres, "GSEA", "reports", "CGP", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
gsea.res <- mclapply(1:length(rank.files), function(x, y, z, ...) { 
	tt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=y[x], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=file.path(saveres, "GSEA", "reports", "CGP", "ALLTISSUES"), replace.res=FALSE, gsea.seed=987654321)
  tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)], "description"=z[rownames(tt)])
  tt[tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (gsea.nperm + 1)
  return(tt)
}, y=rank.files, z=genesets$description, gsea.exec=gsea.exec, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, saveres=saveres)
names(gsea.res) <- names(rank.files)
## save all class comparisons into a multisheets xls file
WriteXLS::WriteXLS("gsea.res", ExcelFileName=file.path(saveres, "ic50_call_cgp_results_gsea_gene_drug.xls"))
gsea.ic50.call.cgp <- lapply(gsea.res, function(x, y) {
  rr <- matrix(NA, nrow=length(y), ncol=3, dimnames=list(y, c("nes", "pvalue", "fdr")))
  rr[as.character(x[ , "geneset"]), ] <- data.matrix(x[ , c("NES", "NOM.p.val", "FDR.q.val")])
  return(rr)
}, y=names(genesets$geneset))
message("")

## consistency between enrichment scores
## all GO terms
message("Correlation enrichment scores based on IC50 sensitivity calling, all genesets:")
pdf(file.path(saveres, "gsea_ic50_call_ccle_cgp_all.pdf"), height=10, width=15)
par(mfrow=c(3, 5))
for(i in 1:length(gsea.ic50.call.ccle)) {
  myx <- intersect(rownames(gsea.ic50.call.ccle[[i]]), rownames(gsea.ic50.call.cgp[[i]]))
  myx <- myx[complete.cases(gsea.ic50.call.ccle[[i]][myx,"nes"], gsea.ic50.call.cgp[[i]][myx,"nes"])]
  if(length(myx) <= 3) {
   cc <- list("estimate"=NA, "p.value"=NA)
  } else {
    cc <- cor.test(gsea.ic50.call.cgp[[i]][myx,"nes"], gsea.ic50.call.ccle[[i]][myx,"nes"], method="spearman", use="complete.obs", alternative="greater")
  }  
  myScatterPlot(x=gsea.ic50.call.cgp[[i]][myx,"nes"], y=gsea.ic50.call.ccle[[i]][myx,"nes"], xlim=c(-3, 3), ylim=c(-3, 3), xlab="Enrichment score IC50 calling (CGP)", ylab="Enrichment score IC50 calling (CCLE)", main=gsub("drugid_", "", names(gsea.ic50.call.ccle)[i]), pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.1, yjust=0.8, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, length(myx)), text.font=2)
  # message("SCC for ", names(gsea.ic50.call.cgp)[i], ": ", cc$estimate)
  correlations[["ic50.call"]][i, "go.drug"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=length(myx), alpha=0.05)
  correlations.gsea.stats[["ic50.call"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, length(myx))
}
message("")
dev.off()
## genes with significant gene-drug association in at least one dataset (FDR)
message("Correlation gene-drug associations based on IC50 sensitivity calling, significant genesets:")
pdf(file.path(saveres, "gsea_ic50_call_ccle_cgp_fdr.pdf"), height=10, width=15)
par(mfrow=c(3, 5))
for(i in 1:length(gsea.ic50.call.ccle)) {
  myx <- intersect(rownames(gsea.ic50.call.ccle[[i]]), rownames(gsea.ic50.call.cgp[[i]]))
  myx <- myx[complete.cases(gsea.ic50.call.ccle[[i]][myx, "fdr"], gsea.ic50.call.cgp[[i]][myx, "fdr"]) & (gsea.ic50.call.ccle[[i]][myx, "fdr"] < myfdr | gsea.ic50.call.cgp[[i]][myx, "fdr"] < myfdr)]
  if(length(myx) <= 3) {
    cc <- list("estimate"=NA, "p.value"=NA)
  } else {
    cc <- cor.test(gsea.ic50.call.cgp[[i]][myx, "nes"], gsea.ic50.call.ccle[[i]][myx, "nes"], method="spearman", use="complete.obs", alternative="greater")
  }
  myScatterPlot(x=gsea.ic50.call.cgp[[i]][myx, "nes"], y=gsea.ic50.call.ccle[[i]][myx, "nes"], xlim=c(-3, 3), ylim=c(-3, 3), xlab="Enrichment score IC50 calling (CGP)", ylab="Enrichment score IC50 calling (CCLE)", main=gsub("drugid_", "", names(gsea.ic50.call.ccle)[i]), pch=16, method="transparent", transparency=0.5)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.1, yjust=0.8, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, length(myx)), text.font=2)
  # message("SCC for ", names(gsea.ic50.call.cgp)[i], ": ", cc$estimate)
  correlations[["ic50.call"]][i, "go.drug.filt"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=length(myx), alpha=0.05)
  correlations.gsea.stats[["ic50.call.filt"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, length(myx))
}
message("")
dev.off()

########################
## AUC per tissue types
for(ii in 1:length(tissuen)) {
  ttt <- toupper(gsub(badchars, "", tissuen[ii]))
  myfn <- file.path(saveres, sprintf("gsea_auc_cgp_ccle_%s.RData", ttt))
  if(!file.exists(myfn)) {
    myfn2 <- file.path(saveres, sprintf("assoc_auc_cgp_ccle_%s.RData", ttt))
    if(file.exists(myfn2)) {
      ## load gene-drug associations for a specific tissue type
      load(myfn2)
      ## CCLE
      message(sprintf("GSEA based on AUC in %s (CCLE)", ttt))
      ## create rankings for each drug
      dir.create(file.path(saveres, "GSEA", "rankings", "CCLE", ttt), recursive=TRUE, showWarnings=FALSE)
      rank.files <- NULL
      for(i in 1:length(assoc.auc.ccle)) {
        ss <- assoc.auc.ccle[[i]][ , "stat"]
        ss <- sort(ss, decreasing=TRUE, na.last=NA)
        ss[ss == Inf] <- .Machine$double.xmax
        ss[ss == -Inf] <- -.Machine$double.xmax
        rankg <- cbind(annot[names(ss), "jetset.EntrezID"], ss)
        fff <- file.path(saveres, "GSEA", "rankings", "CCLE", ttt, sprintf("auc_%s.rnk", gsub("drugid_", "", names(assoc.auc.ccle)[i])))
        write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
        rank.files <- c(rank.files, fff)  
      }
      names(rank.files) <- gsub("drugid_", "", names(assoc.auc.ccle))
      ## GSEA
      dir.create(file.path(saveres, "GSEA", "reports", "CCLE", ttt), recursive=TRUE, showWarnings=FALSE)
      gsea.res <- mclapply(1:length(rank.files), function(x, y, z, ...) { 
      	tt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=y[x], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=file.path(saveres, "GSEA", "reports", "CCLE", ttt), replace.res=FALSE, gsea.seed=987654321)
        tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)], "description"=z[rownames(tt)])
        tt[!is.na(tt[ , "NOM.p.val"]) & tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (gsea.nperm + 1)
        return(tt)
      }, y=rank.files, z=genesets$description, gsea.exec=gsea.exec, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, saveres=saveres)
      names(gsea.res) <- names(rank.files)
      ## save all class comparisons into a multisheets xls file
      WriteXLS::WriteXLS("gsea.res", ExcelFileName=file.path(saveres, sprintf("auc_ccle_results_gsea_gene_drug_%s.xls", ttt)))
      gsea.auc.ccle <- lapply(gsea.res, function(x, y) {
        rr <- matrix(NA, nrow=length(y), ncol=3, dimnames=list(y, c("nes", "pvalue", "fdr")))
        rr[as.character(x[ , "geneset"]), ] <- data.matrix(x[ , c("NES", "NOM.p.val", "FDR.q.val")])
        return(rr)
      }, y=names(genesets$geneset))
      message("")
      ## CGP
      message(sprintf("GSEA based on AUC in %s (CGP)", ttt))
      ## create rankings for each drug
      dir.create(file.path(saveres, "GSEA", "rankings", "CGP", ttt), recursive=TRUE, showWarnings=FALSE)
      rank.files <- NULL
      for(i in 1:length(assoc.auc.cgp)) {
        ss <- assoc.auc.cgp[[i]][ , "stat"]
        ss <- sort(ss, decreasing=TRUE, na.last=NA)
        ss[ss == Inf] <- .Machine$double.xmax
        ss[ss == -Inf] <- -.Machine$double.xmax
        rankg <- cbind(annot[names(ss), "jetset.EntrezID"], ss)
        fff <- file.path(saveres, "GSEA", "rankings", "CGP", ttt, sprintf("auc_%s.rnk", gsub("drugid_", "", names(assoc.auc.cgp)[i])))
        write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
        rank.files <- c(rank.files, fff)  
      }
      names(rank.files) <- gsub("drugid_", "", names(assoc.auc.cgp))
      ## GSEA
      dir.create(file.path(saveres, "GSEA", "reports", "CGP", ttt), recursive=TRUE, showWarnings=FALSE)
      gsea.res <- mclapply(1:length(rank.files), function(x, y, z, ...) { 
      	tt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=y[x], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=file.path(saveres, "GSEA", "reports", "CGP", ttt), replace.res=FALSE, gsea.seed=987654321)
        tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)], "description"=z[rownames(tt)])
         tt[!is.na(tt[ , "NOM.p.val"]) & tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (gsea.nperm + 1)
        return(tt)
      }, y=rank.files, z=genesets$description, gsea.exec=gsea.exec, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, saveres=saveres)
      names(gsea.res) <- names(rank.files)
      ## save all class comparisons into a multisheets xls file
      WriteXLS::WriteXLS("gsea.res", ExcelFileName=file.path(saveres, sprintf("auc_cgp_results_gsea_gene_drug_%s.xls", ttt)))
      gsea.auc.cgp <- lapply(gsea.res, function(x, y) {
        rr <- matrix(NA, nrow=length(y), ncol=3, dimnames=list(y, c("nes", "pvalue", "fdr")))
        rr[as.character(x[ , "geneset"]), ] <- data.matrix(x[ , c("NES", "NOM.p.val", "FDR.q.val")])
        return(rr)
      }, y=names(genesets$geneset))
      message("")
      save(list=c("gsea.auc.ccle", "gsea.auc.cgp"), compress=TRUE, file=myfn)
      rm(list=c("assoc.auc.ccle", "assoc.auc.cgp"))
    }
  }
}

correlations.gsea.stats.auc.tissue <- array(NA, dim=c(length(tissuen), length(drugn), 5), dimnames=list(tissuen, gsub("drugid_", "", drugn), c("rho", "upper", "lower", "p", "n")))
for(ii in 1:length(tissuen)) {
  myfn <- file.path(saveres, sprintf("gsea_auc_cgp_ccle_%s.RData", toupper(gsub(badchars, "", tissuen[ii]))))
  if(file.exists(myfn)) {
    load(myfn)
    for(i in 1:length(gsea.auc.ccle)) {
      nnn <- sum(complete.cases(gsea.auc.ccle[[i]][ , "nes"], gsea.auc.cgp[[i]][ , "nes"]))
      if(nnn >= minsample) {
        cc <- cor.test(gsea.auc.ccle[[i]][ , "nes"], gsea.auc.cgp[[i]][ , "nes"], method="spearman", use="complete.obs", alternative="greater")
      } else {
        cc <- list("estimate"=NA, "p.value"=NA)
      }
      ## correlation statistics
      cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
      correlations.gsea.stats.auc.tissue[ii, i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
    }
  }
}

correlations.gsea.stats.auc.signif.tissue <- array(NA, dim=c(length(tissuen), length(drugn), 5), dimnames=list(tissuen, gsub("drugid_", "", drugn), c("rho", "upper", "lower", "p", "n")))
for(ii in 1:length(tissuen)) {
  myfn <- file.path(saveres, sprintf("gsea_auc_cgp_ccle_%s.RData", toupper(gsub(badchars, "", tissuen[ii]))))
  if(file.exists(myfn)) {
    load(myfn)
    for(i in 1:length(gsea.auc.ccle)) {
      myx <- gsea.auc.ccle[[i]][ , "fdr"] < myfdr | gsea.auc.cgp[[i]][ , "fdr"] < myfdr
      nnn <- sum(complete.cases(gsea.auc.ccle[[i]][myx, "nes"], gsea.auc.cgp[[i]][myx, "nes"]))
      if(nnn >= minsample) {
        cc <- cor.test(gsea.auc.ccle[[i]][myx, "nes"], gsea.auc.cgp[[i]][myx, "nes"], method="spearman", use="complete.obs", alternative="greater")
      } else {
        cc <- list("estimate"=NA, "p.value"=NA)
      }
      ## correlation statistics
      cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
      correlations.gsea.stats.auc.signif.tissue[ii, i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
    }
  }
}

## AUC for all tissues
## CCLE
message("GSEA based on AUC (CCLE)")
myfn2 <- file.path(saveres, "gsea_auc_cgp_ccle.RData")
if(!file.exists(myfn2)) {
  myfn <- file.path(saveres, "assoc_auc_cgp_ccle.RData")
  load(myfn)
  ## create rankings for each drug
  dir.create(file.path(saveres, "GSEA", "rankings", "CCLE", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
  rank.files <- NULL
  for(i in 1:length(assoc.auc.ccle)) {
    ss <- assoc.auc.ccle[[i]][ , "stat"]
    ss <- sort(ss, decreasing=TRUE, na.last=NA)
    ss[ss == Inf] <- .Machine$double.xmax
    ss[ss == -Inf] <- -.Machine$double.xmax
    rankg <- cbind(annot[names(ss), "jetset.EntrezID"], ss)
    fff <- file.path(saveres, "GSEA", "rankings", "CCLE", "ALLTISSUES", sprintf("auc_%s.rnk", gsub("drugid_", "", names(assoc.auc.ccle)[i])))
    write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    rank.files <- c(rank.files, fff)  
  }
  names(rank.files) <- gsub("drugid_", "", names(assoc.auc.ccle))
  ## GSEA
  dir.create(file.path(saveres, "GSEA", "reports", "CCLE", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
  gsea.res <- mclapply(1:length(rank.files), function(x, y, z, ...) { 
  	tt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=y[x], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=file.path(saveres, "GSEA", "reports", "CCLE", "ALLTISSUES"), replace.res=FALSE, gsea.seed=987654321)
    tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)], "description"=z[rownames(tt)])
    tt[tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (gsea.nperm + 1)
    return(tt)
  }, y=rank.files, z=genesets$description, gsea.exec=gsea.exec, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, saveres=saveres)
  names(gsea.res) <- names(rank.files)
  ## save all class comparisons into a multisheets xls file
  WriteXLS::WriteXLS("gsea.res", ExcelFileName=file.path(saveres, "auc_ccle_results_gsea_gene_drug_paper.xls"))
  gsea.auc.ccle <- lapply(gsea.res, function(x, y) {
    rr <- matrix(NA, nrow=length(y), ncol=3, dimnames=list(y, c("nes", "pvalue", "fdr")))
    rr[as.character(x[ , "geneset"]), ] <- data.matrix(x[ , c("NES", "NOM.p.val", "FDR.q.val")])
    return(rr)
  }, y=names(genesets$geneset))
  message("")
  ## CGP
  message("GSEA based on AUC (CGP)")
  ## create rankings for each drug
  dir.create(file.path(saveres, "GSEA", "rankings", "CGP", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
  rank.files <- NULL
  for(i in 1:length(assoc.auc.cgp)) {
    ss <- assoc.auc.cgp[[i]][ , "stat"]
    ss <- sort(ss, decreasing=TRUE, na.last=NA)
    ss[ss == Inf] <- .Machine$double.xmax
    ss[ss == -Inf] <- -.Machine$double.xmax
    rankg <- cbind(annot[names(ss), "jetset.EntrezID"], ss)
    fff <- file.path(saveres, "GSEA", "rankings", "CGP", "ALLTISSUES", sprintf("auc_%s.rnk", gsub("drugid_", "", names(assoc.auc.cgp)[i])))
    write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    rank.files <- c(rank.files, fff)  
  }
  names(rank.files) <- gsub("drugid_", "", names(assoc.auc.cgp))
  ## GSEA
  dir.create(file.path(saveres, "GSEA", "reports", "CGP", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
  gsea.res <- mclapply(1:length(rank.files), function(x, y, z, ...) { 
  	tt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=y[x], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=file.path(saveres, "GSEA", "reports", "CGP", "ALLTISSUES"), replace.res=FALSE, gsea.seed=987654321)
    tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)], "description"=z[rownames(tt)])
    tt[tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (gsea.nperm + 1)
    return(tt)
  }, y=rank.files, z=genesets$description, gsea.exec=gsea.exec, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, saveres=saveres)
  names(gsea.res) <- names(rank.files)
  ## save all class comparisons into a multisheets xls file
  WriteXLS::WriteXLS("gsea.res", ExcelFileName=file.path(saveres, "auc_cgp_results_gsea_gene_drug_paper.xls"))
  gsea.auc.cgp <- lapply(gsea.res, function(x, y) {
    rr <- matrix(NA, nrow=length(y), ncol=3, dimnames=list(y, c("nes", "pvalue", "fdr")))
    rr[as.character(x[ , "geneset"]), ] <- data.matrix(x[ , c("NES", "NOM.p.val", "FDR.q.val")])
    return(rr)
  }, y=names(genesets$geneset))
  message("")
  save(list=c("gsea.auc.ccle", "gsea.auc.cgp"), compress=TRUE, file=myfn2)
} else { load(myfn2) }

## consistency between associations
## all GO terms
message("Correlation enrichment scores based on AUC, all genesets:")
pdf(file.path(saveres, "gsea_auc_ccle_cgp_all.pdf"), height=10, width=15)
par(mfrow=c(3, 5))
for(i in 1:length(gsea.auc.ccle)) {
  myx <- complete.cases(gsea.auc.ccle[[i]][ ,"nes"], gsea.auc.cgp[[i]][ ,"nes"])
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA)
  } else {
    cc <- cor.test(gsea.auc.cgp[[i]][ ,"nes"], gsea.auc.ccle[[i]][ ,"nes"], method="spearman", use="complete.obs", alternative="greater")
  }  
  myScatterPlot(x=gsea.auc.cgp[[i]][ ,"nes"], y=gsea.auc.ccle[[i]][ ,"nes"], xlim=c(-3, 3), ylim=c(-3, 3), xlab="Enrichment score AUC (CGP)", ylab="Enrichment score AUC (CCLE)", main=gsub("drugid_", "", names(gsea.auc.ccle)[i]), pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.1, yjust=0.8, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, sum(myx)), text.font=2)
  # message("SCC for ", names(gsea.auc.cgp)[i], ": ", cc$estimate)
  correlations[["auc"]][i, "go.drug"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=sum(myx), alpha=0.05)
  correlations.gsea.stats[["auc"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, sum(myx))
}
message("")
dev.off()
## genes with significant gene-drug association in at least one dataset (FDR)
message("Correlation gene-drug associations based on AUC, significant genesets:")
pdf(file.path(saveres, "gsea_auc_ccle_cgp_fdr.pdf"), height=10, width=15)
par(mfrow=c(3, 5))
for(i in 1:length(gsea.auc.ccle)) {
  myx <- complete.cases(gsea.auc.ccle[[i]][ ,"fdr"], gsea.auc.cgp[[i]][ ,"fdr"]) & (gsea.auc.ccle[[i]][ ,"fdr"] < myfdr | gsea.auc.cgp[[i]][ ,"fdr"] < myfdr)
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA)
  } else {
    cc <- cor.test(gsea.auc.cgp[[i]][myx, "nes"], gsea.auc.ccle[[i]][myx, "nes"], method="spearman", use="complete.obs", alternative="greater")
  }
  myScatterPlot(x=gsea.auc.cgp[[i]][myx, "nes"], y=gsea.auc.ccle[[i]][myx, "nes"], xlim=c(-3, 3), ylim=c(-3, 3), xlab="Enrichment score AUC (CGP)", ylab="Enrichment score AUC (CCLE)", main=gsub("drugid_", "", names(gsea.auc.ccle)[i]), pch=16, method="transparent", transparency=0.5)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.1, yjust=0.8, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, sum(myx)), text.font=2)
  # message("SCC for ", names(gsea.auc.cgp)[i], ": ", cc$estimate)
  correlations[["auc"]][i, "go.drug.filt"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=sum(myx), alpha=0.05)
  correlations.gsea.stats[["auc.filt"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, sum(myx))
}
message("")
dev.off()


## boxplot for geneset-drug associations based on auc
pdf(file.path(saveres, sprintf("boxplot2_gsea_auc_ccle_cgp_tissue_paper.pdf")), height=8, width=14)
## correlations for each drug per tissue type
ll <- c(list("all_tissues"=correlations.gsea.stats[["auc"]][ , "rho"]), lapply(apply(correlations.gsea.stats.auc.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ll <- lapply(ll, function(x) { return(x[!is.na(x)]) })
# ll <- ll[sapply(ll, length) > 0]
myx <- which(sapply(ll, length) < 3)
wt <- kruskal.test(x=ll[sapply(ll, length) > 0])
ll[sapply(ll, length) < 3] <- NA
par(las=1, mar=c(12, 4, 4, 2) + 0.1, xaxt="n")
mp <- boxplot(ll, outline=FALSE, ylab="Rs", main="Correlations of geneset-drug associations (AUC) across tissue types\nCCLE vs. CGP", border="grey50", ylim=c(-1, 1), col="white", range=0)
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=c("lightgrey", coltissuet))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=1:length(ll))
points(x=jitter(unlist(xx), 1), y=unlist(ll), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", coltissuet))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
dev.off()

## test whether correlations are significantly higher in specific tissue types compared to all tissues combined
cor.diff.test <- NULL
ll <- c(list("all_tissues"=correlations.gsea.stats[["auc"]][ , "rho"]), lapply(apply(correlations.gsea.stats.auc.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ddnn <- names(ll[[1]])
for(i in 2:length(ll)) {
  if(sum(complete.cases(ll[[1]], y=ll[[i]])) < 3) {
   wt <- list("estimate"=NA, "p.value"=NA) 
  } else {
    wt <- wilcox.test(x=ll[[1]], y=ll[[i]], alternative="less", conf.int=TRUE)
  }
  cor.diff.test <- rbind(cor.diff.test, c(wt$estimate, wt$p.value))
}
dimnames(cor.diff.test) <- list(names(ll)[-1], c("difference", "p"))
# cor.diff.test <- cbind(cor.diff.test, "lfdr"=fdrtool::fdrtool(cor.diff.test[ , 2], statistic="pvalue", verbose=FALSE)$lfdr)
cor.diff.test <- cbind(cor.diff.test, "fdr"=p.adjust(cor.diff.test[ , "p"], method="fdr"))
myx <- complete.cases(cor.diff.test) & (cor.diff.test[ , "p"] < 0.05 | cor.diff.test[ , "fdr"] < myfdr)
tt <- cor.diff.test[myx, , drop=FALSE]
tt <- tt[order(tt[ , "p"], decreasing=FALSE), , drop=FALSE]
message("Are correlations between pathway-drug associations (AUC) in some tissue types\n significantly higher than in all tissues combined?")
if(nrow(tt) == 0) { message("\t-> none") } else { print(tt) }

## boxplot for significant geneset-drug associations based on auc
pdf(file.path(saveres, sprintf("boxplot2_gsea_auc_ccle_cgp_signif_tissue_paper.pdf")), height=8, width=14)
## correlations for each drug per tissue type
ll <- c(list("all_tissues"=correlations.gsea.stats[["auc.filt"]][ , "rho"]), lapply(apply(correlations.gsea.stats.auc.signif.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ll <- lapply(ll, function(x) { return(x[!is.na(x)]) })
# ll <- ll[sapply(ll, length) > 0]
myx <- which(sapply(ll, length) < 3)
wt <- kruskal.test(x=ll[sapply(ll, length) > 0])
ll[sapply(ll, length) < 3] <- NA
par(las=1, mar=c(12, 4, 4, 2) + 0.1, xaxt="n")
mp <- boxplot(ll, outline=FALSE, ylab="Rs", main=sprintf("Correlations of geneset-drug associations (AUC) across tissue types (FDR < %i%%)\nCCLE vs. CGP", round(myfdr * 100)), border="grey50", ylim=c(-1, 1), col="white", range=0)
ccol <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=c("lightgrey", coltissuet))
xx <- mapply(FUN=function(x, y) {
    if(x > 0) {
      res <- rep(y, each=x)
    } else { res <- NA }
    return(res)
  }, x=sapply(ll, length), y=1:length(ll))
points(x=jitter(unlist(xx), 1), y=unlist(ll), cex=1.25, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black", coltissuet))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
dev.off()

## test whether correlations are significantly higher in specific tissue types compared to all tissues combined
cor.diff.test <- NULL
ll <- c(list("all_tissues"=correlations.gsea.stats[["auc.filt"]][ , "rho"]), lapply(apply(correlations.gsea.stats.auc.signif.tissue[ , , "rho"], 1, list), function(x) { return(x[[1]]) }))
ddnn <- names(ll[[1]])
for(i in 2:length(ll)) {
  if(sum(complete.cases(ll[[1]], y=ll[[i]])) < 3) {
   wt <- list("estimate"=NA, "p.value"=NA) 
  } else {
    wt <- wilcox.test(x=ll[[1]], y=ll[[i]], alternative="less", conf.int=TRUE)
  }
  cor.diff.test <- rbind(cor.diff.test, c(wt$estimate, wt$p.value))
}
dimnames(cor.diff.test) <- list(names(ll)[-1], c("difference", "p"))
# cor.diff.test <- cbind(cor.diff.test, "lfdr"=fdrtool::fdrtool(cor.diff.test[ , 2], statistic="pvalue", verbose=FALSE)$lfdr)
cor.diff.test <- cbind(cor.diff.test, "fdr"=p.adjust(cor.diff.test[ , "p"], method="fdr"))
myx <- complete.cases(cor.diff.test) & (cor.diff.test[ , "p"] < 0.05 | cor.diff.test[ , "fdr"] < myfdr)
tt <- cor.diff.test[myx, , drop=FALSE]
tt <- tt[order(tt[ , "p"], decreasing=FALSE), , drop=FALSE]
message("Are correlations between significant pathway-drug associations (AUC) in some tissue types\n significantly higher than in all tissues combined?")
if(nrow(tt) == 0) { message("\t-> none") } else { print(tt) }

########################
## AUC sensitivity calling

## CCLE
message("GSEA based on AUC sensitivity calling (CCLE)")
## create rankings for each drug
dir.create(file.path(saveres, "GSEA", "rankings", "CCLE", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
rank.files <- NULL
for(i in 1:length(assoc.auc.call.ccle)) {
  ss <- assoc.auc.call.ccle[[i]][ , "stat"]
  ss <- sort(ss, decreasing=TRUE, na.last=NA)
  ss[ss == Inf] <- .Machine$double.xmax
  ss[ss == -Inf] <- -.Machine$double.xmax
  rankg <- cbind(annot[names(ss), "jetset.EntrezID"], ss)
  fff <- file.path(saveres, "GSEA", "rankings", "CCLE", "ALLTISSUES", sprintf("auc_call_%s.rnk", gsub("drugid_", "", names(assoc.auc.call.ccle)[i])))
  write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
  rank.files <- c(rank.files, fff)  
}
names(rank.files) <- gsub("drugid_", "", names(assoc.auc.call.ccle))
## GSEA
dir.create(file.path(saveres, "GSEA", "reports", "CCLE", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
gsea.res <- mclapply(1:length(rank.files), function(x, y, z, ...) { 
	tt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=y[x], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=file.path(saveres, "GSEA", "reports", "CCLE", "ALLTISSUES"), replace.res=FALSE, gsea.seed=987654321)
  tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)], "description"=z[rownames(tt)])
  tt[tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (gsea.nperm + 1)
  return(tt)
}, y=rank.files, z=genesets$description, gsea.exec=gsea.exec, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, saveres=saveres)
names(gsea.res) <- names(rank.files)
## save all class comparisons into a multisheets xls file
WriteXLS::WriteXLS("gsea.res", ExcelFileName=file.path(saveres, "auc_call_ccle_results_gsea_gene_drug.xls"))
gsea.auc.call.ccle <- lapply(gsea.res, function(x, y) {
  rr <- matrix(NA, nrow=length(y), ncol=3, dimnames=list(y, c("nes", "pvalue", "fdr")))
  rr[as.character(x[ , "geneset"]), ] <- data.matrix(x[ , c("NES", "NOM.p.val", "FDR.q.val")])
  return(rr)
}, y=names(genesets$geneset))
message("")
## CGP
message("GSEA based on AUC sensitivity calling (CGP)")
## create rankings for each drug
dir.create(file.path(saveres, "GSEA", "rankings", "CGP", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
rank.files <- NULL
for(i in 1:length(assoc.auc.call.cgp)) {
  ss <- assoc.auc.call.cgp[[i]][ , "stat"]
  ss <- sort(ss, decreasing=TRUE, na.last=NA)
  ss[ss == Inf] <- .Machine$double.xmax
  ss[ss == -Inf] <- -.Machine$double.xmax
  rankg <- cbind(annot[names(ss), "jetset.EntrezID"], ss)
  fff <- file.path(saveres, "GSEA", "rankings", "CGP", "ALLTISSUES", sprintf("auc_call_%s.rnk", gsub("drugid_", "", names(assoc.auc.call.cgp)[i])))
  write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
  rank.files <- c(rank.files, fff)  
}
names(rank.files) <- gsub("drugid_", "", names(assoc.auc.call.cgp))
## GSEA
dir.create(file.path(saveres, "GSEA", "reports", "CGP", "ALLTISSUES"), recursive=TRUE, showWarnings=FALSE)
gsea.res <- mclapply(1:length(rank.files), function(x, y, z, ...) { 
	tt <- gsea.prerank(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=y[x], gsea.collapse=FALSE, nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.out=file.path(saveres, "GSEA", "reports", "CGP", "ALLTISSUES"), replace.res=FALSE, gsea.seed=987654321)
  tt <- data.frame("geneset"=rownames(tt), tt[ , -c(1, 2, 3, 12)], "description"=z[rownames(tt)])
  tt[tt[ , "NOM.p.val"] == 0, "NOM.p.val"] <- 1 / (gsea.nperm + 1)
  return(tt)
}, y=rank.files, z=genesets$description, gsea.exec=gsea.exec, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size, saveres=saveres)
names(gsea.res) <- names(rank.files)
## save all class comparisons into a multisheets xls file
WriteXLS::WriteXLS("gsea.res", ExcelFileName=file.path(saveres, "auc_call_cgp_results_gsea_gene_drug.xls"))
gsea.auc.call.cgp <- lapply(gsea.res, function(x, y) {
  rr <- matrix(NA, nrow=length(y), ncol=3, dimnames=list(y, c("nes", "pvalue", "fdr")))
  rr[as.character(x[ , "geneset"]), ] <- data.matrix(x[ , c("NES", "NOM.p.val", "FDR.q.val")])
  return(rr)
}, y=names(genesets$geneset))
message("")

## consistency between enrichment scores
## all GO terms
message("Correlation enrichment scores based on AUC sensitivity calling, all genes:")
pdf(file.path(saveres, "gsea_auc_call_ccle_cgp_all.pdf"), height=10, width=15)
par(mfrow=c(3, 5))
for(i in 1:length(gsea.auc.call.ccle)) {
  myx <- intersect(rownames(gsea.auc.call.ccle[[i]]), rownames(gsea.auc.call.cgp[[i]]))
  myx <- myx[complete.cases(gsea.auc.call.ccle[[i]][myx,"nes"], gsea.auc.call.cgp[[i]][myx,"nes"])]
  if(length(myx) <= 3) {
   cc <- list("estimate"=NA, "p.value"=NA)
  } else {
    cc <- cor.test(gsea.auc.call.cgp[[i]][myx,"nes"], gsea.auc.call.ccle[[i]][myx,"nes"], method="spearman", use="complete.obs", alternative="greater")
  }  
  myScatterPlot(x=gsea.auc.call.cgp[[i]][myx,"nes"], y=gsea.auc.call.ccle[[i]][myx,"nes"], xlim=c(-3, 3), ylim=c(-3, 3), xlab="Enrichment score AUC calling (CGP)", ylab="Enrichment score AUC calling (CCLE)", main=gsub("drugid_", "", names(gsea.auc.call.ccle)[i]), pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.1, yjust=0.8, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, length(myx)), text.font=2)
  # message("SCC for ", names(gsea.auc.call.cgp)[i], ": ", cc$estimate)
  correlations[["auc.call"]][i, "go.drug"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=length(myx), alpha=0.05)
  correlations.gsea.stats[["auc.call"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, length(myx))
}
message("")
dev.off()
## genes with significant gene-drug association in at least one dataset (FDR)
message("Correlation enrichment acores based on AUC sensitivity calling, significant genes:")
pdf(file.path(saveres, "gsea_auc_call_ccle_cgp_fdr.pdf"), height=10, width=15)
par(mfrow=c(3, 5))
for(i in 1:length(gsea.auc.call.ccle)) {
  myx <- intersect(rownames(gsea.auc.call.ccle[[i]]), rownames(gsea.auc.call.cgp[[i]]))
  myx <- myx[complete.cases(gsea.auc.call.ccle[[i]][myx, "fdr"], gsea.auc.call.cgp[[i]][myx, "fdr"]) & (gsea.auc.call.ccle[[i]][myx, "fdr"] < myfdr | gsea.auc.call.cgp[[i]][myx, "fdr"] < myfdr)]
  if(length(myx) <= 3) {
   cc <- list("estimate"=NA, "p.value"=NA)
  } else {
    cc <- cor.test(gsea.auc.call.cgp[[i]][myx, "nes"], gsea.auc.call.ccle[[i]][myx, "nes"], method="spearman", use="complete.obs", alternative="greater")
  }
  myScatterPlot(x=gsea.auc.call.cgp[[i]][myx, "nes"], y=gsea.auc.call.ccle[[i]][myx, "nes"], xlim=c(-3, 3), ylim=c(-3, 3), xlab="Enrichment score AUC calling (CGP)", ylab="Enrichment score AUC calling (CCLE)", main=gsub("drugid_", "", names(gsea.auc.call.ccle)[i]), pch=16, method="transparent", transparency=0.5)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.1, yjust=0.8, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, length(myx)), text.font=2)
  # message("SCC for ", names(gsea.auc.call.cgp)[i], ": ", cc$estimate)
  correlations[["auc.call"]][i, "go.drug.filt"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=length(myx), alpha=0.05)
  correlations.gsea.stats[["auc.call.filt"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, length(myx))
}
message("")
dev.off()

## save all correlations statistics
savecor <- NULL
for(i in 1:length(correlations.gsea.stats)) {
  tt <- correlations.gsea.stats[[i]][ , apply(correlations.gsea.stats[[i]], 2, function(x) { return(!all(is.na(x))) }), drop=FALSE]
  savecor <- c(savecor, list(data.frame(tt)))
}
names(savecor) <- names(correlations.gsea.stats)
WriteXLS::WriteXLS("savecor", ExcelFileName=file.path(saveres, "correlations_gsea_stats.xls"), row.names=TRUE)


save(list=c("correlations", "correlations.gsea.stats"), compress=TRUE, file=file.path(saveres, "correlations_gsea_stats.RData"))
save(list=c("correlations.gsea.stats.ic50.tissue", "correlations.gsea.stats.auc.tissue", "correlations.gsea.stats.ic50.signif.tissue", "correlations.gsea.stats.auc.signif.tissue"), compress=TRUE, file=file.path(saveres, "correlations_gsea_stats_tissue.RData"))

## barplot for geneset-drug associations across tissue types
pdf(file.path(saveres, "barplot_gsea_ic50_auc_tissue_cgp_ccle.pdf"), height=6, width=6)
for(ii in 1:length(tissuen)) {
  ncelline <- sum(tissue == tissuen[ii])
  if(max(correlations.gsea.stats.ic50.tissue[ii, , "n"], na.rm=TRUE) >= minsample && max(correlations.gsea.stats.ic50.tissue[ii, , "rho"], na.rm=TRUE) > 0) {
    ## IC50
    xx <- correlations.gsea.stats.ic50.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.gsea.stats.ic50.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.gsea.stats.ic50.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.gsea.stats.ic50.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.gsea.stats.ic50.tissue)[[2]]
    cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    ## AUC
    xx <- correlations.gsea.stats.auc.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.gsea.stats.auc.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.gsea.stats.auc.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.gsea.stats.auc.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.gsea.stats.auc.tissue)[[2]]
    cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    par(las=1, mar=c(6, 4, 6, 0), xaxt="n")
    # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
    yylim <- c(0, 1)
    mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main=sprintf("Pathway-drug associations (all)\n%s (%i cell lines)", gsub("_", " ", tissuen[ii]), ncelline))
    # legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
    # plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
    axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
    text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
    text(x=mp[1, ], y=cors.ic50$rho, pos=ifelse(is.na(cors.ic50$rho) | cors.ic50$rho < 0, 1, 3), labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
    text(x=mp[2, ], y=cors.auc$rho, pos=ifelse(is.na(cors.auc$rho) | cors.auc$rho < 0, 1, 3), labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
  }
}
dev.off()

## multiple barplots per page
count <- 0
nf <- 6
for(ii in 1:length(tissuen)) {
  ncelline <- sum(tissue == tissuen[ii])
  if(max(correlations.gsea.stats.ic50.tissue[ii, , "n"], na.rm=TRUE) >= minsample && max(correlations.gsea.stats.ic50.tissue[ii, , "rho"], na.rm=TRUE) > 0) {
    count <- count + 1
    ## IC50
    xx <- correlations.gsea.stats.ic50.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.gsea.stats.ic50.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.gsea.stats.ic50.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.gsea.stats.ic50.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.gsea.stats.ic50.tissue)[[2]]
    cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    ## AUC
    xx <- correlations.gsea.stats.auc.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.gsea.stats.auc.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.gsea.stats.auc.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.gsea.stats.auc.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.gsea.stats.auc.tissue)[[2]]
    cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    if((count %% nf) == 1) {
      pdf(file.path(saveres, sprintf("barplot_gsea_ic50_auc_tissue%i_cgp_ccle_paper.pdf", ceiling(count / nf))), height=14, width=10)
      par(las=1, mar=c(6, 4, 6, 0), xaxt="n", mfrow=c(3, 2))
    }
    # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
    yylim <- c(0, 1)
    mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main=sprintf("Pathway-drug associations (all)\n%s (%i cell lines)", gsub("_", " ", tissuen[ii]), ncelline))
    legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
    # plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
    axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
    text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
    text(x=mp[1, ], y=cors.ic50$rho, pos=ifelse(is.na(cors.ic50$rho) | cors.ic50$rho < 0, 1, 3), labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
    text(x=mp[2, ], y=cors.auc$rho, pos=ifelse(is.na(cors.auc$rho) | cors.auc$rho < 0, 1, 3), labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
    if((count %% nf) == 0) { dev.off() }
  }
}
if((count %% nf) != 0) { dev.off() }

## barplot for significant geneset-drug associations across tissue types
pdf(file.path(saveres, "barplot_gsea_ic50_auc_signif_tissue_cgp_ccle.pdf"), height=6, width=6)
for(ii in 1:length(tissuen)) {
  ncelline <- sum(tissue == tissuen[ii])
  if(max(correlations.gsea.stats.ic50.signif.tissue[ii, , "n"], na.rm=TRUE) >= minsample && max(correlations.gsea.stats.ic50.signif.tissue[ii, , "rho"], na.rm=TRUE) > 0) {
    ## IC50
    xx <- correlations.gsea.stats.ic50.signif.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.gsea.stats.ic50.signif.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.gsea.stats.ic50.signif.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.gsea.stats.ic50.signif.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.gsea.stats.ic50.signif.tissue)[[2]]
    cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    ## AUC
    xx <- correlations.gsea.stats.auc.signif.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.gsea.stats.auc.signif.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.gsea.stats.auc.signif.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.gsea.stats.auc.signif.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.gsea.stats.auc.signif.tissue)[[2]]
    cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    par(las=1, mar=c(6, 4, 6, 0), xaxt="n")
    # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
    yylim <- c(0, 1)
    mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main=sprintf("Pathway-drug associations (all)\n%s (%i cell lines)", gsub("_", " ", tissuen[ii]), ncelline))
    # legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
    # plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
    axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
    text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
    text(x=mp[1, ], y=cors.ic50$rho, pos=ifelse(is.na(cors.ic50$rho) | cors.ic50$rho < 0, 1, 3), labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
    text(x=mp[2, ], y=cors.auc$rho, pos=ifelse(is.na(cors.auc$rho) | cors.auc$rho < 0, 1, 3), labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
  }
}
dev.off()

## multiple barplots per page
count <- 0
nf <- 6
for(ii in 1:length(tissuen)) {
  ncelline <- sum(tissue == tissuen[ii])
  if(max(correlations.gsea.stats.ic50.signif.tissue[ii, , "n"], na.rm=TRUE) >= minsample && max(correlations.gsea.stats.ic50.signif.tissue[ii, , "rho"], na.rm=TRUE) > 0) {
    count <- count + 1
    ## IC50
    xx <- correlations.gsea.stats.ic50.signif.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.gsea.stats.ic50.signif.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.gsea.stats.ic50.signif.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.gsea.stats.ic50.signif.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.gsea.stats.ic50.signif.tissue)[[2]]
    cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    ## AUC
    xx <- correlations.gsea.stats.auc.signif.tissue[ii, , "rho"]
    xx[!is.na(xx) & xx < 0] <- 0
    ll <- correlations.gsea.stats.auc.signif.tissue[ii, , "lower"]
    ll[!is.na(ll) & ll < 0] <- 0
    uu <- correlations.gsea.stats.auc.signif.tissue[ii, , "upper"]
    uu[!is.na(uu) & uu < 0] <- 0
    uu[!is.na(uu) & uu > 1] <- 1
    pp <- correlations.gsea.stats.auc.signif.tissue[ii, , "p"]
    pp[xx == 0] <- 1
    names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.gsea.stats.auc.signif.tissue)[[2]]
    cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
    if((count %% nf) == 1) {
      pdf(file.path(saveres, sprintf("barplot_gsea_ic50_auc_signif_tissue%i_cgp_ccle_paper.pdf", ceiling(count / nf))), height=14, width=10)
      par(las=1, mar=c(6, 4, 6, 0), xaxt="n", mfrow=c(3, 2))
    }
    # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
    yylim <- c(0, 1)
    mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main=sprintf("Pathway-drug associations (all)\n%s (%i cell lines)", gsub("_", " ", tissuen[ii]), ncelline))
    legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
    # plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
    axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
    text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
    text(x=mp[1, ], y=cors.ic50$rho, pos=ifelse(is.na(cors.ic50$rho) | cors.ic50$rho < 0, 1, 3), labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
    text(x=mp[2, ], y=cors.auc$rho, pos=ifelse(is.na(cors.auc$rho) | cors.auc$rho < 0, 1, 3), labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
    if((count %% nf) == 0) { dev.off() }
  }
}
if((count %% nf) != 0) { dev.off() }

## barpots for correlations
## Spearman
for (i in 1:length(correlations.gsea.stats)) {
  fnnn <- sprintf("barplot_gsea_%s_cgp_ccle.pdf", gsub(badchars, "_", names(correlations.gsea.stats)[i]))
  pdf(file.path(saveres, fnnn), height=6, width=6)
  par(las=3, mar=c(6, 4, 3, 0), xaxt="n")
  xx <- correlations.gsea.stats[[i]][ , "rho"]
  xx[!is.na(xx) & xx < 0] <- 0
  ll <- correlations.gsea.stats[[i]][ , "lower"]
  ll[!is.na(ll) & ll < 0] <- 0
  uu <- correlations.gsea.stats[[i]][ , "upper"]
  uu[!is.na(uu) & uu < 0] <- 0
  uu[!is.na(uu) & uu > 1] <- 1
  pp <- correlations.gsea.stats[[i]][ , "p"]
  names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.gsea.stats[[i]])
  # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
  yylim <- c(0, 0.71)
  mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
  axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
  text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
  title(main=sprintf("Pathway-drug associations %s\nCGP vs. CCLE", toupper(gsub(badchars, " ", names(correlations.gsea.stats)[i]))))
  text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
  dev.off()
}

## gsea IC50 with scatter plot
pdf(file.path(saveres, "scatterbarplot_gsea_ic50_cgp_ccle_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
for(i in 1:length(gsea.ic50.ccle)) {
  myx <- complete.cases(gsea.ic50.ccle[[i]][ ,"nes"], gsea.ic50.cgp[[i]][ ,"nes"])
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(gsea.ic50.cgp[[i]][ ,"nes"], gsea.ic50.ccle[[i]][ ,"nes"], method="spearman", use="complete.obs", alternative="greater")
  }
  par(mar=c(6, 4, 3, 1) + 0.1, las=1)
  # llim <- round(range(c(gsea.ic50.cgp[[i]][ ,"nes"], gsea.ic50.ccle[[i]][ ,"nes"]), na.rm=TRUE) * 10) / 10
  llim <- c(-3, 3)
  myScatterPlot(x=gsea.ic50.cgp[[i]][ ,"nes"], y=gsea.ic50.ccle[[i]][ ,"nes"], xlab=ifelse(i > 11, "Enrichment score IC50 (CGP)", ""), ylab=ifelse((i %% 4) == 1, "Enrichment score IC50 (CCLE)", ""), main=gsub("drugid_", "", names(gsea.ic50.ccle)[i]), xlim=llim, ylim=llim, pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
}
par(mar=c(6, 4, 3, 1) + 0.1, xaxt="n", las=1)
xx <- correlations.gsea.stats[["ic50"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.gsea.stats[["ic50"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.gsea.stats[["ic50"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.gsea.stats[["ic50"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.gsea.stats[["ic50"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()

## gsea IC50 fdr with scatter plot
pdf(file.path(saveres, "scatterbarplot_gsea_ic50_filt_cgp_ccle_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
for(i in 1:length(gsea.ic50.ccle)) {
  myx <- complete.cases(gsea.ic50.ccle[[i]][ ,"fdr"], gsea.ic50.cgp[[i]][ ,"fdr"]) & (gsea.ic50.ccle[[i]][ ,"fdr"] < myfdr | gsea.ic50.cgp[[i]][ ,"fdr"] < myfdr)
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(gsea.ic50.cgp[[i]][myx, "nes"], gsea.ic50.ccle[[i]][myx, "nes"], method="spearman", use="complete.obs", alternative="greater")
  }
  par(mar=c(6, 4, 3, 1) + 0.1, las=1)
  # llim <- round(range(c(gsea.ic50.cgp[[i]][ ,"nes"], gsea.ic50.ccle[[i]][ ,"nes"]), na.rm=TRUE) * 10) / 10
  llim <- c(-3, 3)
  myScatterPlot(x=gsea.ic50.cgp[[i]][myx, "nes"], y=gsea.ic50.ccle[[i]][myx, "nes"], xlab=ifelse(i > 11, "Enrichment score IC50 (CGP)", ""), ylab=ifelse((i %% 4) == 1, "Enrichment score IC50 (CCLE)", ""), main=gsub("drugid_", "", names(gsea.ic50.ccle)[i]), xlim=llim, ylim=llim, pch=16, method="transparent", transparency=0.5)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
}
par(mar=c(6, 4, 3, 1) + 0.1, xaxt="n", las=1)
xx <- correlations.gsea.stats[["ic50.filt"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.gsea.stats[["ic50.filt"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.gsea.stats[["ic50.filt"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.gsea.stats[["ic50.filt"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.gsea.stats[["ic50.filt"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()
  
## gsea AUC with scatter plot
pdf(file.path(saveres, "scatterbarplot_gsea_auc_cgp_ccle_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
for(i in 1:length(gsea.auc.ccle)) {
  myx <- complete.cases(gsea.auc.ccle[[i]][ ,"nes"], gsea.auc.cgp[[i]][ ,"nes"])
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(gsea.auc.cgp[[i]][ ,"nes"], gsea.auc.ccle[[i]][ ,"nes"], method="spearman", use="complete.obs", alternative="greater")
  }
  par(mar=c(6, 4, 3, 1) + 0.1, las=1)
  # llim <- round(range(c(gsea.auc.cgp[[i]][ ,"nes"], gsea.auc.ccle[[i]][ ,"nes"]), na.rm=TRUE) * 10) / 10
  llim <- c(-3, 3)
  myScatterPlot(x=gsea.auc.cgp[[i]][ ,"nes"], y=gsea.auc.ccle[[i]][ ,"nes"], xlab=ifelse(i > 11, "Enrichment score AUC (CGP)", ""), ylab=ifelse((i %% 4) == 1, "Enrichment score AUC (CCLE)", ""), main=gsub("drugid_", "", names(gsea.auc.ccle)[i]), xlim=llim, ylim=llim, pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
}
par(mar=c(6, 4, 3, 1) + 0.1, xaxt="n", las=1)
xx <- correlations.gsea.stats[["auc"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.gsea.stats[["auc"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.gsea.stats[["auc"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.gsea.stats[["auc"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.gsea.stats[["auc"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()

## gsea AUC fdr with scatter plot
pdf(file.path(saveres, "scatterbarplot_gsea_auc_filt_cgp_ccle_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
for(i in 1:length(gsea.auc.ccle)) {
  myx <- complete.cases(gsea.auc.ccle[[i]][ ,"fdr"], gsea.auc.cgp[[i]][ ,"fdr"]) & (gsea.auc.ccle[[i]][ ,"fdr"] < myfdr | gsea.auc.cgp[[i]][ ,"fdr"] < myfdr)
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(gsea.auc.cgp[[i]][myx, "nes"], gsea.auc.ccle[[i]][myx, "nes"], method="spearman", use="complete.obs", alternative="greater")
  }
  par(mar=c(6, 4, 3, 1) + 0.1, las=1)
  # llim <- round(range(c(gsea.auc.cgp[[i]][ ,"nes"], gsea.auc.ccle[[i]][ ,"nes"]), na.rm=TRUE) * 10) / 10
  llim <- c(-3, 3)
  myScatterPlot(x=gsea.auc.cgp[[i]][myx, "nes"], y=gsea.auc.ccle[[i]][myx, "nes"], xlab=ifelse(i > 11, "Enrichment score AUC (CGP)", ""), ylab=ifelse((i %% 4) == 1, "Enrichment score AUC (CCLE)", ""), main=gsub("drugid_", "", names(gsea.auc.ccle)[i]), xlim=llim, ylim=llim, pch=16, method="transparent", transparency=0.5)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
}
par(mar=c(6, 4, 3, 1) + 0.1, xaxt="n", las=1)
xx <- correlations.gsea.stats[["auc.filt"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.gsea.stats[["auc.filt"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.gsea.stats[["auc.filt"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.gsea.stats[["auc.filt"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.gsea.stats[["auc.filt"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()



## statistically compare correlations
## IC50 vs filetered IC50
message("IC50 vs IC50 filtered by max concentration:")
wilcox.test(x=correlations[["ic50"]][ ,"drug.sensitivity.filt"], y=correlations[["ic50"]][ ,"drug.sensitivity"], alternative="two.sided", paired=TRUE)
message("")
## AUC vs IC50
message("AUC vs IC50:")
wilcox.test(x=correlations[["auc"]][ ,"drug.sensitivity"], y=correlations[["ic50"]][ ,"drug.sensitivity"], alternative="two.sided", paired=TRUE)
message("")
## AUC vs IC50 calls
message("AUC vs IC50 calls:")
wilcox.test(x=correlations[["auc.call"]][ ,"drug.sensitivity"], y=correlations[["ic50.call"]][ ,"drug.sensitivity"], alternative="two.sided", paired=TRUE)
message("")
## gene-drug associations, AUC vs IC50
message("AUC vs IC50 for gene-drug associations:")
wilcox.test(x=correlations[["auc"]][ ,"gene.drug"], y=correlations[["ic50"]][ ,"gene.drug"], alternative="two.sided", paired=TRUE)
message("")
## gene-drug associations fdr, AUC vs IC50
message("AUC vs IC50 for significant (FDR) gene-drug associations")
wilcox.test(x=correlations[["auc"]][ ,"gene.drug.filt"], y=correlations[["ic50"]][ ,"gene.drug.filt"], alternative="two.sided", paired=TRUE)
message("")
## geneset-drug associations, AUC vs IC50
message("AUC vs IC50 for geneset-drug associations")
wilcox.test(x=correlations[["auc"]][ ,"go.drug"], y=correlations[["ic50"]][ ,"go.drug"], alternative="two.sided", paired=TRUE)
message("")
## geneset-drug associations fdr, AUC vs IC50
message("AUC vs IC50 for significant (FDR) geneset-drug associations")
wilcox.test(x=correlations[["auc"]][ ,"go.drug.filt"], y=correlations[["ic50"]][ ,"go.drug.filt"], alternative="two.sided", paired=TRUE)
message("")

## generate LaTeX code for supplementary tables
## IC50
xt <- correlations[["ic50"]][ , -2]
colnames(xt) <- c("(A) IC50 measures", "(B) Gene-drug associations", sprintf("(C) Significant (FDR < %i%%) gene-drug associations", round(myfdr * 100)), "(D) Pathway-drug associations", sprintf("(E) Significant (FDR < %i%%) pathway-drug associations", round(myfdr * 100)))
xt <- xtable::xtable(xt, digits=2, align=c("r", rep("p{2cm}", ncol(xt))))
xtable::print.xtable(x=xt, file=file.path(saveres, "correlations_ic50.tex"), append=FALSE)
## p-values
xt <- cbind(correlations.stats[["ic50"]][ , "rho"], correlations.assoc.stats[["ic50"]][ , "rho"], correlations.assoc.stats[["ic50.filt"]][ , "rho"], correlations.gsea.stats[["ic50"]][ , "rho"], correlations.gsea.stats[["ic50.filt"]][ , "rho"])
tt <- cbind(correlations.stats[["ic50"]][ , "p"], correlations.assoc.stats[["ic50"]][ , "p"], correlations.assoc.stats[["ic50.filt"]][ , "p"], correlations.gsea.stats[["ic50"]][ , "p"], correlations.gsea.stats[["ic50.filt"]][ , "p"])
colnames(xt) <- c("(A) IC50 measures", "(B) Gene-drug associations", sprintf("(C) Significant (FDR < %i%%) gene-drug associations", round(myfdr * 100)), "(D) Pathway-drug associations", sprintf("(E) Significant (FDR < %i%%) pathway-drug associations", round(myfdr * 100)))
xt2 <- tt <- matrix(p.adjust(tt, method="none"), ncol=ncol(tt), nrow=nrow(tt), dimnames=dimnames(tt))
xt2[!is.na(tt) & tt < 0.001] <- "***"
xt2[!is.na(tt) & tt >= 0.001 & tt < 0.01] <- "**"
xt2[!is.na(tt) & tt >= 0.01 & tt < 0.05] <- "*"
xt2[!is.na(tt) & tt >= 0.05] <- "[NS]"
xt3 <- matrix(paste(round(xt, digits=2), xt2), ncol=ncol(xt), nrow=nrow(xt), dimnames=dimnames(xt))
xt3 <- gsub("NA", "", xt3)
xt3 <- xtable::xtable(xt3, digits=0, align=c("r", rep("p{2cm}", ncol(xt3))))
xtable::print.xtable(x=xt3, file=file.path(saveres, "correlations_ic50_pvalue_paper.tex"), append=FALSE)
## AUC
xt <- correlations[["auc"]][ , -2]
colnames(xt) <- c("(A) AUC measures", "(B) Gene-drug associations", sprintf("(C) Significant (FDR < %i%%) gene-drug associations", round(myfdr * 100)), "(D) Pathway-drug associations", sprintf("(E) Significant (FDR < %i%%) pathway-drug associations", round(myfdr * 100)))
xt <- xtable::xtable(xt, digits=2, align=c("r", rep("p{2cm}", ncol(xt))))
xtable::print.xtable(x=xt, file=file.path(saveres, "correlations_auc_paper.tex"), append=FALSE)
## p-values
xt <- cbind(correlations.stats[["auc"]][ , "rho"], correlations.assoc.stats[["auc"]][ , "rho"], correlations.assoc.stats[["auc.filt"]][ , "rho"], correlations.gsea.stats[["auc"]][ , "rho"], correlations.gsea.stats[["auc.filt"]][ , "rho"])
tt <- cbind(correlations.stats[["auc"]][ , "p"], correlations.assoc.stats[["auc"]][ , "p"], correlations.assoc.stats[["auc.filt"]][ , "p"], correlations.gsea.stats[["auc"]][ , "p"], correlations.gsea.stats[["auc.filt"]][ , "p"])
colnames(xt) <- c("(A) AUC measures", "(B) Gene-drug associations", sprintf("(C) Significant (FDR < %i%%) gene-drug associations", round(myfdr * 100)), "(D) Pathway-drug associations", sprintf("(E) Significant (FDR < %i%%) pathway-drug associations", round(myfdr * 100)))
xt2 <- tt <- matrix(p.adjust(tt, method="none"), ncol=ncol(tt), nrow=nrow(tt), dimnames=dimnames(tt))
xt2[!is.na(tt) & tt < 0.001] <- "***"
xt2[!is.na(tt) & tt >= 0.001 & tt < 0.01] <- "**"
xt2[!is.na(tt) & tt >= 0.01 & tt < 0.05] <- "*"
xt2[!is.na(tt) & tt >= 0.05] <- "[NS]"
xt3 <- matrix(paste(round(xt, digits=2), xt2), ncol=ncol(xt), nrow=nrow(xt), dimnames=dimnames(xt))
xt3 <- gsub("NA", "", xt3)
xt3 <- xtable::xtable(xt3, digits=0, align=c("r", rep("p{2cm}", ncol(xt3))))
xtable::print.xtable(x=xt3, file=file.path(saveres, "correlations_auc_pvalue_paper.tex"), append=FALSE)



########################






# ## drug sensitivity measures per tissue types
# pdf(file.path(saveres, "barplot_ic50_auc_tissue_cgp_ccle_paper2.pdf"), height=6, width=6)
# for(ii in 1:length(tissuen)) {
#   ncelline <- sum(tissue == tissuen[ii])
#   if(max(correlations.stats.ic50.tissue[ii, , "n"], na.rm=TRUE) >= minsample) {
#     ## IC50
#     xx <- correlations.stats.ic50.tissue[ii, , "rho"]
#     xx[!is.na(xx) & xx < 0] <- 0
#     ll <- correlations.stats.ic50.tissue[ii, , "lower"]
#     ll[!is.na(ll) & ll < 0] <- 0
#     uu <- correlations.stats.ic50.tissue[ii, , "upper"]
#     uu[!is.na(uu) & uu < 0] <- 0
#     uu[!is.na(uu) & uu > 1] <- 1
#     pp <- correlations.stats.ic50.tissue[ii, , "p"]
#     pp[xx == 0] <- 1
#     names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.stats.ic50.tissue)[[2]]
#     cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
#     ## AUC
#     xx <- correlations.stats.auc.tissue[ii, , "rho"]
#     xx[!is.na(xx) & xx < 0] <- 0
#     ll <- correlations.stats.auc.tissue[ii, , "lower"]
#     ll[!is.na(ll) & ll < 0] <- 0
#     uu <- correlations.stats.auc.tissue[ii, , "upper"]
#     uu[!is.na(uu) & uu < 0] <- 0
#     uu[!is.na(uu) & uu > 1] <- 1
#     pp <- correlations.stats.auc.tissue[ii, , "p"]
#     pp[xx == 0] <- 1
#     names(xx) <- names(ll) <- names(uu) <- names(pp) <- dimnames(correlations.stats.auc.tissue)[[2]]
#     cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
#     par(las=1, mar=c(6, 4, 6, 0), xaxt="n")
#     # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
#     yylim <- c(0, 1)
#     mycolls <- rep(rainbow(length(xx), v=0.9)
#     
#     mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=mycolls, each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main=sprintf("Drug sensitivity measures\n%s (%i cell lines)", gsub("_", " ", tissuen[ii]), ncelline))
#     # legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
#     # plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
#     axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
#     text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
#     text(x=mp[1, ], y=cors.ic50$rho, pos=ifelse(is.na(cors.ic50$rho) | cors.ic50$rho < 0, 1, 3), labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
#     text(x=mp[2, ], y=cors.auc$rho, pos=ifelse(is.na(cors.auc$rho) | cors.auc$rho < 0, 1, 3), labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5, offset=0)
#   }
# }
# dev.off()


## load results for all tissues combined
load(file.path(saveres, "correlations_stats.RData"))
load(file.path(saveres, "correlations_assoc_stats.RData"))
load(file.path(saveres, "correlations_gsea_stats.RData"))


pdf(file.path(saveres, "barplot_ic50_auc_call_cgp_ccle_paper.pdf"), height=6, width=6)
par(las=3, mar=c(6, 4, 1, 0), xaxt="n")
## IC50
xx <- correlations.stats[["ic50.call"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.stats[["ic50.call"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.stats[["ic50.call"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.stats[["ic50.call"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.stats[["ic50.call"]])
cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
## AUC
xx <- correlations.stats[["auc.call"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.stats[["auc.call"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.stats[["auc.call"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.stats[["auc.call"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.stats[["auc.call"]])
cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.5)
# mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Kappa", ylim=yylim, angle=c(45, -45), density=c(100, 40), legend.text=c("IC50 calls", "AUC calls"), args.legend=list(x="topleft", col=c("darkgrey", "darkgrey"), bty="n"), main="Drug sensitivity calling")
  mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Kappa", ylim=yylim, angle=c(45, -45), density=c(100, 40), main="Drug sensitivity calling")
  legend("topleft", legend=c("IC50 calls", "AUC calls"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
# plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
text(x=mp[1, ] + 1.45, y=rbind(cors.ic50$rho, cors.auc$rho)[1, ], pos=2, labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
text(x=mp[2, ] + 1.45, y=rbind(cors.ic50$rho, cors.auc$rho)[2, ], pos=2, labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()

# legend("topleft", legend=c("IC50 calls", "AUC calls", "rgerg", "wrgetg"), fill=rep("red", 4), density=c(200, 40, 25, 10), bty="n", cex=2)

## gene/geneset-drug associations
pdf(file.path(saveres, "barplot_assoc_ic50_auc_cgp_ccle_paper.pdf"), height=12, width=12)
par(mfrow=c(2, 2))

## all gene-drug associations
par(las=3, mar=c(6, 4, 4, 1), xaxt="n")
yylim <- c(0, 0.76)
## IC50
xx <- correlations.assoc.stats[["ic50"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.assoc.stats[["ic50"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.assoc.stats[["ic50"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.assoc.stats[["ic50"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.assoc.stats[["ic50"]])
cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
## AUC
xx <- correlations.assoc.stats[["auc"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.assoc.stats[["auc"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.assoc.stats[["auc"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.assoc.stats[["auc"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.assoc.stats[["auc"]])
cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main="Gene-drug associations (all)")
legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
# plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
text(x=mp[1, ] + 1.45, y=rbind(cors.ic50$rho, cors.auc$rho)[1, ], pos=2, labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
text(x=mp[2, ] + 1.45, y=rbind(cors.ic50$rho, cors.auc$rho)[2, ], pos=2, labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)

## significant gene-drug associations
par(las=3, mar=c(6, 4, 4, 1), xaxt="n")
## IC50
xx <- correlations.assoc.stats[["ic50.filt"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.assoc.stats[["ic50.filt"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.assoc.stats[["ic50.filt"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.assoc.stats[["ic50.filt"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.assoc.stats[["ic50.filt"]])
cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
## AUC
xx <- correlations.assoc.stats[["auc.filt"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.assoc.stats[["auc.filt"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.assoc.stats[["auc.filt"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.assoc.stats[["auc.filt"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.assoc.stats[["auc.filt"]])
cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main=sprintf("Gene-drug associations (FDR < %i%%)", round(myfdr * 100)))
legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
# plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
text(x=mp[1, ] + 1.45, y=rbind(cors.ic50$rho, cors.auc$rho)[1, ], pos=2, labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
text(x=mp[2, ] + 1.45, y=rbind(cors.ic50$rho, cors.auc$rho)[2, ], pos=2, labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)

## all geneset-drug associations
par(las=3, mar=c(6, 4, 4, 1), xaxt="n")
## IC50
xx <- correlations.gsea.stats[["ic50"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.gsea.stats[["ic50"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.gsea.stats[["ic50"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.gsea.stats[["ic50"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.gsea.stats[["ic50"]])
cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
## AUC
xx <- correlations.gsea.stats[["auc"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.gsea.stats[["auc"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.gsea.stats[["auc"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.gsea.stats[["auc"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.gsea.stats[["auc"]])
cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main="Pathway-drug associations (all)")
legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
# plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
text(x=mp[1, ] + 1.45, y=rbind(cors.ic50$rho, cors.auc$rho)[1, ], pos=2, labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
text(x=mp[2, ] + 1.45, y=rbind(cors.ic50$rho, cors.auc$rho)[2, ], pos=2, labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)

## significant geneset-drug associations
par(las=3, mar=c(6, 4, 4, 1), xaxt="n")
## IC50
xx <- correlations.gsea.stats[["ic50.filt"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.gsea.stats[["ic50.filt"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.gsea.stats[["ic50.filt"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.gsea.stats[["ic50.filt"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.gsea.stats[["ic50.filt"]])
cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
## AUC
xx <- correlations.assoc.stats[["auc.filt"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.assoc.stats[["auc.filt"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.assoc.stats[["auc.filt"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.assoc.stats[["auc.filt"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.assoc.stats[["auc.filt"]])
cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main=sprintf("Pathway-drug associations (FDR < %i%%)", round(myfdr * 100)))
legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
# plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
text(x=mp[1, ] + 1.45, y=rbind(cors.ic50$rho, cors.auc$rho)[1, ], pos=2, labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
text(x=mp[2, ] + 1.45, y=rbind(cors.ic50$rho, cors.auc$rho)[2, ], pos=2, labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()


## save all correlations
savecor <- NULL
for(i in 1:length(correlations)) {
  tt <- correlations[[i]][ , apply(correlations[[i]], 2, function(x) { return(!all(is.na(x))) }), drop=FALSE]
  savecor <- c(savecor, list(data.frame(tt)))
}
names(savecor) <- names(correlations)
WriteXLS::WriteXLS("savecor", ExcelFileName=file.path(saveres, "correlations.xls"), row.names=TRUE)



###############################################################################
###############################################################################
## Correlation of mutation-drug associations
###############################################################################
###############################################################################

## binarize the mutation matrix as in CGP and CCLE
mutation.cgp <- (mutation.cgp != "wt") * 1
mutation.ccle <- (mutation.ccle != "wt") * 1

## correlation statistics
tt <- matrix(NA, nrow=length(drugsn), ncol=5, dimnames=list(drugsn, c("rho", "lower", "upper", "p", "n")))
correlations.assoc.mut.stats <- list("ic50"=tt, "auc"=tt)


## IC50 all tissues
myfn <- file.path(saveres, "assoc_mut_ic50_cgp_ccle.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Mutation-drug association based on IC50 (CCLE)")
  assoc.mut.ic50.ccle <- NULL
  for(i in 1:ncol(ic50.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(ic50.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(mutation.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, ic50, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=ic50, z=tissue, method=genedrugm)
      return(res)      
    }, data=mutation.ccle, ic50=ic50.ccle[ ,i], tissue=tissue.ccle)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(mutation.ccle), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.mut.ic50.ccle <- c(assoc.mut.ic50.ccle, list(mcres))
  }
  message("")
  names(assoc.mut.ic50.ccle) <- colnames(ic50.ccle)
  ## CGP
  message("Gene-drug association based on IC50 (CGP)")
  assoc.mut.ic50.cgp <- NULL
  for(i in 1:ncol(ic50.cgp)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(ic50.cgp)[i]))
    splitix <- parallel::splitIndices(nx=ncol(mutation.cgp), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, ic50, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=ic50, z=tissue, method=genedrugm)
      return(res)      
    }, data=mutation.cgp, ic50=ic50.cgp[ ,i], tissue=tissue.cgp)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(mutation.cgp), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.mut.ic50.cgp <- c(assoc.mut.ic50.cgp, list(mcres))
  }
  message("")
  names(assoc.mut.ic50.cgp) <- colnames(ic50.cgp)
  ## save all associations
  ## CCLE
  rr <- NULL
  for(i in 1:length(assoc.mut.ic50.ccle)) {
    tt <- cbind(assoc.mut.ic50.ccle[[i]])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.mut.ic50.ccle)
  WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "ic50_ccle_results_mutation_drug_paper.xls"), row.names=TRUE)
  ## CGP
  rr <- NULL
  for(i in 1:length(assoc.mut.ic50.cgp)) {
    tt <- cbind(assoc.mut.ic50.cgp[[i]])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.mut.ic50.cgp)
  WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "ic50_cgp_results_mutation_drug_paper.xls"), row.names=TRUE)
  save(list=c("assoc.mut.ic50.cgp", "assoc.mut.ic50.ccle"), compress=TRUE, file=myfn)
} else { load(myfn) }

## consistency between associations
## all genes
## IC50
message("Correlation mutation-drug associations based on IC50, all genes:")
for(i in 1:length(assoc.mut.ic50.ccle)) {
  tt <- min(round(quantile(abs(c(assoc.mut.ic50.cgp[[i]][ , "estimate"], assoc.mut.ic50.ccle[[i]][ , "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 5)
  xxlim <- yylim <- c(-tt, tt)
  myx <- complete.cases(assoc.mut.ic50.ccle[[i]][ ,"estimate"], assoc.mut.ic50.cgp[[i]][ ,"estimate"])
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(assoc.mut.ic50.cgp[[i]][ ,"estimate"], assoc.mut.ic50.ccle[[i]][ ,"estimate"], method="spearman", use="complete.obs", alternative="greater")
  }
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=sum(myx), alpha=0.05)
  correlations.assoc.mut.stats[["ic50"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, sum(myx))
}
message("")

## AUC all tissues
myfn <- file.path(saveres, "assoc_mut_auc_cgp_ccle.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Mutation-drug association based on AUC (CCLE)")
  assoc.mut.auc.ccle <- NULL
  for(i in 1:ncol(auc.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(mutation.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm)
      return(res)      
    }, data=mutation.ccle, auc=auc.ccle[ ,i], tissue=tissue.ccle)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(mutation.ccle), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.mut.auc.ccle <- c(assoc.mut.auc.ccle, list(mcres))
  }
  message("")
  names(assoc.mut.auc.ccle) <- colnames(auc.ccle)
  ## CGP
  message("Gene-drug association based on AUC (CGP)")
  assoc.mut.auc.cgp <- NULL
  for(i in 1:ncol(auc.cgp)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.cgp)[i]))
    splitix <- parallel::splitIndices(nx=ncol(mutation.cgp), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      res <- apply(X=data[ ,x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm)
      return(res)      
    }, data=mutation.cgp, auc=auc.cgp[ ,i], tissue=tissue.cgp)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(mutation.cgp), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.mut.auc.cgp <- c(assoc.mut.auc.cgp, list(mcres))
  }
  message("")
  names(assoc.mut.auc.cgp) <- colnames(auc.cgp)
  ## save all associations
  ## CCLE
  rr <- NULL
  for(i in 1:length(assoc.mut.auc.ccle)) {
    tt <- cbind(assoc.mut.auc.ccle[[i]])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.mut.auc.ccle)
  WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "auc_ccle_results_mutation_drug_paper.xls"), row.names=TRUE)
  ## CGP
  rr <- NULL
  for(i in 1:length(assoc.mut.auc.cgp)) {
    tt <- cbind(assoc.mut.auc.cgp[[i]])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.mut.auc.cgp)
  WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "auc_cgp_results_mutation_drug_paper.xls"), row.names=TRUE)
  save(list=c("assoc.mut.auc.cgp", "assoc.mut.auc.ccle"), compress=TRUE, file=myfn)
} else { load(myfn) }

## consistency between associations
## all genes
## AUC
message("Correlation mutation-drug associations based on AUC, all genes:")
for(i in 1:length(assoc.mut.auc.ccle)) {
  tt <- min(round(quantile(abs(c(assoc.mut.auc.cgp[[i]][ , "estimate"], assoc.mut.auc.ccle[[i]][ , "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 5)
  xxlim <- yylim <- c(-tt, tt)
  myx <- complete.cases(assoc.mut.auc.ccle[[i]][ ,"estimate"], assoc.mut.auc.cgp[[i]][ ,"estimate"])
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(assoc.mut.auc.cgp[[i]][ ,"estimate"], assoc.mut.auc.ccle[[i]][ ,"estimate"], method="spearman", use="complete.obs", alternative="greater")
  }
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=sum(myx), alpha=0.05)
  correlations.assoc.mut.stats[["auc"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, sum(myx))
}
message("")

pdf(file.path(saveres, "barplot_assoc_mut_ic50_auc_cgp_ccle_paper.pdf"), height=6, width=6)
## all mutation-drug associations
par(las=3, mar=c(6, 4, 4, 1), xaxt="n")
yylim <- c(0, 0.76)
## IC50
xx <- correlations.assoc.mut.stats[["ic50"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.assoc.mut.stats[["ic50"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.assoc.mut.stats[["ic50"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.assoc.mut.stats[["ic50"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.assoc.mut.stats[["ic50"]])
cors.ic50 <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
## AUC
xx <- correlations.assoc.mut.stats[["auc"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.assoc.mut.stats[["auc"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.assoc.mut.stats[["auc"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.assoc.mut.stats[["auc"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.assoc.mut.stats[["auc"]])
cors.auc <- list("rho"=xx, "lower"=ll, "upper"=uu, "p"=pp)
mp <- barplot(height=rbind(cors.ic50$rho, cors.auc$rho), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab="Rs", ylim=yylim, angle=c(45, -45), density=c(100, 40), main="Mutation-drug associations (all)")
legend("topleft", legend=c("IC50", "AUC"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
# plotrix::plotCI(x=mp, y=rbind(cors.ic50$rho, cors.auc$rho), li=rbind(cors.ic50$lower, cors.auc$lower), ui=rbind(cors.ic50$upper, cors.auc$upper), err="y", pch=".", add=TRUE, sfrac=0.0025)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
text(x=mp[1, ] + 1.75, y=rbind(cors.ic50$rho, cors.auc$rho)[1, ], pos=2, labels=ifelse(cors.ic50$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
text(x=mp[2, ] + 1.75, y=rbind(cors.ic50$rho, cors.auc$rho)[2, ], pos=2, labels=ifelse(cors.auc$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()

## mutation-drug associations IC50 with scatter plot
pdf(file.path(saveres, "scatterbarplot_assoc_mut_ic50_cgp_ccle_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
for(i in 1:length(assoc.mut.ic50.ccle)) {
  tt <- min(round(quantile(abs(c(assoc.mut.ic50.cgp[[i]][ , "estimate"], assoc.mut.ic50.ccle[[i]][ , "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 5)
  xxlim <- yylim <- c(-tt, tt)
  myx <- complete.cases(assoc.mut.ic50.ccle[[i]][ ,"estimate"], assoc.mut.ic50.cgp[[i]][ ,"estimate"])
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(assoc.mut.ic50.cgp[[i]][ ,"estimate"], assoc.mut.ic50.ccle[[i]][ ,"estimate"], method="spearman", use="complete.obs", alternative="greater")
  }
  par(mar=c(6, 4, 3, 1) + 0.1, las=1)
  myScatterPlot(x=assoc.mut.ic50.cgp[[i]][ ,"estimate"], y=assoc.mut.ic50.ccle[[i]][ ,"estimate"], xlab=ifelse(i > 11, "Mutation-drug association IC50 (CGP)", ""), ylab=ifelse((i %% 4) == 1, "Mutation-drug association IC50 (CCLE)", ""), main=gsub("drugid_", "", names(assoc.mut.ic50.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
}
par(mar=c(6, 4, 3, 1) + 0.1, xaxt="n", las=1)
xx <- correlations.assoc.mut.stats[["ic50"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.assoc.mut.stats[["ic50"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.assoc.mut.stats[["ic50"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.assoc.mut.stats[["ic50"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.assoc.mut.stats[["ic50"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()


## mutation-drug associations AUC with scatter plot
pdf(file.path(saveres, "scatterbarplot_assoc_mut_auc_cgp_ccle_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
for(i in 1:length(assoc.mut.auc.ccle)) {
  tt <- min(round(quantile(abs(c(assoc.mut.auc.cgp[[i]][ , "estimate"], assoc.mut.auc.ccle[[i]][ , "estimate"])), probs=0.999, na.rm=TRUE) * 10) / 10, 5)
  xxlim <- yylim <- c(-tt, tt)
  myx <- complete.cases(assoc.mut.auc.ccle[[i]][ ,"estimate"], assoc.mut.auc.cgp[[i]][ ,"estimate"])
  if(sum(myx) < minsample) {
   cc <- list("estimate"=NA, "p.value"=NA) 
  } else {
    cc <- cor.test(assoc.mut.auc.cgp[[i]][ ,"estimate"], assoc.mut.auc.ccle[[i]][ ,"estimate"], method="spearman", use="complete.obs", alternative="greater")
  }
  par(mar=c(6, 4, 3, 1) + 0.1, las=1)
  myScatterPlot(x=assoc.mut.auc.cgp[[i]][ ,"estimate"], y=assoc.mut.auc.ccle[[i]][ ,"estimate"], xlab=ifelse(i > 11, "Mutation-drug association AUC (CGP)", ""), ylab=ifelse((i %% 4) == 1, "Mutation-drug association AUC (CCLE)", ""), main=gsub("drugid_", "", names(assoc.mut.auc.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.3)
  abline(v=0, col="red", lty=2, lwd=0.25)
  abline(h=0, col="red", lty=2, lwd=0.25)
}
par(mar=c(6, 4, 3, 1) + 0.1, xaxt="n", las=1)
xx <- correlations.assoc.mut.stats[["auc"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.assoc.mut.stats[["auc"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.assoc.mut.stats[["auc"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.assoc.mut.stats[["auc"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.assoc.mut.stats[["auc"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()







## end


