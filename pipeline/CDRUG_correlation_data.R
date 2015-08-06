########################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## February 11, 2013
########################


require(amap) || stop("Library amap is not available!")
require(vcd) || stop("Library vcd is not available!")
require(gplots) || stop("Library gplots is not available!")
require(WriteXLS) || stop("Library WriteXLS is not available!")
require(xtable) || stop("Library xtable is not available!")

################################################
## correlation across replicated cell lines
################################################

myfn <- file.path(saveres, sprintf("ge_var%i_cellines_correlations_replicates_cgp.RData", topvar))
if(!file.exists(myfn)) {
  ## load full cgp dataset (with replicates)
  load(file.path(saveres, "cgp_full_frma.RData"))
  ## gene centric data
  myx <- which(annot.full.cgp[ ,"best"])
  myx <- myx[!duplicated(annot.full.cgp[myx, "jetset.EntrezID"])]
  data.full.cgp <- data.full.cgp[ , myx,drop=FALSE]
  annot.full.cgp <- annot.full.cgp[myx, ,drop=FALSE]
  colnames(data.full.cgp) <- rownames(annot.full.cgp) <- paste("geneid", annot.full.cgp[ , "jetset.EntrezID"], sep="_")
  
  ## select the top most variant genes
  varg <- unique(c(colnames(data.full.cgp)[order(apply(data.full.cgp, 2, var, na.rm=TRUE), decreasing=TRUE)[1:topvar]], colnames(data.full.cgp)[order(apply(data.full.cgp, 2, var, na.rm=TRUE), decreasing=TRUE)[1:topvar]]))

  ## list of duplicated cell lines
  duplix <- duplicated(sampleinfo.full.cgp[ , "cellid"])
  dcelln <- as.character(unique(sampleinfo.full.cgp[duplix, "cellid"]))
  cor.reps <- NULL
  for(i in 1:length(dcelln)) {
    iix <- which(!is.na(sampleinfo.full.cgp[ , "cellid"]) & sampleinfo.full.cgp[ , "cellid"] == dcelln[i])
    iixn <- combn(length(iix), 2)
    ccor <- apply(iixn, 2, function(x, y, z) {
      return(cor(x=z[y[x[1]], ], y=z[y[x[2]], ], method="spearman", use="complete.obs"))
    }, y=iix, z=data.full.cgp[ , varg, drop=FALSE])
    cor.reps <- c(cor.reps, list(ccor))
  }
  names(cor.reps) <- dcelln

  ## list of duplicated cell lines
  cor.diffs <- cor(x=t(data.full.cgp[!duplix, varg, drop=FALSE]), method="spearman", use="pairw")
  diag(cor.diffs) <- NA
  
  ## boxplot
  pdf(file.path(saveres, sprintf("ge_var%i_cellines_boxplot_replicates_cgp_paper.pdf", topvar)), height=7, width=7)
  ll <- list("Replicates"=unlist(cor.reps), "Different"=as.vector(cor.diffs[lower.tri(cor.diffs)]))
  wt <- wilcox.test(x=ll[[1]], y=ll[[2]])
  par(xaxt="s", yaxt="s", mar=c(5, 4, 4, 2) + 0.1)
  boxplot(ll, outline=FALSE, ylab="Rs", xlab="Cell lines", main="Gene expression profiles across cell lines\nCGP", border="black", col="lightgrey", ylim=c(-0.10, 1))
  legend("bottomleft", legend=sprintf("Wilcoxon test p-value=%.1E", wt$p.value), bty="n")
  dev.off()
  message("")
  save(list=c("cor.diffs", "cor.reps"), compress=TRUE, file=myfn)
} else { load(myfn) }

################################################
## correlation across cell lines in CGP and CCLE
################################################

message("Correlation of gene expressions across cell lines (CGP vs CCLE)")

load(file.path(saveres, "CDRUG_cgp_ccle.RData"))
tissue.cgp <- sampleinfo.cgp[ ,"tissue.type"]
tissue.ccle <- sampleinfo.ccle[ ,"tissue.type"]

myfn <- file.path(saveres, sprintf("ge_var%i_cellines_correlations.RData", topvar))
if(!file.exists(myfn)) {
  ## the top 1000 most variant genes in either CCLE or CGP
  varg <- unique(c(colnames(data.ccle)[order(apply(data.ccle, 2, var, na.rm=TRUE), decreasing=TRUE)[1:topvar]], colnames(data.cgp)[order(apply(data.cgp, 2, var, na.rm=TRUE), decreasing=TRUE)[1:topvar]]))
  cor.all <- cor(x=t(data.cgp[ , varg, drop=FALSE]), y=t(data.ccle[ , varg, drop=FALSE]), method="spearman", use="pairw")
  save(list="cor.all", compress=TRUE, file=myfn)
} else { load(myfn) }
oo <- amap::hcluster(cor.all, link="complete")$order
pdf(file.path(saveres, sprintf("ge_var%i_cellines_boxplot_ccle_cgp_paper.pdf", topvar)), height=7, width=7)
par(xaxt="n", yaxt="n")
image(cor.all[oo, rev(oo)], col=blueyellow(255), xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="Correlation of gene expression profiles across cell lines\nCCLE vs. CGP")
dev.off()
pdf(file.path(saveres, sprintf("ge_var%i_cellines_heatmap2_ccle_cgp_paper.pdf", topvar)), height=8, width=8)
heatmap.3(x=cor.all[oo, oo], Colv=FALSE, Rowv=FALSE, scale="none", dendrogram="none", col=blueyellow(255), trace="none", density.info="none", keysize=1.5, xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="Gene expressions across cell lines\nCCLE vs. CGP", labRow=NA, labCol=NA, key=TRUE, KeyValueName="", keyTitle="Rs", breaks=seq(0, 1, length.out=256), mar=c(2, 2) + 0.1, na.color="lightgray")
dev.off()


cor.allVar <- diag(cor.all)
median(cor.allVar)
cor.allALL <- cor(x=t(data.cgp[ , , drop=FALSE]), y=t(data.ccle[ , , drop=FALSE]), method="spearman", use="pairw")
cor.all_genesVar <- cor(x=data.cgp[ , varg, drop=FALSE], y=data.ccle[ , varg, drop=FALSE], method="spearman", use="pairw")
cor.all_genesVarDiag <- diag(cor.all_genesVar)
cor.all_genesNoVar <- cor(x=data.cgp[ , , drop=FALSE], y=data.ccle[ , , drop=FALSE], method="spearman", use="pairw")


######################################################################
############## Begin Geeleher et. al Code ############################
######################################################################

## *Across* samples Spearman correlations for all genes (gene expression data)
geneCorsAllGenes <- numeric()
for(i in 1:ncol(data.cgp))
{
  geneCorsAllGenes[i] <- cor(data.cgp[, i], data.ccle[, i], method="spearman")
}
median(geneCorsAllGenes)
# Result:
# [1] 0.5608107

## *Between* samples Spearman correlations for all genes.
sampCors <- numeric()
for(i in 1:nrow(data.cgp))
{
  sampCors[i] <- cor(data.cgp[i, ], data.ccle[i, ], method="spearman")
}
median(sampCors)
# Result:
# [1] 0.8805196


####################################################################################
##################### End Geeleher et. al Code #####################################
####################################################################################




## are the identical cell lines the mostly correlated?
cc <- sapply(1:nrow(cor.all), function(x, d) {
  return(sum(!is.na(d[x, -x]) & d[x, x] <= d[x, -x]))
}, d=cor.all)
names(cc) <- rownames(cor.all)
## cell lines which do not cocor.allelate most with their counterpart
message("These cell lines are not the most cocor.allelated with their counterpart (CCLE vs. CGP):\n", paste(names(cc)[!is.na(cc) & cc > 0], collapse="\t"))
## boxplot
pdf(file.path(saveres, sprintf("ge_var%i_cellines_boxplot_ccle_cgp_paper.pdf", topvar)), height=7, width=7)
ll <- list("Identical"=diag(cor.all), "Different"=cor.all[lower.tri(cor.all)])
wt <- wilcox.test(x=ll[[1]], y=ll[[2]])
par(xaxt="s", yaxt="s", mar=c(5, 4, 4, 2) + 0.1)
boxplot(ll, outline=FALSE, ylab="Rs", xlab="Cell lines", main="Gene expression profiles across cell lines\nCCLE vs. CGP", border="black", col="lightgrey", ylim=c(-0.10, 1))
legend("bottomleft", legend=sprintf("Wilcoxon test p-value=%.1E", wt$p.value), bty="n")
dev.off()
## with replicates in CGP
pdf(file.path(saveres, sprintf("ge_var%i_cellines_boxplot_ccle_cgp_replicates_paper.pdf", topvar)), height=7, width=7)
ll <- list("Replicates"=unlist(cor.reps), "Identical"=diag(cor.all), "Different"=cor.all[lower.tri(cor.all)])
wt0 <- wilcox.test(x=ll[[1]], y=ll[[2]])
wt1 <- wilcox.test(x=ll[[1]], y=ll[[3]])
wt2 <- wilcox.test(x=ll[[2]], y=ll[[3]])
par(xaxt="s", yaxt="s", mar=c(5, 4, 4, 2) + 0.1)
boxplot(ll, outline=FALSE, ylab="Rs", xlab="Cell lines", main="Gene expression profiles across cell lines\nReplicates in CGP + CCLE vs. CGP", border="black", col="lightgrey", ylim=c(-0.10, 1), las=1)
# legend("bottomleft", legend=sprintf("Wilcoxon test p-value=%.1E", wt$p.value), bty="n")
dev.off()
message("")

## per tissue type
message("Correlation of gene expressions across cell lines per tissue type (CGP vs CCLE)")

utissue <- table(as.character(tissue.cgp))
coltissue <- coltissuet <- factor(names(utissue))
# levels(coltissue) <- gplots::rich.colors(length(utissue))
levels(coltissue) <- rainbow(n=length(utissue), s=0.5, v=1)
levels(coltissuet) <- rainbow(n=length(utissue), s=0.75, v=1)
coltissue <- as.character(coltissue)
coltissuet <- as.character(coltissuet)
names(coltissue) <- names(coltissuet) <- names(utissue)

myfn <- file.path(saveres, sprintf("ge_var%i_cellines_correlations_tissue.RData", topvar))
if(!file.exists(myfn)) { 
  cor.tissue <- NULL
  for(kk in 1:length(utissue)) {
    message("Correlations for tissue ", names(utissue[kk]))
    ## the top 1000 most variant genes in either CCLE or CGP
    iix.ccle <- !is.na(tissue.ccle) & tissue.ccle == names(utissue)[kk]
    iix.cgp <- !is.na(tissue.cgp) & tissue.cgp == names(utissue)[kk]
    varg <- unique(c(colnames(data.ccle)[order(apply(data.ccle[iix.ccle, , drop=FALSE], 2, var, na.rm=TRUE), decreasing=TRUE)[1:topvar]], colnames(data.cgp)[order(apply(data.cgp[iix.cgp, , drop=FALSE], 2, var, na.rm=TRUE), decreasing=TRUE)[1:topvar]]))
    rr <- cor(x=t(data.cgp[iix.cgp, varg, drop=FALSE]), y=t(data.ccle[iix.ccle, varg, drop=FALSE]), method="spearman", use="pairw")
    cor.tissue <- c(cor.tissue, list(rr))
  }
  names(cor.tissue) <- names(utissue)
  save(list="cor.tissue", compress=TRUE, file=myfn) 
} else { load(myfn) }

## boxplot
pdf(file.path(saveres, sprintf("ge_var%i_cellines_boxplot_ccle_cgp_tissue.pdf", topvar)), height=7, width=14)
ll <- c(list("all_tissues"=diag(cor.all)), lapply(cor.tissue, diag))
wt <- kruskal.test(x=ll)
myx <- which(lapply(ll, length) < 3)
pp <- NULL
for(j in myx) { for(jj in ll[[j]]) { pp <- rbind(pp, c(j, jj)) } }
ll[sapply(ll, length) < 3] <- NA
par(las=2, mar=c(12, 4, 4, 2) + 0.1, xaxt="n")
mp <- boxplot(ll, outline=FALSE, ylab="Rs", main="Correlations of gene expression profiles across tissue types\nCCLE vs. CGP", border="black", ylim=c(0.5, 1), col=c("lightgrey", coltissue))
points(x=pp[ , 1], y=pp[ , 2], col=coltissue[pp[ , 1]], pch=20)
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.02), pos=2, labels=mp$names, srt=45, xpd=TRUE, font=2, col=c("black", coltissuet))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
dev.off()

## boxplot
pdf(file.path(saveres, sprintf("ge_var%i_cellines_boxplot2_ccle_cgp_tissue_paper.pdf", topvar)), height=8, width=14)
ll <- c(list("all_tissues"=diag(cor.all)), lapply(cor.tissue, diag))
wt <- kruskal.test(x=ll)
ll[sapply(ll, length) < 3] <- NA
par(las=2, mar=c(12, 4, 4, 2) + 0.1, xaxt="n")
mp <- boxplot(ll, outline=FALSE, ylab="Rs", main="Correlations of gene expression profiles across tissue types\nCCLE vs. CGP", border="grey50", ylim=c(0, 1), col="white", range=0)
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
# xx$all_tissues <- rep(NA, length(xx$all_tissues))
points(x=jitter(unlist(xx), 1), y=unlist(ll), cex=0.75, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.02), pos=2, labels=mp$names, srt=45, xpd=TRUE, font=2, col=c("black", coltissuet))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
dev.off()

## test whether correlations are significantly higher in specific tissue types compared to all tissues combined
cor.diff.test <- NULL
ll <- c(list("all_tissues"=diag(cor.all)), lapply(cor.tissue, diag))
ddnn <- names(ll[[1]])
for(i in 2:length(ll)) {
  wt <- wilcox.test(x=ll[[1]], y=ll[[i]], alternative="two.sided", conf.int=TRUE)
  cor.diff.test <- rbind(cor.diff.test, c(wt$estimate, wt$p.value))
}
dimnames(cor.diff.test) <- list(names(ll)[-1], c("difference", "p"))
# cor.diff.test <- cbind(cor.diff.test, "lfdr"=fdrtool::fdrtool(cor.diff.test[ , 2], statistic="pvalue", verbose=FALSE)$lfdr)
cor.diff.test <- cbind(cor.diff.test, "fdr"=p.adjust(cor.diff.test[ , "p"], method="fdr"))
myx <- complete.cases(cor.diff.test) & (cor.diff.test[ , "p"] < 0.05 | cor.diff.test[ , "fdr"] < myfdr)
tt <- cor.diff.test[myx, , drop=FALSE]
tt <- tt[order(tt[ , "p"], decreasing=FALSE), , drop=FALSE]
message("Are correlations between gene expressions of identical cell lines\nin some tissue types significantly higher than in all tissues combined?")
print(tt)


################################################
## correlation of mutation data across cell lines in CGP and CCLE
################################################

message("Correlation of (missense) mutation data across cell lines (CGP vs CCLE)")

## how many cell lines with mutation data
message(sprintf("%i cel lines with mutation data both in CGP and CCLE", sum(apply(mutation.cgp, 1, function(x) { return(any(!is.na(x))) }) & apply(mutation.ccle, 1, function(x) { return(any(!is.na(x))) }))))

## find the common matches when multiple mutations weer detected
iix <- unique(c(grep("///", mutation.cgp), grep("///", mutation.ccle)))
ttt <- NULL
for(iii in iix) {
  ttt <- rbind(ttt, c(mutation.cgp[iii], mutation.ccle[iii]))
  x <- sort(unique(unlist(strsplit(mutation.cgp[iii], split="///"))))
  if(length(x) == 0) { x <- NA }
  y <- sort(unique(unlist(strsplit(mutation.ccle[iii], split="///"))))
  if(length(y) == 0) { y <- NA }
  tt <- intersect(x, y)  
  if(length(tt) != 0) {
    mutation.cgp[iii] <- mutation.ccle[iii] <- paste(tt, collapse="///")
  }
}

mutcor <- function(x, mutation1, mutation2) {
  require(epibasix)
  ff1 <- factor(mutation1[x[1], ], levels=c(0, 1))
  ff2 <- factor(mutation2[x[2], ], levels=c(0, 1))
  tt <- table(ff1, ff2)
  err <- try(rr <- epibasix::epiKappa(tt, k0=0)$kappa, silent=TRUE)
  if(class(err) == "try-error") { rr <- NA }
  return(rr)
}

## compare all mutations
myx <- (!is.na(mutation.cgp) & !is.na(mutation.ccle)) & (mutation.cgp != "wt" | mutation.ccle != "wt")
table(mutation.cgp[myx] == mutation.ccle[myx])
table(as.numeric(mutation.cgp != "wt")[myx] == as.numeric(mutation.ccle != "wt")[myx])
## compare mutations identified in CGP
myx <- (!is.na(mutation.cgp) & !is.na(mutation.ccle)) & (mutation.cgp != "wt")
table(mutation.cgp[myx] == mutation.ccle[myx])
table(as.numeric(mutation.cgp != "wt")[myx] == as.numeric(mutation.ccle != "wt")[myx])

## binarize the mutation matrix as in CGP and CCLE
mutation.cgp <- (mutation.cgp != "wt") * 1
mutation.ccle <- (mutation.ccle != "wt") * 1

myfn <- file.path(saveres, sprintf("mut_cellines_correlations.RData"))
if(!file.exists(myfn)) {
  ## compute Kappa between all pairs of cell lines
  tt <- apply(cbind(1:nrow(mutation.cgp), 1:nrow(mutation.ccle)), 1, list)
  tt <- lapply(tt, function(x) { return(x[[1]]) })
  combix <- c(tt, combn(1:nrow(mutation.cgp), 2, simplify=FALSE))
  mcres <- mclapply(combix, mutcor, mutation1=mutation.cgp, mutation2=mutation.ccle) 
  cor.all <- matrix(NA, nrow=nrow(mutation.cgp), ncol=nrow(mutation.ccle), dimnames=list(rownames(mutation.cgp), rownames(mutation.ccle)))
  iix <- do.call(rbind, combix)
  cor.all[iix] <- unlist(mcres)
  cor.all[iix[ , c(2, 1)]] <- unlist(mcres)
  save(list="cor.all", compress=TRUE, file=myfn)
} else { load(myfn) }
tt <- imputation::kNNImpute(x=cor.all, k=10, verbose=FALSE)$x
oo <- amap::hcluster(tt, link="complete")$order
pdf(file.path(saveres, sprintf("mut_cellines_boxplot_ccle_cgp_paper.pdf", topvar)), height=7, width=7)
par(xaxt="n", yaxt="n")
image(cor.all[oo, rev(oo)], col=blueyellow(255), xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="Correlation of mutation profiles across cell lines\nCCLE vs. CGP")
dev.off()
pdf(file.path(saveres, sprintf("mut_cellines_heatmap2_ccle_cgp.pdf", topvar)), height=8, width=8)
heatmap.3(x=cor.all[oo, oo], Colv=FALSE, Rowv=FALSE, scale="none", dendrogram="none", col=blueyellow(255), trace="none", density.info="none", keysize=1.5, xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="Gene expressions across cell lines\nCCLE vs. CGP", labRow=NA, labCol=NA, key=TRUE, KeyValueName="", keyTitle="Rs", breaks=seq(0, 1, length.out=256), mar=c(2, 2) + 0.1, na.color="lightgray")
dev.off()

## are the identical cell lines the mostly correlated?
cc <- sapply(1:nrow(cor.all), function(x, d) {
  return(sum(!is.na(d[x, -x]) & d[x, x] <= d[x, -x]))
}, d=cor.all)
names(cc) <- rownames(cor.all)
## cell lines which do not cocor.allelate most with their counterpart
message("These cell lines are not the most cocor.allelated with their counterpart (CCLE vs. CGP):\n", paste(names(cc)[!is.na(cc) & cc > 0], collapse="\t"))
## boxplot
pdf(file.path(saveres, sprintf("mut_cellines_boxplot_ccle_cgp_paper.pdf", topvar)), height=7, width=7)
ll <- list("Identical"=diag(cor.all), "Different"=cor.all[lower.tri(cor.all)])
wt <- wilcox.test(x=ll[[1]], y=ll[[2]])
par(xaxt="s", yaxt="s", mar=c(5, 4, 4, 2) + 0.1)
boxplot(ll, outline=FALSE, ylab="Kappa", xlab="Cell lines", main="Missense mutation profiles across cell lines\nCCLE vs. CGP", border="black", col="lightgrey", ylim=c(-0.20, 1))
legend("topright", legend=sprintf("Wilcoxon test p-value=%.1E", wt$p.value), bty="n")
dev.off()
message("")

## per tissue type
message("Correlation of gene expressions across cell lines per tissue type (CGP vs CCLE)")

myfn <- file.path(saveres, sprintf("mut_cellines_correlations_tissue.RData", topvar))
if(!file.exists(myfn)) { 
  cor.tissue <- NULL
  for(kk in 1:length(utissue)) {
    message("Correlations for tissue ", names(utissue[kk]))
    ## the top 1000 most variant genes in either CCLE or CGP
    iix.ccle <- !is.na(tissue.ccle) & tissue.ccle == names(utissue)[kk]
    iix.cgp <- !is.na(tissue.cgp) & tissue.cgp == names(utissue)[kk]
    rr <- matrix(NA, nrow=sum(iix.cgp), ncol=sum(iix.ccle), dimnames=list(rownames(mutation.cgp)[iix.cgp], rownames(mutation.ccle)[iix.ccle]))
    if(sum(iix.cgp) >= minsample && sum(iix.ccle) >= minsample) {
      ## compute Kappa between all pairs of cell lines
      tt <- apply(cbind(1:sum(iix.cgp), 1:sum(iix.ccle)), 1, list)
      tt <- lapply(tt, function(x) { return(x[[1]]) })
      combix <- c(tt, combn(1:sum(iix.cgp), 2, simplify=FALSE))
      mcres <- mclapply(combix, mutcor, mutation1=mutation.cgp[iix.cgp, , drop=FALSE], mutation2=mutation.ccle[iix.ccle, , drop=FALSE])
      iix <- do.call(rbind, combix)
      rr[iix] <- unlist(mcres)
      rr[iix[ , c(2, 1)]] <- unlist(mcres)      
    }
    cor.tissue <- c(cor.tissue, list(rr))
  }
  names(cor.tissue) <- names(utissue)
  save(list="cor.tissue", compress=TRUE, file=myfn) 
} else { load(myfn) }

## boxplot
pdf(file.path(saveres, sprintf("mut_cellines_boxplot_ccle_cgp_tissue.pdf", topvar)), height=7, width=14)
ll <- c(list("all_tissues"=diag(cor.all)), lapply(cor.tissue, diag))
wt <- kruskal.test(x=ll[sapply(ll, function(x) { return(sum(!is.na(x)) >= 3) })])
myx <- which(lapply(ll, length) < 3)
pp <- NULL
for(j in myx) { for(jj in ll[[j]]) { pp <- rbind(pp, c(j, jj)) } }
ll[sapply(ll, length) < 3] <- NA
par(las=2, mar=c(12, 4, 6, 2) + 0.1, xaxt="n")
mp <- boxplot(ll, outline=FALSE, ylab="Rs", main="Agreement of missense mutation profiles across tissue types\nCCLE vs. CGP", border="black", ylim=c(-0.1, 1.1), col=c("lightgrey", coltissue))
points(x=pp[ , 1], y=pp[ , 2], col=coltissue[pp[ , 1]], pch=20)
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.02), pos=2, labels=mp$names, srt=45, xpd=TRUE, font=2, col=c("black", coltissuet))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
dev.off()

## boxplot
pdf(file.path(saveres, sprintf("mut_cellines_boxplot2_ccle_cgp_tissue_paper.pdf", topvar)), height=8, width=14)
ll <- c(list("all_tissues"=diag(cor.all)), lapply(cor.tissue, diag))
wt <- kruskal.test(x=ll[sapply(ll, function(x) { return(sum(!is.na(x)) >= 3) })])
ll[sapply(ll, length) < 3] <- NA
par(las=2, mar=c(12, 4, 4, 2) + 0.1, xaxt="n")
mp <- boxplot(ll, outline=FALSE, ylab="Rs", main="Agreement of missense mutation profiles across tissue types\nCCLE vs. CGP", border="grey50", ylim=c(-0.1, 1.1), col="white", range=0)
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
# xx$all_tissues <- rep(NA, length(xx$all_tissues))
points(x=jitter(unlist(xx), 1), y=unlist(ll), cex=0.75, pch=20, col=unlist(ccol))
axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.02), pos=2, labels=mp$names, srt=45, xpd=TRUE, font=2, col=c("black", coltissuet))
legend("topright", legend=sprintf("Kruskal-Wallis test p-value=%.1E", wt$p.value), bty="n")
dev.off()

## test whether correlations are significantly higher in specific tissue types compared to all tissues combined
cor.diff.test <- NULL
ll <- c(list("all_tissues"=diag(cor.all)), lapply(cor.tissue, diag))
ddnn <- names(ll[[1]])
for(i in 2:length(ll)) {
  if(sum(!is.na(ll[[i]])) >= 3) {
    wt <- wilcox.test(x=ll[[1]], y=ll[[i]], alternative="two.sided", conf.int=TRUE)
  } else { wt <- list("estimate"=NA, "p.value"=NA) }
  cor.diff.test <- rbind(cor.diff.test, c(wt$estimate, wt$p.value))
}
dimnames(cor.diff.test) <- list(names(ll)[-1], c("difference", "p"))
# cor.diff.test <- cbind(cor.diff.test, "lfdr"=fdrtool::fdrtool(cor.diff.test[ , 2], statistic="pvalue", verbose=FALSE)$lfdr)
cor.diff.test <- cbind(cor.diff.test, "fdr"=p.adjust(cor.diff.test[ , "p"], method="fdr"))
myx <- complete.cases(cor.diff.test) & (cor.diff.test[ , "p"] < 0.05 | cor.diff.test[ , "fdr"] < myfdr)
tt <- cor.diff.test[myx, , drop=FALSE]
tt <- tt[order(tt[ , "p"], decreasing=FALSE), , drop=FALSE]
message("Are correlations between gene expressions of identical cell lines\nin some tissue types significantly higher than in all tissues combined?")
print(tt)


#################################################
## Correlation of ic50s across cell lines
#################################################

message("Correlation of IC50s across cell lines (CGP vs CCLE)")

myfn <- file.path(saveres, "ic50_cellines_correlations.RData")
if(!file.exists(myfn)) {
  rr <- cor(t(ic50.cgp), t(ic50.ccle), method="spearman", use="pairw")
  save(list="rr", compress=TRUE, file=myfn)
} else { load(myfn) }
## use the same order than the sone computed from gene expression data
pdf(file.path(saveres, "ic50_cellines_heatmap_ccle_cgp.pdf"), height=7, width=7)
par(xaxt="n", yaxt="n")
image(rr[oo, rev(oo)], col=blueyellow(255), xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="IC50s across cell lines\nCCLE vs. CGP")
dev.off()
pdf(file.path(saveres, "ic50_cellines_heatmap2_ccle_cgp_paper.pdf"), height=8, width=8)
heatmap.3(x=rr[oo, oo], Colv=FALSE, Rowv=FALSE, scale="none", dendrogram="none", col=blueyellow(255), trace="none", density.info="none", keysize=1.5, xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="IC50s across cell lines\nCCLE vs. CGP", labRow=NA, labCol=NA, key=TRUE, KeyValueName="", keyTitle="Rs", breaks=seq(0, 1, length.out=256), mar=c(2, 2) + 0.1, na.color="lightgray")
dev.off()
## remove most of the missing values
naix <- apply(rr, 1, function(x, y) { return((sum(is.na(x)) / length(x)) >= y) }, y=0.5) | apply(rr, 2, function(x, y) { return((sum(is.na(x)) / length(x)) >= y) }, y=0.5)
rr2 <- rr[!naix, !naix]
oo <- amap::hcluster(rr2, link="complete")$order
## use the same order than the sone computed from gene expression data
pdf(file.path(saveres, "ic50_cellines_heatmap_lessna_ccle_cgp.pdf"), height=7, width=7)
par(xaxt="n", yaxt="n")
image(rr2[oo, oo], col=blueyellow(255), xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="IC50s across cell lines\nCCLE vs. CGP")
dev.off()
pdf(file.path(saveres, "ic50_cellines_heatmap2_lessna_ccle_cgp.pdf"), height=8, width=8)
heatmap.3(x=rr2[oo, oo], Colv=FALSE, Rowv=FALSE, scale="none", dendrogram="none", col=blueyellow(255), trace="none", density.info="none", keysize=1.5, xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="IC50s across cell lines\nCCLE vs. CGP", labRow=NA, labCol=NA, key=TRUE, KeyValueName="", keyTitle="Rs", breaks=seq(0, 1, length.out=256), mar=c(2, 2) + 0.1, na.color="lightgray")
dev.off()

## are the identical cell lines the mostly correlated?
cc <- sapply(1:nrow(rr), function(x, d) {
  return(sum(!is.na(d[x, -x]) & d[x, x] <= d[x, -x]))
}, d=rr)
names(cc) <- rownames(rr)
## cell lines which do not correlate most with their counterpart
message("These cell lines are not the most correlated with their counterpart (CCLE vs. CGP):\n", paste(names(cc)[!is.na(cc) & cc > 0], collapse="\t"))
## boxplot
pdf(file.path(saveres, "ic50_cellines_boxplot_ccle_cgp.pdf"), height=7, width=7)
ll <- list("Identical"=diag(rr), "Different"=rr[lower.tri(rr)])
wt <- wilcox.test(x=ll[[1]], y=ll[[2]])
par(xaxt="s", yaxt="s", mar=c(5, 4, 4, 2) + 0.1)
boxplot(ll, outline=FALSE, ylab="Rs", xlab="Cell lines", main="IC50s across cell lines\nCCLE vs. CGP", border="black", col="lightgrey", ylim=c(-0.10, 1))
legend("bottomleft", legend=sprintf("Wilcoxon test p-value=%.1E", wt$p.value), bty="n")
dev.off()
message("")

########################
## Correlation of IC50 calls across cell lines

message("Correlation of IC50 calls across cell lines (CGP vs CCLE)")

myfn <- file.path(saveres, "ic50_call_cellines_correlations.RData")
if(!file.exists(myfn)) {
  iix <- t(combn(nrow(ic50.call.ccle), 2, simplify=TRUE))
  iix <- rbind(iix, iix[ ,2:1], cbind(1:nrow(ic50.call.ccle), 1:nrow(ic50.call.ccle)))
  splitix <- parallel::splitIndices(nx=nrow(iix), ncl=nbcore)
  mcres <- parallel::mclapply(splitix, function(x, y, d1, d2) {
    res <- apply(y[x, , drop=FALSE], 1, function(xx, dd1, dd2) { 
      ff1 <- factor(x=dd1[xx[1], ], levels=c("resistant", "intermediate", "sensitive"))
      ff2 <- factor(x=dd2[xx[2], ], levels=c("resistant", "intermediate", "sensitive"))
      tt <- table("CCLE"=ff1, "CGP"=ff2)
      ## check if all elements are on the diagonal
      ttt <- tt
      diag(ttt) <- 0
      if(sum(ttt) == 0) {
        rr <- 1
      } else {
        err <- try(rr <- epibasix::epiKappa(tt, k0=0)$kappa, silent=TRUE)
        if(class(err) == "try-error") { rr <- NA }
      }
      return(rr)
    }, dd1=d1, dd2=d2)
    return(res)
  }, y=iix, d1=ic50.call.ccle, d2=ic50.call.cgp)
  mcres <- do.call(c, mcres)
  rr <- matrix(NA, nrow=nrow(ic50.call.cgp), ncol=nrow(ic50.call.ccle), dimnames=list(rownames(ic50.call.cgp), rownames(ic50.call.ccle)))
  rr[iix] <- mcres
  rr <- t(rr)
  save(list="rr", compress=TRUE, file=myfn)
} else { load(myfn) }

## use the same order than the sone computed from gene expression data
pdf(file.path(saveres, "ic50_call_cellines_heatmap_ccle_cgp.pdf"), height=7, width=7)
par(xaxt="n", yaxt="n")
image(rr[oo, rev(oo)], col=blueyellow(255), xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="IC50 calls across cell lines\nCCLE vs. CGP")
dev.off()
pdf(file.path(saveres, "ic50_call_cellines_heatmap2_ccle_cgp.pdf"), height=8, width=8)
heatmap.3(x=rr[oo, oo], Colv=FALSE, Rowv=FALSE, scale="none", dendrogram="none", col=blueyellow(255), trace="none", density.info="none", keysize=1.5, xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="IC50 calls across cell lines\nCCLE vs. CGP", labRow=NA, labCol=NA, key=TRUE, KeyValueName="", keyTitle="Kappa", breaks=seq(0, 1, length.out=256), mar=c(2, 2) + 0.1, na.color="lightgray")
dev.off()
## remove most of the missing values
naix <- apply(rr, 1, function(x, y) { return((sum(is.na(x)) / length(x)) >= y) }, y=0.5) | apply(rr, 2, function(x, y) { return((sum(is.na(x)) / length(x)) >= y) }, y=0.5)
rr2 <- rr[!naix, !naix]
oo <- amap::hcluster(rr2, link="complete")$order
## use the same order than the sone computed from gene expression data
pdf(file.path(saveres, "ic50_call_cellines_heatmap_lessna_ccle_cgp.pdf"), height=7, width=7)
par(xaxt="n", yaxt="n")
image(rr2[oo, oo], col=blueyellow(255), xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="IC50 calls across cell lines\nCCLE vs. CGP")
dev.off()
pdf(file.path(saveres, "ic50_call_cellines_heatmap2_lessna_ccle_cgp.pdf"), height=8, width=8)
heatmap.3(x=rr2[oo, oo], Colv=FALSE, Rowv=FALSE, scale="none", dendrogram="none", col=blueyellow(255), trace="none", density.info="none", keysize=1.5, xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="IC50 calls across cell lines\nCCLE vs. CGP", labRow=NA, labCol=NA, key=TRUE, KeyValueName="", keyTitle="Rs", breaks=seq(0, 1, length.out=256), mar=c(2, 2) + 0.1, na.color="lightgray")
dev.off()

## are the identical cell lines the mostly correlated?
cc <- sapply(1:nrow(rr), function(x, d) {
  return(sum(!is.na(d[x, -x]) & d[x, x] <= d[x, -x]))
}, d=rr)
names(cc) <- rownames(rr)
## cell lines which do not correlate most with their counterpart
message("These cell lines are not the most correlated with their counterpart (CCLE vs. CGP):\n", paste(names(cc)[!is.na(cc) & cc > 0], collapse="\t"))
## boxplot
pdf(file.path(saveres, "ic50_call_cellines_boxplot_ccle_cgp.pdf"), height=7, width=7)
ll <- list("Identical"=diag(rr), "Different"=rr[lower.tri(rr)])
wt <- wilcox.test(x=ll[[1]], y=ll[[2]])
par(xaxt="s", yaxt="s", mar=c(5, 4, 4, 2) + 0.1)
boxplot(ll, outline=FALSE, ylab="Kappa", xlab="Cell lines", main="IC50 calls across cell lines\nCCLE vs. CGP", border="black", col="lightgrey", ylim=c(-0.10, 1))
legend("bottomleft", legend=sprintf("Wilcoxon test p-value=%.1E", wt$p.value), bty="n")
dev.off()
message("")

########################
## Correlation of AUCs across cell lines

message("Correlation of AUCs across cell lines (CGP vs CCLE)")

myfn <- file.path(saveres, "auc_cellines_correlations.RData")
if(!file.exists(myfn)) {
  rr <- cor(t(auc.cgp), t(auc.ccle), method="spearman", use="pairw")
  save(list="rr", compress=TRUE, file=myfn)
} else { load(myfn) }
## use the same order than the sone computed from gene expression dataf
pdf(file.path(saveres, "auc_cellines_heatmap_ccle_cgp.pdf"), height=7, width=7)
par(xaxt="n", yaxt="n")
image(rr[oo, rev(oo)], col=blueyellow(255), xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="AUCs across cell lines\nCCLE vs. CGP")
dev.off()
pdf(file.path(saveres, "auc_cellines_heatmap2_ccle_cgp.pdf"), height=8, width=8)
heatmap.3(x=rr[oo, oo], Colv=FALSE, Rowv=FALSE, scale="none", dendrogram="none", col=blueyellow(255), trace="none", density.info="none", keysize=1.5, xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="AUCs across cell lines\nCCLE vs. CGP", labRow=NA, labCol=NA, key=TRUE, KeyValueName="", keyTitle="Rs", breaks=seq(0, 1, length.out=256), mar=c(2, 2) + 0.1, na.color="lightgray")
dev.off()
## remove most of the missing values
naix <- apply(rr, 1, function(x, y) { return((sum(is.na(x)) / length(x)) >= y) }, y=0.5) | apply(rr, 2, function(x, y) { return((sum(is.na(x)) / length(x)) >= y) }, y=0.5)
rr2 <- rr[!naix, !naix]
oo <- amap::hcluster(rr2, link="complete")$order
## use the same order than the sone computed from gene expression data
pdf(file.path(saveres, "auc_cellines_heatmap_lessna_ccle_cgp.pdf"), height=7, width=7)
par(xaxt="n", yaxt="n")
image(rr2[oo, oo], col=blueyellow(255), xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="AUCs across cell lines\nCCLE vs. CGP")
dev.off()
pdf(file.path(saveres, "auc_cellines_heatmap2_lessna_ccle_cgp.pdf"), height=8, width=8)
heatmap.3(x=rr2[oo, oo], Colv=FALSE, Rowv=FALSE, scale="none", dendrogram="none", col=blueyellow(255), trace="none", density.info="none", keysize=1.5, xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="AUCs across cell lines\nCCLE vs. CGP", labRow=NA, labCol=NA, key=TRUE, KeyValueName="", keyTitle="Rs", breaks=seq(0, 1, length.out=256), mar=c(2, 2) + 0.1, na.color="lightgray")
dev.off()

## are the identical cell lines the mostly correlated?
cc <- sapply(1:nrow(rr), function(x, d) {
  return(sum(!is.na(d[x, -x]) & d[x, x] <= d[x, -x]))
}, d=rr)
names(cc) <- rownames(rr)
## cell lines which do not correlate most with their counterpart
message("These cell lines are not the most correlated with their counterpart (CCLE vs. CGP):\n", paste(names(cc)[!is.na(cc) & cc > 0], collapse="\t"))
## boxplot
pdf(file.path(saveres, "auc_cellines_boxplot_ccle_cgp.pdf"), height=7, width=7)
ll <- list("Identical"=diag(rr), "Different"=rr[lower.tri(rr)])
wt <- wilcox.test(x=ll[[1]], y=ll[[2]])
par(xaxt="s", yaxt="s", mar=c(5, 4, 4, 2) + 0.1)
boxplot(ll, outline=FALSE, ylab="Rs", xlab="Cell lines", main="AUCs across cell lines\nCCLE vs. CGP", border="black", col="lightgrey", ylim=c(-0.10, 1))
legend("bottomleft", legend=sprintf("Wilcoxon test p-value=%.1E", wt$p.value), bty="n")
dev.off()
message("")

########################
## Correlation of AUC calls across cell lines

message("Correlation of AUC calls across cell lines (CGP vs CCLE)")

myfn <- file.path(saveres, "auc_call_cellines_correlations.RData")
if(!file.exists(myfn)) {
  iix <- t(combn(nrow(auc.call.ccle), 2, simplify=TRUE))
  iix <- rbind(iix, iix[ ,2:1], cbind(1:nrow(auc.call.ccle), 1:nrow(auc.call.ccle)))
  splitix <- parallel::splitIndices(nx=nrow(iix), ncl=nbcore)
  mcres <- parallel::mclapply(splitix, function(x, y, d1, d2) {
    res <- apply(y[x, , drop=FALSE], 1, function(xx, dd1, dd2) { 
      ff1 <- factor(x=dd1[xx[1], ], levels=c("resistant", "intermediate", "sensitive"))
      ff2 <- factor(x=dd2[xx[2], ], levels=c("resistant", "intermediate", "sensitive"))
      tt <- table("CCLE"=ff1, "CGP"=ff2)
      ## check if all elements are on the diagonal
      ttt <- tt
      diag(ttt) <- 0
      if(sum(ttt) == 0) {
        rr <- 1
      } else {
        err <- try(rr <- epibasix::epiKappa(tt, k0=0)$kappa, silent=TRUE)
        if(class(err) == "try-error") { rr <- NA }
      }
      return(rr)
    }, dd1=d1, dd2=d2)
    return(res)
  }, y=iix, d1=auc.call.ccle, d2=auc.call.cgp)
  mcres <- do.call(c, mcres)
  rr <- matrix(NA, nrow=nrow(auc.call.cgp), ncol=nrow(auc.call.ccle), dimnames=list(rownames(auc.call.cgp), rownames(auc.call.ccle)))
  rr[iix] <- mcres
  rr <- t(rr)
  save(list="rr", compress=TRUE, file=myfn)
} else { load(myfn) }

## use the same order than the sone computed from gene expression data
pdf(file.path(saveres, "auc_call_cellines_heatmap_ccle_cgp.pdf"), height=7, width=7)
par(xaxt="n", yaxt="n")
image(rr[oo, rev(oo)], col=blueyellow(255), xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="AUC calls across cell lines\nCCLE vs. CGP")
dev.off()
pdf(file.path(saveres, "auc_call_cellines_heatmap2_ccle_cgp.pdf"), height=8, width=8)
heatmap.3(x=rr[oo, oo], Colv=FALSE, Rowv=FALSE, scale="none", dendrogram="none", col=blueyellow(255), trace="none", density.info="none", keysize=1.5, xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="AUC calls across cell lines\nCCLE vs. CGP", labRow=NA, labCol=NA, key=TRUE, KeyValueName="", keyTitle="Kappa", breaks=seq(0, 1, length.out=256), mar=c(2, 2) + 0.1, na.color="lightgray")
dev.off()
## use the same order than the sone computed from gene expression data
pdf(file.path(saveres, "auc_call_cellines_heatmap_lessna_ccle_cgp.pdf"), height=7, width=7)
par(xaxt="n", yaxt="n")
image(rr2[oo, oo], col=blueyellow(255), xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="AUC calls across cell lines\nCCLE vs. CGP")
dev.off()
pdf(file.path(saveres, "auc_call_cellines_heatmap2_lessna_ccle_cgp.pdf"), height=8, width=8)
heatmap.3(x=rr2[oo, oo], Colv=FALSE, Rowv=FALSE, scale="none", dendrogram="none", col=blueyellow(255), trace="none", density.info="none", keysize=1.5, xlab="Cell lines (CGP)", ylab="Cell lines (CCLE)", main="AUC calls across cell lines\nCCLE vs. CGP", labRow=NA, labCol=NA, key=TRUE, KeyValueName="", keyTitle="Rs", breaks=seq(0, 1, length.out=256), mar=c(2, 2) + 0.1, na.color="lightgray")
dev.off()

## are the identical cell lines the mostly correlated?
cc <- sapply(1:nrow(rr), function(x, d) {
  return(sum(!is.na(d[x, -x]) & d[x, x] <= d[x, -x]))
}, d=rr)
names(cc) <- rownames(rr)
## cell lines which do not correlate most with their counterpart
message("These cell lines are not the most correlated with their counterpart (CCLE vs. CGP):\n", paste(names(cc)[!is.na(cc) & cc > 0], collapse="\t"))
## boxplot
pdf(file.path(saveres, "auc_call_cellines_boxplot_ccle_cgp.pdf"), height=7, width=7)
ll <- list("Identical"=diag(rr), "Different"=rr[lower.tri(rr)])
wt <- wilcox.test(x=ll[[1]], y=ll[[2]])
par(xaxt="s", yaxt="s", mar=c(5, 4, 4, 2) + 0.1)
boxplot(ll, outline=FALSE, ylab="Kappa", xlab="Cell lines", main="AUC calls across cell lines\nCCLE vs. CGP", border="black", col="lightgrey", ylim=c(-0.10, 1))
legend("bottomleft", legend=sprintf("Wilcoxon test p-value=%.1E", wt$p.value), bty="n")
dev.off()
message("")
