
## NB, first change working directory to the directory that contains this script!

## NB to run this file, you must download these three RData files below from XXXXXXXX
## and be sure to store them in the same directory as this script, and set this 
## directory as the current working directory!
## Load the data files that were prepared from Haibe-Kains et al original code...
load("ge_var1000_cellines_correlations_replicates_cgp.RData")
load("ge_var1000_cellines_correlations.RData")
load("CDRUG_cgp_ccle.RData")


############################################################################
######## Code below copied from "CDRUG_pipeline.R" (lines 21-68) ###########
############################################################################

## define functions required for the analysis pipeline
source(file.path("CDRUG_foo.R"))

########################
## global parameters

## set random seed to ensuer reproducibility of the resuls
set.seed(54321)


## list of characters to be removed from row and column names
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

## directory where all the analysis results will be stored

saveres <- "." ## Geeleher et al.: In order to make this run, this line has been added and the two lines below have been commented out.
# saveres <- "saveres"
# if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE, recursive=TRUE) }

## minimum number of samples to compute correlation
minsample <- 10

## max fdr, threshold used to identify genes with significant concordance index
myfdr <- 0.20

## method to estimate the association gene-drug, controlled for tissue type
# genedrugm <- c("lm", "cindex")
genedrugm <- "lm"

## number of most variant genes to consider for correlation
topvar <- 1000


############################################################################
######## End Code copied from "CDRUG_pipeline.R" (lines 17-68) #############
############################################################################


######################################################################
############## Begin Geeleher et. al Code ############################
######################################################################

require(survcomp) || stop("Library survcomp is not available!")
require(ggplot2) || stop("Library ggplot2 is not available!")


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


#####################################################################
######## Code below copied from "CDRUG_analysis.R" ##################
#####################################################################

###############################################################################
###############################################################################
## load data
###############################################################################
###############################################################################

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
# xtable::print.xtable(xtable::xtable(mm), include.rownames=FALSE, floating=FALSE, file=file.path(saveres, "tissue_type_cgp_ccle_paper.tex"), append=FALSE)

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
## AUC

message("Correlation AUC measures")
# ## consistency between AUC with spearman correlation
# pdf(file.path(saveres, "auc_spearman_ccle_cgp_pres.pdf"), height=9, width=16)
# par(mfrow=c(3, 5), cex=0.8, las=1)
for(i in 1:ncol(auc.ccle)) {
  xxlim <- yylim <- round(range(c(auc.cgp[, i], auc.ccle[, i]), na.rm=TRUE) * 10) / 10
  nnn <- sum(complete.cases(auc.ccle[, i], auc.cgp[, i]))
  if(nnn >= minsample) {
    cc <- cor.test(auc.ccle[, i], auc.cgp[, i], method="spearman", use="complete.obs", alternative="greater")
  } else {
    cc <- list("estimate"=NA, "p.value"=NA)
  }
#   par(mar=c(4, 4, 3, 1) + 0.1)
#   myScatterPlot(x=auc.cgp[, i], y=auc.ccle[, i], xlab=ifelse(i > 10, "AUC (CGP)", ""), ylab=ifelse((i %% 5) == 1, "AUC (CCLE)", ""), main=gsub("drugid_", "", colnames(auc.ccle)[i]), xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75)
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("R=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, sum(complete.cases(auc.ccle[ ,i], auc.cgp[ ,i]))), text.font=2)
  correlations[["auc"]][i, "drug.sensitivity"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
  correlations.stats[["auc"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
}
# dev.off()


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
ggplot(data=dat, aes(x=x, y=y)) + theme_bw() + geom_point(aes(color=varCpg), size=I(3)) + geom_text(aes(label=Drug), vjust=-.5, hjust=-.24, size=2.5, angle=15) + xlab("Correlation of AUC between CCLE and CGP") + ylab("Variance of AUC in CCLE") + scale_color_continuous(low="steelblue4",high="tomato2", name=expression(sigma^2 ~ "CPG")) + theme(legend.position=c(.1,.82))
dev.off()

## The Spearman correlation tests below show the obvious relationship between lack of drug response (i.e. high median AUC) and lack of variability.
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

###################################################################
######## End code copied from "CDRUG_analysis.R" ##################
###################################################################












