runOverallMultivarJobTMT <- function(...) {

##############################################
# main function for overall analysis
##############################################

library(heatmap3)
library(lattice)
library(openxlsx)
library(XLConnect) # only for printStyle
library(VennDiagram)
library(scatterplot3d) 
library(TMTPrePro)


############
# parameters
############

# zipfname = "files.zip"
# designfname = "designRuns.xlsx"
FCCutoff = 1.5
ZScoreCutoff = 2
CountsCutoff = 0

args <- list(...)

for(i in 1:length(args)) {
flag <- substring(args[[i]], 0, 2)
value <- substring(args[[i]], 3, nchar(args[[i]]))


if(flag=='-f') zipfname <- value;
if(flag=='-d') designfname <- value;
if(flag=='-r') FCCutoff <- as.numeric(value);
if(flag=='-z') ZScoreCutoff <- as.numeric(value);
if(flag=='-c') CountsCutoff <- as.numeric(value);

} 

	

# unzip files and check zip
files <- unzip(zipfname, junkpaths=TRUE)
if (length(files) < 1) stop("The zip file is empty")

	
# Error checking, have two tabs etc?	
designAll = try( readWorkbook(designfname, 1) )
if (inherits(designAll, "try-error")) stop("Could not extract the first tab from the design file.")

references = readWorkbook(designfname, 2)
if (inherits(references, "try-error")) stop("Could not extract the second tab from the design file.")

files = references[,1]
refLabels = references[,2]


NN = nrow(references)
if (NN < 1) stop("The second tab of the design has no files.")

if (ncol(designAll) < (2*NN + 1)) stop("The design file does not have the right format.")

design.list = list()

for (i in 1:NN) {
	design.list[[i]] = designAll[, c(1, 2*i, 2*i+1)]
}
			
# match reference groups using design
refGroups = try(lapply(1:NN, FUN=function(nn){design=design.list[[nn]];design[match(refLabels[nn], design$Label),"Group"]}))
if (inherits(refGroups, "try-error")) stop("Could not parse design to extract the reference; check Group columns exist.");

if (length(unique(unlist(refGroups))) > 1) cat("Warning: the references belong to different groups; data interpretation may not make sense.");

file.list = try(lapply(files, FUN=read.delim))
if (inherits(file.list,  "try-error")) stop("Could not read TMT files or files do not exist; check correct names in second tab of design and files present in zip.")

ratios.list = list()


for (fidx in 1:length(file.list)) { 

dat1 = file.list[[fidx]]
design = design.list[[fidx]]
referenceLabel = refLabels[[fidx]]
referenceGroup = refGroups[[fidx]]

counts.idx = grep("Count$", names(dat1))
var.idx = setdiff(grep("Variability", names(dat1)), grep("^New", names(dat1)))

Counts1 = dat1[,counts.idx][,1]
Var1 = dat1[,var.idx]

# extract ratios only
# exclude counts, variability, "New" or AA to end, if exists
endcols = NULL
newcol = grep("^(New|X..AA)", names(dat1))
if (length(newcol) > 1) {
	endcols = newcol[1]:ncol(dat1)
} 
ratios1 = as.matrix(dat1[,-c(1:(-1+grep("^X12", names(dat1))[1]), counts.idx, var.idx, endcols)])


# extract top and bottom ratios
top1 = gsub("X(.*)\\.(.*)", "\\1", colnames(ratios1))
bot1 = gsub("X(.*)\\.(.*)", "\\2", colnames(ratios1))

GroupTop = as.factor(design$Group[match(top1, design$Label)])
GroupBot = as.factor(design$Group[match(bot1, design$Label)])

# select only ratios TO the reference; 

ref.idx = ( bot1 == referenceLabel)  #              

ratiosRef = ratios1[,ref.idx]
varRef = Var1[, ref.idx]
topRef = top1[ref.idx]

# keep accession only, or all?
ratiosRef = data.frame(dat1[,"Accession"], ratiosRef)
varRef = data.frame(dat1[,"Accession"],varRef)
GroupRef = GroupTop[ref.idx]

# Select same-same ratios, Reference to Reference (e.g. Ctr/Ctr, for FDR)
ratiosRefRef = (ratiosRef[,-1][ ,GroupRef == referenceGroup, drop=FALSE])
varRefRef = (varRef[,-1][ ,GroupRef == referenceGroup, drop=FALSE])

NSame = ncol(ratiosRefRef)

if (NSame > 0) {
for (jj in 1:NSame) {

# for each same-same ratio (could be several), generate images and FDR
Var = varRefRef[,jj,drop=T]
Rat = ratiosRefRef[,jj]
zscore = 100*log(Rat)/Var

imgName = paste("File", fidx, "Ref", referenceLabel, "Ratio", jj, ".png", sep="")


FDRRat = sum( !is.na(Rat) & !is.na(Counts1) & (abs(log(Rat)) > log(FCCutoff))) /length(Rat)
FDRRatCounts =sum( !is.na(Rat) & !is.na(Counts1) & (abs(log(Rat)) > log(FCCutoff)) & ( Counts1 > 1)) /length(Rat)
FDRRatZ = sum( !is.na(Rat) & !is.na(zscore) & (abs(log(Rat)) > log(FCCutoff)) & (log(abs(zscore)) > log(ZScoreCutoff))) /length(Rat)

png(imgName, width =3000, height=2000, res=300)

layout(matrix(1:2, nrow=1))
plot(log(Rat), log(Counts1), pch = 20, main=paste("Ratios", names(ratiosRefRef)[jj], "and counts"), 
	xlab=paste("FDR Ratio cutoff", FCCutoff, "=", format(FDRRat, digits=3)),
	sub=paste("FDR Ratio cutoff", FCCutoff, "Counts >1", "=", format(FDRRatCounts, digits=3)))
abline(v=log(FCCutoff), col="blue")
abline(v=-log(FCCutoff), col="blue")

plot(log(Rat), log(abs(zscore)), pch=20, main=paste("Ratios", names(ratiosRefRef)[jj], "and Z-score"),
	xlab=paste("FDR Ratio cutoff", FCCutoff, "=", format(FDRRat, digits=3)),
	sub=paste("FDR Ratio cutoff", FCCutoff, "ZScore >", ZScoreCutoff, "=", format(FDRRatZ, digits=3)))
abline(v=log(FCCutoff), col="blue")
abline(v=-log(FCCutoff), col="blue")
abline(h=log(ZScoreCutoff), col="red")

dev.off()


} # end for all same-same ratios

}

# add variability and counts

ratios.list[[fidx]] = list(Ratios = ratios1, ratiosRef = ratiosRef, Group=GroupRef, GroupTop=GroupTop, GroupBot=GroupBot, Var=varRef, Counts=data.frame(dat1[,"Accession"],Counts1))

} # end for all files


# merge ratios to selected reference from all files

mg.rat = ratios.list[[1]]$ratiosRef
mg.var = ratios.list[[1]]$Var
mg.Counts = ratios.list[[1]]$Counts
Group = as.factor(as.vector(ratios.list[[1]]$Group))
colnames(mg.rat) = paste("R1", colnames(mg.rat))


if (length(ratios.list) > 1) {

for (xx in 2:length(ratios.list)) {

nextrun = ratios.list[[xx]]$ratiosRef
colnames(nextrun) = paste("R", xx, colnames(nextrun), sep="")

nextVar = ratios.list[[xx]]$Var
colnames(nextVar) = paste("R", xx, colnames(nextVar), sep="")

nextCount = ratios.list[[xx]]$Counts
colnames(nextCount) = paste("R", xx, colnames(nextCount), sep="")

  mg = merge(mg.rat, nextrun , by.x=1, by.y=1, all=TRUE)
  vvar = merge(mg.var, nextVar , by.x=1, by.y=1, all=TRUE)
  counts = merge(mg.Counts, nextCount, by.x=1, by.y=1, all=TRUE)
  Gp = as.factor(c(as.vector(Group), as.vector(ratios.list[[xx]]$Group) ) )
  
  mg.rat = mg
  mg.var = vvar
  mg.Counts = counts
  Group = Gp
  
}
}


# check ordering same and not modified by merging
if (sum(mg.rat[,1] != mg.var[,1]) > 0) stop("Error in merging runs: ratios and var mismatch");
if (sum(mg.rat[,1] != mg.Counts[,1]) > 0) stop("Error in merging runs: ratios and counts mismatch");


# remove identifier column from ratios
rownames(mg.rat) = mg.rat[,1]
mg.rat = as.matrix(mg.rat[,-1])

full.ratios.only = data.frame(mg.Counts, mg.rat[,order(Group)])
full.res = data.frame(mg.Counts, mg.var[, c(1,1+order(Group))])

####################################################
# Some overall data look metrics
# -- correlation
# -- boxplots and density plots
# -- within group correlation for each level of the group
# -- PCA
# -- heatmap
# -- Anova, followed by clustering of DE proteins
# -- Diff exp proteins _to the reference_, e.g. Mod/Control, ..., etc, via 1-sample t-test,
# -- Also combine z-scores
# -- Venn diagrams of overlap, and also barplots
####################################################


# good correlation across runs based on ratios to indicated reference
png("Correlation heatmap.png", 2000, 2000,res=300)
heatmap3(cor(log(na.omit(mg.rat)), use="pairwise.complete.obs"), distfun=cordist,
     col=colorRampPalette(c("green","black", "red"))(1024),
	ColSideColors=rainbow(nlevels(Group))[Group], margins=c(10,10))
legend("topright", fill=rainbow(nlevels(Group)), legend=levels(Group))
dev.off()

png("BoxplotDensity.png", width=3500, height=1700,res=300)
layout(matrix(1:2, nrow=1))

plotDensity(t(log(na.omit(mg.rat))), Group)

# boxplots and density plots
boxplot(log(mg.rat[, order(Group)]), las=2, col=rainbow(nlevels(Group))[Group[order(Group)]], ylim=c(-1.5,1.5), cex.axis=0.6)
legend("topright", fill=rainbow(nlevels(Group)), legend=levels(Group))


dev.off()


# Correlations for all levels of the group
for (lev in levels(Group)) {

cmat = cor(mg.rat[,Group == lev, drop=FALSE],method="spearman", use="pairwise.complete.obs")
NN = nrow(cmat)
main = paste(lev, "; min cor:", format(min( cmat[upper.tri(matrix(NA, NN, NN))] ), digits=3), ", max cor:", 
				format(max( cmat[upper.tri(matrix(NA, NN, NN))] ), digits=3), sep="")
				
png(paste("Cor", lev, ".png", sep=""), 2000, 2000,res=300)
plot(data.frame(log(mg.rat[,Group == lev, drop=FALSE])), main=main, pch=20)
dev.off()

}

############################################
# unsupervised analysis: clustering and PCA
############################################

png("HeatmapAll.png", 2000, 2000,res=300)
heatmap(as.matrix(log(na.omit(mg.rat))), distfun=cordist,col=colorRampPalette(c("green","black","red"))(100),
	ColSideColors=rainbow(nlevels(Group))[Group], margins=c(10,10))
legend("topright", fill=rainbow(nlevels(Group)), legend=levels(Group))
dev.off()



pca.res <- PCA(log(t(na.omit(mg.rat))), Group, k=5)
z <- pca.res$componentScores

ld = pca.res$componentLoadings
# can also extract the proteins with extreme loadings

mycols= rainbow(nlevels(Group))  # cc[seq(1, length(cc), floor(length(cc)/nlevels(Group)))]


png("PCA3dPlot.png", 2000, 2000, res=300)
s3d <- scatterplot3d(z[,1:3], color = mycols[Group], col.axis=gray(0.85),col.grid="lightblue",
	box = T, angle = 26, pch=(15+(1:nlevels(Group)))[Group] )
s3d$points3d(z[,1:3])
legend("topright", fill=mycols, legend=levels(Group), cex=.8)
text(s3d$xyz.convert(3+z[,1], 3+z[,2], z[,3]), labels = colnames(mg.rat), cex=0.4)
dev.off()

ord.list = list()

png("PCATopLoadings.png", width=2000, height=700, res=300)
par(oma=c(2,1,1,1))
layout(matrix(1:3, nrow=1))

for (ii in 1:3) {
 ord = order(abs(ld[,ii]), decreasing=TRUE)[1:5]
 barplot(sort(ld[ord, ii]), las=2, main=paste("Top loadings PC", ii))
 ord.list[[ii]]=ord
}
dev.off()

png("PCATopLoadingsProteinPatterns.png", width=2500, height=2500, res=300)
par(mar=c(5,2,3,1))
layout(matrix(1:15, nrow=3, byrow=T))
for (ii in 1:3) {
ord = ord.list[[ii]]
for (xx in 1:5) {
boxplot(mg.rat[match(rownames(ld)[ord[xx]], rownames(mg.rat)),] ~ Group, boxwex=0.5, main=rownames(ld)[ord[xx]], col="gray", las=2)
}
}
dev.off()



############################################
# End unsupervised analysis
############################################


################
# ANOVA
################

Anova = rep(NA, nrow(mg.rat))

# compute Group means (in log space, geometric)
data.ag = aggregate(t(mg.rat), by=list(Group=Group), FUN=function(v){exp(mean(log(na.omit(v))))} )
Means = t( data.ag[,-1])
colnames(Means) = paste("Means",data.ag[,1])
MaxFC = apply(Means, 1, FUN=function(v){max(v)/min(v)})


for (i in 1:nrow(mg.rat)) {

v= t(mg.rat[i,])
nna.idx = !is.na(v)

an.res =  try(anova(lm( log(v[nna.idx]) ~ Group[nna.idx, drop=TRUE] ))[1,"Pr(>F)"])

if (!inherits(an.res, "try-error")) Anova[i] = an.res;

}



Anova.adj  = p.adjust(Anova, method = "fdr")
# Anova.idx = ( MaxFC > FCCutoff ) & !is.na(Anova) & (Anova < 0.05)  # modified 15 Feb
Anova.idx =  !is.na(MaxFC) & ( MaxFC > FCCutoff ) & !is.na(Anova) & (Anova < 0.05)



png("Heatmap - Anova DE.png", 2000, 2000, res=300)
hm1 <- heatmap3(as.matrix(na.omit(log(mg.rat[Anova.idx,]))),  margins=c(10,15), cexRow=1,
	 col=colorRampPalette(c("green", "black", "red"))(120), 
	ColSideColors=rainbow(nlevels(as.factor(Group)))[as.factor(Group)]  )
legend("topright", fill=rainbow(nlevels(as.factor(Group))), legend=levels(as.factor(Group)))
dev.off()



cluster.data = na.omit(log(mg.rat[Anova.idx,]))
NotOmitted = setdiff(1:sum(Anova.idx), attr(cluster.data, "na.action"))

gp = Group



res1 <- try(HClust(cluster.data, metric="pearsonCorrelation", method="complete", cutNumber=4))
clustID <- res1$clustID

Cluster = rep(NA, sum(Anova.idx))
Cluster[NotOmitted] = clustID

noClusters <- 4

r.temp <- aggregate(t(cluster.data), by=list(gp=gp), FUN=mean)
ag.sample <- r.temp[,-1]
rownames(ag.sample) <- r.temp[,1]
ag.genes <- aggregate(t(ag.sample), by=list(Cluster=clustID), FUN=mean)
ag.sd <- aggregate(t(ag.sample), by=list(Cluster=clustID), FUN=sd)
ag.matrix <- as.matrix(ag.genes[,-1])
ag.counts <- summary(as.factor(clustID))
ag.bars <- as.matrix(ag.sd[,-1])

png("ClusterPatterns.png", 2000, 2000, res=300)

par(bg=gray(.95), fg=gray(0.3), oma= c(5, 2, 2, 1) + 0.1, col.main="black", col.sub="black", col.lab="black", col.axis="black")
layout(matrix(1:4, ncol=2, byrow=TRUE))
NSig <- noClusters
for(i in 1:NSig) {
cols <- rainbow(4) 
# cols <- rep("gray", 6)
gname <-  paste("Cluster", i, "(", ag.counts[i], "proteins )")
lines <- ag.sample[, clustID==i, drop=FALSE]
plotErrorBarsLines(ag.matrix[i,], 2*ag.bars[i,], lines, 
labels=1:ncol(ag.matrix), col=cols[i],  main=gname, # bgcol="gray", split=split,
xlab="", ylab="Average log ratio", 
ylim=c(min(ag.matrix), max(ag.matrix)))
axis(1,at=1:ncol(ag.matrix), las=2, labels=colnames(ag.matrix), col="black")
abline(h=0, lty="dotted")
}

dev.off()


Clusters = rep(NA, nrow(full.ratios.only))
Clusters[Anova.idx] = Cluster

################
# End ANOVA
################

TotalCounts = apply(mg.Counts[,-1, drop=FALSE], 1, FUN=function(x){sum(na.omit(x))})


full.ratios.only = data.frame(full.ratios.only, Anova, Anova.adj, MaxFC, Clusters, TotalCounts)

#########################################################
# all ratios changing to reference, with 1-sample t-test
# and composite z-score
#########################################################

# need to bring in z-score and counts
# Outputs should be for each group level - mean ratio; mean z-score; protein p-value where it can be constructed; Total Counts

CompToRef.list = list()

for (levname in levels(Group)) {


which.idx = (Group == levname)

rat.comp = mg.rat[,which.idx, drop=FALSE]
var.comp = mg.var[,-1][, which.idx, drop=FALSE]


# analyses: 1.t-test, 2.cutoff-based, 3.cutoff-count, 4.cutoff-zcore

# FC
FC = apply(rat.comp, 1, FUN=function(v){ w = log(na.omit(v)); exp(mean(w)) })

# z-scores anc composite
zscoreMatrix = 100*log(rat.comp)/(var.comp + 1)  # avoid 0 variance
AverageZ = apply(zscoreMatrix, 1, FUN=function(v){mean(na.omit(v))})
data.frame(rat.comp,var.comp, FC, AverageZ)[order(AverageZ)[1:10],]


# OneSplTTest
OneSplTTest = apply(rat.comp, 1, FUN=function(v){ res = NA; w = log(na.omit(v)); 
								if (length(w) > 1) { tp = try(t.test(w)$p.value); if (!inherits(tp, "try-error")) res=tp }; 
								res })
OneSplTTestAdj = p.adjust(OneSplTTest, method="fdr")


if (sum(!is.na(OneSplTTest)) > 1) {  # addded 15 Feb

png(paste("Volcano ", levname, ".png", sep=""))
plot(log(FC), -log(OneSplTTest), pch=20, main=levname)
abline(h=-log(0.05), col="red")
abline(v=log(1.5), col="blue")
abline(v=-log(1.5), col="blue")
dev.off()

}

sum( !is.na(OneSplTTest) & (abs(log(FC)) > log(FCCutoff)) & (OneSplTTest < 0.05))
sum( !is.na(OneSplTTest) & (abs(log(FC)) > log(FCCutoff)) & (OneSplTTest < 0.05)& !is.na(AverageZ) & (abs(AverageZ) > 2))

sum( !is.na(OneSplTTest) & (abs(log(FC)) > log(FCCutoff)) & (OneSplTTest < 0.05))/nrow(rat.comp)
sig.idx = !is.na(OneSplTTest) & (abs(log(FC)) > log(FCCutoff)) & (OneSplTTest < 0.05)

if (sum(sig.idx) > 1) {
png(paste("SigCor ", levname, ".png", sep=""))
plot(data.frame(log(rat.comp[sig.idx,])), main=levname)
dev.off()
}


CompToRef.list[[levname]] = list(data.frame(FC, OneSplTTest,OneSplTTestAdj, AverageZ), sig.idx)

}


#########################################################
# Venn Diagram overlap or barplot
#########################################################

IDlist = lapply(CompToRef.list, FUN=function(x){mg.Counts[x[[2]],1]})
listLen = sapply(CompToRef.list, FUN=function(x){sum(x[[2]])})

if (length(listLen) > 3) {
	listVen = IDlist[order(listLen, decreasing=T)[1:3]]
} else {
	listVen = IDlist
}


v1 = venn.diagram(listVen[c(1,2)], fill=c("green", "yellow"),  filename=NULL, cex = 1.5, cat.cex = 1.2)
# select 3 largest, plot all overlaps
# single condition check
if (length(listVen) == 2) {

	png("VennOverlap.png")
	grid.draw(v1)
	dev.off() }

if (length(listVen) == 3) {
v2 = venn.diagram(listVen[c(1,3)], fill=c("green", "cornflowerblue"),  filename=NULL, cex = 1.5, cat.cex = 1.2)
v3= venn.diagram(listVen[c(2,3)], fill=c("yellow", "cornflowerblue"),  filename=NULL, cex = 1.5, cat.cex = 1.2)

png("VennOverlap3.png", width=3000, height=1000, res=300)

gl <- grid.layout(nrow=1, ncol=3)

# setup viewports
vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1) 
vp.2 <- viewport(layout.pos.col=2, layout.pos.row=1) 
vp.3 <- viewport(layout.pos.col=3, layout.pos.row=1) 

# init layout
pushViewport(viewport(layout=gl))
# access the first position
pushViewport(vp.1)
grid.draw(v1)
popViewport()

pushViewport(vp.2)
grid.draw(v2)
popViewport()

pushViewport(vp.3)

grid.draw(v3)
popViewport()
dev.off()

}

#######################################
# Output all results to a spreadsheet
#######################################

pairwise.res = data.frame(lapply(CompToRef.list, FUN=function(x){x[[1]]}))
pairwise.res = data.frame(mg.Counts, pairwise.res, mg.rat)


full.res = data.frame(full.ratios.only, Means, full.res)


# including ratios in order of levels
# include all data
# Anova results and clusters
# separate comparison results in each tab


tabnames = names(CompToRef.list)
# remove any odd characters and spaces

wb <- loadWorkbook("ResultsPairwise.xlsx", create = TRUE)

for (jjj in 1:length(CompToRef.list)) {
printStyle(pairwise.res[CompToRef.list[[jjj]][[2]],], ratios=grep("FC", names(pairwise.res)), 
				pvals=c(grep("Anova", names(pairwise.res)),grep("OneSplTTest", names(pairwise.res))), 
				wb = wb, 
                tabName = tabnames[jjj], hiCutoff=FCCutoff, lowCutoff=1/FCCutoff, pvalCutoff=0.05)

}
saveWorkbook(wb)


wb <- loadWorkbook("ResultsOverall.xlsx", create = TRUE)

printStyle(full.res, ratios=c(grep("MaxFC", names(full.res)), grep("Means", names(full.res))), 
				pvals=c(grep("Anova", names(full.res)),grep("OneSplTTest", names(full.res))), 
				wb = wb, 
                tabName = "AllData", hiCutoff=FCCutoff, lowCutoff=1/FCCutoff, pvalCutoff=0.05)

printStyle(data.frame(rownames(pca.res$componentScores),pca.res$componentScores), ratios=NULL, pvals=NULL, wb = wb, tabName = "PCAScores")

printStyle(data.frame(rownames(ld), ld), ratios=NULL, pvals=NULL, wb = wb, tabName = "PCALoadings")


saveWorkbook(wb)






} # end main function

# test
if (FALSE) {

setwd( "X:/Projects/TMT/BookChapter/multiple quality results 126 ref")
runOverallMultivarJobTMT("-ffiles.zip", "-ddesignRuns.xlsx")

}



