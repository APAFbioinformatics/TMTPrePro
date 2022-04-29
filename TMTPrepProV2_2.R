
####################
# TMTPrePro for PD2.2
# JW 28/08/2018
####################



TMTPrepProV2_GP <- function(...) {

	runOverallAndTargetMultivarJobTMT(...);
	
	
}

runOverallAndTargetMultivarJobTMT <- function(...) {

	##############################################
	# main function for overall and targeted analysis
	##############################################

	library(heatmap3)
	library(lattice)
	library(openxlsx)
	library(scatterplot3d) 
	library(TMTPrePro)

	############
	# parameters
	############

	# zipfname = "Raw data.zip"
	# designfname = "Design_Runs1-3.xlsx"
		FCCutoff = 1.2
		PvalCutoff = 0.05
		CountsCutoff = 0
		Clean <- TRUE
		MasterFilter = 'IsMasterProtein'
		
		args <- list(...)

	for(i in 1:length(args)) {
		flag <- substring(args[[i]], 0, 2)
		value <- substring(args[[i]], 3, nchar(args[[i]]))

		if(flag=='-f') zipfname <- value;
		if(flag=='-d') designfname <- value;
		if(flag=='-r') FCCutoff <- as.numeric(value);
		if(flag=='-c') CountsCutoff <- as.numeric(value);
		if(flag=='-p') PvalCutoff <- as.numeric(value);
		if(flag=='-l') Clean <- ifelse(tolower(value)=='yes', TRUE, FALSE)
		if(flag=='-m') MasterFilter <- value;

	} 

	dat.para = data.frame(Parameter=c('P value', 'Fold change', 'Clean'), Cutoff=c(PvalCutoff, FCCutoff, as.character(Clean) ))
	
	

	# unzip files and check zip
	files <- unzip(zipfname, junkpaths=TRUE)
	if (length(files) < 1) stop("The zip file is empty")


	# designfname = "designBradjw_onerun.xlsx"  #"designBradjw.xlsx"
		
	# Error checking, have two or more tabs etc?	

	designSheets = try(loadWorkbook(designfname) )
	if(inherits(designSheets, 'try-error')) stop("Error with loading design file.")

	if(length(names(designSheets)) < 2) stop("Design file needs to have at least 2 tabs.")

	designAll = try( readWorkbook(designfname, 1) )
	if (inherits(designAll, "try-error")) stop("Could not extract the first tab from the design file.")

	references = readWorkbook(designfname, 2)
	if (inherits(references, "try-error")) stop("Could not extract the second tab from the design file.")

	files = readWorkbook(designfname, 2)[,1]
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

	file.list = lapply(files, FUN=read.delim, as.is=T)

	file.list.clean = file.list

	ratios.list = list()

	#######
	# for each file
	#######
	for (fidx in 1:length(file.list)) {


		dat1 = file.list[[fidx]]
		design = design.list[[fidx]]
		referenceLabel = refLabels[[fidx]]
		referenceGroup = refGroups[[fidx]]

		dat1 = dat1[!duplicated(dat1$Accession),]
		
		# data check & cleaning
		
		if(!'Master' %in% colnames(dat1) || !sum(grepl('Protein.FDR.Confidence', colnames(dat1)) ) )
			stop('Must have Master and Protein.FDR.Confidence columns')
		
		if(Clean==TRUE) dat1 = dat1[dat1$Master == MasterFilter & 
		                              dat1[,grep('Protein.FDR.Confidence', colnames(dat1))] == "High", ]
		
		
		file.list.clean[[fidx]] = dat1
		
		if(nrow(dat1)==0) stop("No data after cleaning! Please check cleaning condition!")
		
		
		### No "count" in PD2.2 - use #Peptide as alternative
		
		counts.idx = grep("X\\.+Peptides$", names(dat1))
		# var.idx = setdiff(grep("Variability", names(dat1)), grep("^New", names(dat1)))

		# plot(dat1[,counts.idx])
		
		
		Counts1 = dat1[,counts.idx] 

		# plot(dat1[,counts.idx])

		# extract ratios only
		# exclude counts, variability, "New" or AA to end, if exists
		
		iidx.rat = grepl("Abundance.Ratio", colnames(dat1)) & ! grepl("Adj", colnames(dat1))

		ratios1 = (dat1[,iidx.rat])


		for (i in 1:ncol(ratios1)) {ratios1[,i] = as.numeric(ratios1[,i])}

		# median normalize 
		med = apply(ratios1, 2, FUN=median, na.rm=T)
		rat.norm = sweep(ratios1, 2, med, FUN="/")
		
		# layout(matrix(1:2, ncol=1)) 
		# boxplot(log(ratios1), main='Plain log ratios')
		# boxplot(log(rat.norm), main='Normalised log ratios (medium)' )

		raw.ratios1 = ratios1
		ratios1 = rat.norm

		# substitute multiple dots to single dot in ratio names
		names(ratios1) = gsub("(\\.)+", "\\.", names(ratios1))
		names(ratios1) = gsub("\\.$", "", names(ratios1))

		# remove F? in the ratio names -- added 8 May 2018
		names(ratios1) = gsub("(F\\d+\\.)", "", names(ratios1))
		
		
		
		# remove "_" from design
		# extract labels

		# extract top and bottom ratios
		top1 = gsub("Abundance.Ratio.(.*)\\.(.*)", "\\1", colnames(ratios1))
		bot1 = gsub("Abundance.Ratio.(.*)\\.((.*))", "\\2", colnames(ratios1))

		GroupTop = as.factor(design$Group[match(top1, gsub("_", "", design$Label))])
		GroupBot = as.factor(design$Group[match(bot1, gsub("_", "", design$Label))])
		
		
		ref.idx = ( bot1 == gsub("_", "", referenceLabel))  #              

		ratiosRef = ratios1[,ref.idx]
		# varRef = Var1[, ref.idx]
		topRef = top1[ref.idx]

	
		
		ratiosRef = data.frame(dat1[,"Accession"], ratiosRef)
		
		ratios1 = data.frame(dat1[,"Accession"], ratios1) # added on 12/04

		# varRef = data.frame(dat1[,"Accession"],varRef)
		GroupRef = GroupTop[ref.idx]

		
		# add variability and counts

		ratios.list[[fidx]] = list(Ratios = ratios1, ratiosRef = ratiosRef, Group=GroupRef, GroupTop=GroupTop, GroupBot=GroupBot, Counts=data.frame(dat1[,"Accession"],Counts1))


	}
	 # end for all files


	# merge ratios to selected reference from all files

	all.fullrat = ratios.list[[1]]$Ratios
	colnames(all.fullrat) = paste("R1", colnames(all.fullrat),sep="")

	mg.rat = ratios.list[[1]]$ratiosRef
	# mg.var = ratios.list[[1]]$Var
	mg.Counts = ratios.list[[1]]$Counts
	Group = as.factor(as.vector(ratios.list[[1]]$Group))
	colnames(mg.rat) = paste("R1", colnames(mg.rat))


	if (length(ratios.list) > 1) {

		for (xx in 2:length(ratios.list)) {

			nextrun.fullrat = ratios.list[[xx]]$Ratios
			colnames(nextrun.fullrat) = paste("R", xx, colnames(nextrun.fullrat), sep="")
			
			nextrun.rat = ratios.list[[xx]]$ratiosRef
			colnames(nextrun.rat) = paste("R", xx, colnames(nextrun.rat), sep="")

			nextCount = ratios.list[[xx]]$Counts
			colnames(nextCount) = paste("R", xx, colnames(nextCount), sep="")

			# merge runs by accession; commented out on 12/04/18
			# fullrat = merge(all.fullrat, nextrun.fullrat, by.x=1, by.y=1, all=TRUE)
			
			mg = merge(mg.rat, nextrun.rat , by.x=1, by.y=1, all=TRUE)
			  # vvar = merge(mg.var, nextVar , by.x=1, by.y=1, all=TRUE)
			counts = merge(mg.Counts, nextCount, by.x=1, by.y=1, all=TRUE)
			Gp = as.factor(c(as.vector(Group), as.vector(ratios.list[[xx]]$Group) ) )
			
			#  commented out on 12/04/18
			# all.fullrat = fullrat
			mg.rat = mg
			 #  mg.var = vvar
			mg.Counts = counts
			Group = Gp
		  
		}
	}


	# check ordering same and not modified by merging
	# if (sum(mg.rat[,1] != mg.var[,1]) > 0) stop("Error in merging runs: ratios and var mismatch");
	if (sum(mg.rat[,1] != mg.Counts[,1]) > 0) stop("Error in merging runs: ratios and counts mismatch");


	# remove identifier column from ratios

	mg.rat[,1] = as.vector(mg.rat[,1])


	rownames(mg.rat) = mg.rat[,1]

	###! mg.rat.copy = mg.rat
	mg.rat = as.matrix(mg.rat[,-1])

	###! for (i in 2:ncol(mg.rat.copy)) {mg.rat.copy[,i] = as.numeric(mg.rat.copy[,i])};
	###! mg.rat = as.matrix(mg.rat.copy[,-1])

	cat(paste('ncol(mg.rat)=', ncol(mg.rat), '\n' ) )

	for(ii in 1:ncol(mg.rat)) mg.rat[,ii] = as.numeric(mg.rat[,ii])
	
	# reorder ratio columns by order of Group
	full.ratios.only = data.frame(mg.Counts, mg.rat[,order(Group)])
	# full.res = data.frame(mg.Counts, mg.var[, c(1,1+order(Group))])
	full.res = data.frame(mg.Counts)
	
	# Modified on 22/08/18 - output all samples and corresponding groups in the new order
	write.csv(data.frame(Sample=colnames(mg.rat[,order(Group)]), Group=Group[order(Group)]), "samplegroup.csv")
	
	
	
	

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
	legend("topright", fill=rainbow(nlevels(Group)), legend=levels(Group), xpd=TRUE, cex=.6 )
	dev.off()


	png("BoxplotDensity.png", width=3500, height=1700,res=300)
	layout(matrix(1:2, nrow=1))

	par(mar=c(13,4,4,2)+.1)
	plotDensity(t(log(na.omit(mg.rat))), Group)

	# boxplots and density plots
	boxplot(log(mg.rat[, order(Group)]), las=2, col=rainbow(nlevels(Group))[Group[order(Group)]], 
		ylim=c(-1.5,1.5), cex.axis=0.6)
	

	dev.off()


	cat('Begin correlations for all groups\n')
	
	
	# Correlations for all levels of the group
	for (lev in levels(Group)) {
	
		dd = mg.rat[,Group == lev, drop=FALSE]
		
		if(ncol(dd) > 1) {
			cmat = try(cor(dd,method="spearman", use="pairwise.complete.obs"))
			NN = nrow(cmat)
			main = paste(lev, "; min cor:", format(min( cmat[upper.tri(matrix(NA, NN, NN))] ), digits=3), ", max cor:", 
							format(max( cmat[upper.tri(matrix(NA, NN, NN))] ), digits=3), sep="")
							
			png(paste("Cor", lev, ".png", sep=""), 2000, 2000,res=300)
			plot(data.frame(log(mg.rat[,Group == lev, drop=FALSE])), main=main, pch=20)
			dev.off()
		}
	}

	cat('End correlations for all groups\n')
	
	############################################
	# unsupervised analysis: clustering and PCA
	############################################

		cat('Begin unsupervised\n')
		
	png("HeatmapAll.png", 2000, 2000,res=300)
	heatmap3(as.matrix(log(na.omit(mg.rat))), distfun=cordist,col=colorRampPalette(c("green","black","red"))(100),
		ColSideColors=rainbow(nlevels(Group))[Group], margins=c(10,10))
	legend("topright", fill=rainbow(nlevels(Group)), legend=levels(Group),  xpd=TRUE, cex=.6 )
	dev.off()



	pca.res <- PCA(log(t(na.omit(mg.rat))), Group, k=5)
	z <- pca.res$componentScores
	
	vars = apply(z, 2, var) 
	
	props = round(100*vars/sum(vars), 1)
	
	

	ld = pca.res$componentLoadings
	# can also extract the proteins with extreme loadings

	mycols= rainbow(nlevels(Group))  # cc[seq(1, length(cc), floor(length(cc)/nlevels(Group)))]


	png("PCA3dPlot.png", 2000, 2000, res=300)
	s3d <- scatterplot3d(z[,1:3], color = mycols[Group], col.axis=gray(0.85),col.grid="lightblue",
		box = T, angle = 26, pch=20 )
	s3d$points3d(z[,1:3], pch=21)
	legend("topright", fill=mycols, legend=levels(Group), cex=.8)
	text(s3d$xyz.convert(3+z[,1], 3+z[,2], z[,3]), labels = colnames(mg.rat), cex=0.4)
	dev.off()

	ord.list = list()



	png("PCA12.png", 2000, 2000, res=300)

	plot(z[,1], z[,2], col=rainbow(nlevels(as.factor(Group)))[as.factor(Group)], pch=20, 
	     xlab=paste0("PC1 (", props[1], "%)"),
	     ylab=paste0("PC2 (", props[2], "%)"))
	points(z[,1], z[,2], pch=21, cex=1.1, lwd=1.3)
	text(z[,1], z[,2], colnames(mg.rat), pos=3, cex=.5)
	legend("bottomright", fill=rainbow(nlevels(as.factor(Group))), legend=levels(as.factor(Group)),
		 xpd=TRUE, inset=c(-.13,0), cex=.6 )

	dev.off()

	png("PCA13.png", 2000, 2000, res=300)

	plot(z[,1], z[,3], col=rainbow(nlevels(as.factor(Group)))[as.factor(Group)], pch=20, 
	     xlab=paste0("PC1 (", props[1], "%)"),
	     ylab=paste0("PC3 (", props[3], "%)"))
	points(z[,1], z[,3], pch=21, cex=1.1, lwd=1.3)
	text(z[,1], z[,3], colnames(mg.rat), pos=3, cex=.5)
	legend("topright", fill=rainbow(nlevels(as.factor(Group))), legend=levels(as.factor(Group)),
		 xpd=TRUE, inset=c(-.13,0), cex=.6 )

	dev.off()

	png("PCA23.png", 2000, 2000, res=300)

	plot(z[,2], z[,3], col=rainbow(nlevels(as.factor(Group)))[as.factor(Group)], pch=20, 
	     xlab=paste0("PC2 (", props[2], "%)"),
	     ylab=paste0("PC3 (", props[3], "%)"))
	points(z[,2], z[,3], pch=21, cex=1.1, lwd=1.3)
	text(z[,2], z[,3], colnames(mg.rat), pos=3, cex=.5)
	legend("topright", fill=rainbow(nlevels(as.factor(Group))), legend=levels(as.factor(Group)),
	 xpd=TRUE, inset=c(-.13,0), cex=.6 )

	dev.off()
	
	
	png("PCA2DAll.png", 2000, 2000, res=300)
	layout(matrix(1:4, ncol=2))
	plot(z[,1], z[,2], col=rainbow(nlevels(as.factor(Group)))[as.factor(Group)], pch=20, xlab=paste0("PC1 (", props[1], "%)"),
	     ylab=paste0("PC2 (", props[2], "%)"))
	points(z[,1], z[,2], pch=21, cex=1.1, lwd=1.3)
	text(z[,1], z[,2], colnames(mg.rat), pos=3, cex=.5)
	
	plot(z[,1], z[,3], col=rainbow(nlevels(as.factor(Group)))[as.factor(Group)], pch=20, xlab=paste0("PC1 (", props[1], "%)"),
	     ylab=paste0("PC3 (", props[3], "%)"))
	points(z[,1], z[,3], pch=21, cex=1.1, lwd=1.3)
	text(z[,1], z[,3], colnames(mg.rat), pos=3, cex=.5)
	
	plot(z[,2], z[,3], col=rainbow(nlevels(as.factor(Group)))[as.factor(Group)], pch=20, xlab=paste0("PC2 (", props[2], "%)"),
	     ylab=paste0("PC3 (", props[3], "%)"))
	points(z[,2], z[,3], pch=21, cex=1.1, lwd=1.3)
	text(z[,2], z[,3], colnames(mg.rat), pos=3, cex=.5)
	
	legend("topright", fill=rainbow(nlevels(as.factor(Group))), legend=levels(as.factor(Group)),
	       xpd=TRUE, inset=c(-.13,0), cex=.6 )
	
	dev.off()
	
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


	cat('End unsupervised \n')
	############################################
	# End unsupervised analysis
	############################################


	################
	# ANOVA
	################

		cat('Begin anova\n')
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
	# Anova.idx = ( MaxFC > FCCutoff ) & !is.na(Anova) & (Anova < 0.05)  
	Anova.idx =  !is.na(MaxFC) & ( MaxFC > FCCutoff ) & !is.na(Anova) & (Anova < PvalCutoff)


	if(nrow(mg.rat[Anova.idx,]) > 3) {
	png("Heatmap - Anova DE.png", 2000, 2000, res=300)
	hm1 <- heatmap3(as.matrix(na.omit(log(mg.rat[Anova.idx,]))),  margins=c(15,15), cexRow=1,
		 col=colorRampPalette(c("green", "black", "red"))(120), 
		ColSideColors=rainbow(nlevels(as.factor(Group)))[as.factor(Group)]  )
	legend("topright", fill=rainbow(nlevels(as.factor(Group))), legend=levels(as.factor(Group)),
	 xpd=TRUE,cex=.6 )
	dev.off()
	}
		cat('End anova\n')

				cat('Begin clustering\n')
	cluster.data = na.omit(log(mg.rat[Anova.idx,]))
	NotOmitted = setdiff(1:sum(Anova.idx), attr(cluster.data, "na.action"))

	gp = Group


	Cluster = rep(NA, sum(Anova.idx))
	
	res1 <- try(HClust(cluster.data, metric="pearsonCorrelation", method="complete", cutNumber=4))
	
	if(!inherits(res1, "try-error")) {
	clustID <- res1$clustID

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
		labels=1:ncol(ag.matrix), 
		col=cols[i],  main=gname, # bgcol="gray", split=split,
		xlab="", ylab="Average log ratio", 
		ylim=c(min(ag.matrix), max(ag.matrix)))
		axis(1,at=1:ncol(ag.matrix), las=2, labels=colnames(ag.matrix), col="black")
		abline(h=0, lty="dotted")
	}

	dev.off()
	}

	Clusters = rep(NA, nrow(full.ratios.only))
	Clusters[Anova.idx] = Cluster

	cat('end clustering\n')
	################
	# End ANOVA
	################

	# TotalCounts = apply(mg.Counts[,-1, drop=FALSE], 1, FUN=function(x){sum(na.omit(x))})


	full.ratios.only = data.frame(full.ratios.only, Anova, Anova.adj, MaxFC, Clusters)
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
		# var.comp = mg.var[,-1][, which.idx, drop=FALSE]


		# analyses: 1.t-test, 2.cutoff-based, 3.cutoff-count, 4.cutoff-zcore

		# FC
		FC = apply(rat.comp, 1, FUN=function(v){ w = log(na.omit(v)); exp(mean(w)) })

		# z-scores anc composite
		# zscoreMatrix = 100*log(rat.comp)/(var.comp + 1)  # avoid 0 variance
		# AverageZ = apply(zscoreMatrix, 1, FUN=function(v){mean(na.omit(v))})
		# data.frame(rat.comp,var.comp, FC, AverageZ)[order(AverageZ)[1:10],]


		# OneSplTTest
		OneSplTTest = apply(rat.comp, 1, FUN=function(v){ res = NA; w = log(na.omit(v));
										if (length(w) > 1) { tp = try(t.test(w)$p.value); if (!inherits(tp, "try-error")) res=tp }; 
										res })
		
										
		OneSplTTestAdj = p.adjust(OneSplTTest, method="fdr")


		if (sum(!is.na(OneSplTTest)) > 1) { 

			png(paste("Volcano ", levname, ".png", sep=""))
			plot(log(FC), -log(OneSplTTest), pch=20, main=levname)
			abline(h=-log(PvalCutoff), col="red")
			abline(v=log(FCCutoff), col="blue")
			abline(v=-log(FCCutoff), col="blue")
			dev.off()

		}

		# sum( !is.na(OneSplTTest) & (abs(log(FC)) > log(FCCutoff)) & (OneSplTTest < 0.05))
		# sum( !is.na(OneSplTTest) & (abs(log(FC)) > log(FCCutoff)) & (OneSplTTest < 0.05)& !is.na(AverageZ) & (abs(AverageZ) > 2))

		# sum( !is.na(OneSplTTest) & (abs(log(FC)) > log(FCCutoff)) & (OneSplTTest < 0.05))/nrow(rat.comp)
		sig.idx = !is.na(OneSplTTest) & (abs(log(FC)) > log(FCCutoff)) & (OneSplTTest < PvalCutoff)

		if (sum(sig.idx) > 1) {
			png(paste("SigCor ", levname, ".png", sep=""))
			plot(data.frame(log(rat.comp[sig.idx,])), main=levname)
			dev.off()
		}


		CompToRef.list[[levname]] = list(data.frame(FC, OneSplTTest,OneSplTTestAdj), sig.idx)
	}

	#######################################
	# Output all results to a spreadsheet
	#######################################

	pairwise.res = data.frame(lapply(CompToRef.list, FUN=function(x){x[[1]]}))


	# add annotation columns from raw data

	keep.colums = c("Description", 
	                "Biological.Process",
	                "Cellular.Component",
	                "Molecular.Function",
	                "Pfam.IDs",
	                "Entrez.Gene.ID",
	                "Ensembl.Gene.ID",
	                "Gene.Symbol",
	                "Chromosome",
	                "KEGG.Pathways",
	                "WikiPathways")
	
	dat.all = data.frame(mg.Counts)
	
	for(kk in 1:length(keep.colums)) {
	  
	  if(sum(grepl(keep.colums[kk], colnames(file.list.clean[[1]]))) ==1) { 
	    
	    temp = getAnnotation(dat.all[,1], file.list.clean, keep.colums[kk])
	  
	    dat.all = cbind(dat.all, temp[,-1, drop=FALSE])
	  }
	  
	    
	}
	
	

	###! pairwise.res = data.frame(dat.all, mg.Counts, pairwise.res, mg.rat)
	###! pairwise.res = data.frame(dat.all, mg.Counts[,-1], pairwise.res)
	
	pairwise.res = data.frame(dat.all, pairwise.res)
	
	colnames(pairwise.res)[1] = 'Accession'

	###! full.res = data.frame(dat.all, full.ratios.only, Means, full.res)

	full.res = data.frame(dat.all, full.ratios.only[,-1], Means)
	colnames(full.res)[1] = 'Accession'
	
	# including ratios in order of levels
	# include all data
	# Anova results and clusters
	# separate comparison results in each tab


	tabnames = names(CompToRef.list)
	# remove any odd characters and spaces

	wb <- createWorkbook()

	# add parameter cutoffs
	addWorksheet(wb, sheet='Parameter')
	writeData(wb, 'Parameter', dat.para)
	
	for (jjj in 1:length(CompToRef.list)) {

		# only the DE proteins?
		dd = pairwise.res[CompToRef.list[[jjj]][[2]],]		

		if(nrow(dd) > 1) {
			dd[dd=='NaN'] = NA

			ps <- try(printOpenxlsxStyle(dd, ratios=grep("FC", names(pairwise.res)), 
							pvals=c(grep("Anova", names(pairwise.res)),grep("OneSplTTest", names(pairwise.res))), 
							wb = wb, 
							tabName = tabnames[jjj], hiCutoff=FCCutoff, lowCutoff=1/FCCutoff, pvalCutoff=PvalCutoff))
			if(inherits(ps, 'try-error') ) warning('Error with print comparison to reference xlsx file')
		}
	}
	saveWorkbook(wb, "ResultsPairwise.xlsx", overwrite=TRUE)


	wb <- createWorkbook()
	full.res[full.res=='NaN'] = NA

	# add parameter cutoffs
	addWorksheet(wb, sheet='Parameter')
	writeData(wb, 'Parameter', dat.para)
		
	
	ps <- try(printOpenxlsxStyle(full.res, ratios=c(grep("MaxFC", names(full.res)), grep("Means", names(full.res))), 
					pvals=c(grep("Anova", names(full.res)),grep("OneSplTTest", names(full.res))), 
					wb = wb, 
					tabName = "AllData", hiCutoff=FCCutoff, lowCutoff=1/FCCutoff, pvalCutoff=PvalCutoff) )

	if(inherits(ps, 'try-error') ) warning('Error with print overall xlsx file')
					
	ps <- try(printOpenxlsxStyle(data.frame(rownames(pca.res$componentScores),pca.res$componentScores), ratios=NULL, pvals=NULL, wb = wb, tabName = "PCAScores"))

	if(inherits(ps, 'try-error') ) warning('Error with print overall PCA component scores tab')

	ps <- try(printOpenxlsxStyle(data.frame(rownames(ld), ld), ratios=NULL, pvals=NULL, wb = wb, tabName = "PCALoadings"))

	if(inherits(ps, 'try-error') ) warning('Error with print overall PCA loading tab')


	saveWorkbook(wb, file="ResultsOverall.xlsx", overwrite=TRUE)


	#output design
	saveWorkbook(designSheets, file='Design.xlsx', overwrite=TRUE)
	
	

	########################################
	###!!! Targeted pairwise comparison 
	########################################

	if(sum(grepl("comp", tolower(names(designSheets)) )) == 1 ) {

		# use ratios to reference		
		
		dat.comparisons = readWorkbook(designfname, 3)

		tarcompres.list = list()

		for(idx.comp in 1:nrow(dat.comparisons) ) {

			comp = dat.comparisons[idx.comp,-1]

			mat.rat =  mg.rat
			
			idx.group1 = which(dat.comparisons[idx.comp,2] == as.character(Group))

			idx.group2 = which(dat.comparisons[idx.comp,3] == as.character(Group))

			idx.mean1 = which(dat.comparisons[idx.comp,2] == gsub('Means ', '', colnames(Means), fixed=TRUE) ) 
			  
			idx.mean2 = which(dat.comparisons[idx.comp,3] == gsub('Means ', '', colnames(Means), fixed=TRUE) )

			if(length(idx.group1) < 2 || length(idx.group2) < 2) {
				cat('Warning: not enough samples for targted comparison\n')
			} else {
				FC = Means[,idx.mean1]/Means[,idx.mean2]

				TwoSplTTest = rep(NA, nrow(mg.rat))
				TwoSplTTestAdj = rep(NA, nrow(mg.rat))
				
				for(ii in 1:nrow(mat.rat)) {
					
					v= t(mat.rat[ii,])				

					# two sample t-test
					
					ttest.res =  try(t.test(log(na.omit(v[idx.group1])), 
									log(na.omit(v[idx.group2])) , var.equal=TRUE ) )

					if (!inherits(ttest.res, "try-error")) TwoSplTTest[ii] = ttest.res$p.value;
					
					
				}
				
				TwoSplTTestAdj = p.adjust(TwoSplTTest, method="fdr")
				
				idx.sig = (!is.na(TwoSplTTest)) & (abs(log(FC)) > log(FCCutoff)) & (TwoSplTTest < PvalCutoff)


				# volcano plot

				png(paste('Volcano plot for targeted', dat.comparisons[idx.comp,1], '.png', sep=''), 2000, 2000, res=300)
				plot(log(FC), -log(TwoSplTTest), xlab='log FC', ylab='-log p value', 
							main=paste('Protein volcano plot', dat.comparisons[idx.comp,1]))
				abline(h=-log(PvalCutoff), col="red")
				abline(v=log(FCCutoff), col="blue")
				abline(v=-log(FCCutoff), col="blue")
				
				if(sum(idx.sig) > 0) 			
					points(log(FC[idx.sig]), -log(TwoSplTTest[idx.sig]), col='red', pch=20)

				dev.off()
				

				# Get annotation from raw
				
				dat.all = data.frame(Accession=rownames(mg.rat))
				
				for(kk in 1:length(keep.colums)) {
				  
				  if(sum(grepl(keep.colums[kk], colnames(file.list.clean[[1]]))) ==1) { 
				    
				    temp = getAnnotation(dat.all[,1], file.list.clean, keep.colums[kk])
				    
				    dat.all = cbind(dat.all, temp[,-1, drop=FALSE])
				  }
				  
				  
				}
				
				

				dat.all = merge(dat.all, mg.Counts, by=1, all=TRUE)
				
				dat.comp = data.frame(Accession=rownames(mg.rat), mat.rat[,c(idx.group1, idx.group2)], Means[,c(idx.mean1, idx.mean2)],
							FC, TwoSplTTest, TwoSplTTestAdj, Significant = idx.sig)
				
				dat.comp = merge(dat.all, dat.comp, by=1, all=TRUE)
				
				tarcompres.list[[idx.comp]] = list(Result=dat.comp, ComparisonName=dat.comparisons[idx.comp,1] )
			}
		
		}

		
		wb <- createWorkbook()
		
		for(idx.comp in 1:length(tarcompres.list)) {
		
			dd = tarcompres.list[[idx.comp]][[1]]
			dd[dd=='NaN'] = NA
			
			colnames(dd)[1] = 'Accession'
			
			# sort by significant
			dd = dd[order(dd$Significant, decreasing=TRUE), ]
			
			printOpenxlsxStyle(dd, ratios=grep('FC', names(tarcompres.list[[idx.comp]][[1]])),
					pvals=grep('TwoSplTTest', names(tarcompres.list[[idx.comp]][[1]])), wb=wb,
					tabName = tarcompres.list[[idx.comp]]$ComparisonName, hiCutoff=FCCutoff, lowCutoff=1/FCCutoff, pvalCutoff=PvalCutoff)
		}
		
		
		# images
		addWorksheet(wb, sheet='images')	
		startCol = 1
		for(idx.comp in 1:length(tarcompres.list)) {
			if(file.exists(paste('Volcano plot for targeted', dat.comparisons[idx.comp,1], '.png', sep=''))) {			
				insertImage(wb, sheet='images', 
					file=paste('Volcano plot for targeted', dat.comparisons[idx.comp,1], '.png', sep=''), 
					width=12, height=15, startRow=2, startCol=startCol, units='cm')
				
				startCol = startCol + 10
			}
		}
		
		# comparisons
		addWorksheet(wb, sheet='comparisons')
		writeData(wb, sheet='comparisons', dat.comparisons)
		
		saveWorkbook(wb, file="ResultsTargetedPaiwise.xlsx", overwrite=TRUE)
	}

###!!! END TARGETED

} # END runOverallTargetJob




# function printOpenxlsxStyle

printOpenxlsxStyle  <- function (dat, ratios, pvals, wb, tabName = "results", hiCutoff = 1.5, lowCutoff=0.67, pvalCutoff=0.05) 
{

	addWorksheet(wb, sheet=tabName)
	
	upReg <- createStyle(fgFill = "violet")
	downReg <- createStyle(fgFill = "forestgreen")
	sigStyle <- createStyle(fgFill = "khaki1")
	
	writeData(wb, tabName, dat, keepNA=FALSE)
	
    for (rat in ratios) {
        up.idx <- which(!is.na(dat[, rat]) & (dat[, rat] > hiCutoff))
        if (length(up.idx) > 1) 
			addStyle(wb, tabName, style=upReg, rows = 1 + up.idx, cols = rat)

        down.idx <- which(!is.na(dat[, rat]) & (dat[, rat] < 
            lowCutoff))
        if (length(down.idx) > 1) 
			addStyle(wb, tabName, style=downReg, rows = 1 + down.idx, cols = rat)
    }
	
    for (pval in pvals) {
        sig.idx <- which(!is.na(dat[, pval]) & (dat[, pval] < 
            pvalCutoff))
        if (length(sig.idx) > 1) 
			addStyle(wb, tabName, style=sigStyle, rows = 1 + sig.idx, cols = pval)
    }
}


# Added on 22/08/18 to fix the overlapped labels
plotErrorBarsLines <- function (v, barSizes, lines, labels = NULL, col = "blue", ylim = c(min(lines), 
    max(lines)), ...) 
{
    barSizes[is.na(barSizes)] <- 0
    topBars <- v + 0.5 * barSizes
    bottomBars <- v - 0.5 * barSizes
    N <- length(v)
    if (is.null(labels)) 
        labels <- 1:N
    ylims <- c(min(bottomBars, ylim[1], min(lines)), max(topBars, 
        ylim[2], max(lines)))
    par(pch = 19, xaxt = "n")
    plot(as.numeric(labels), v, ylim = ylims, col = col, type = "b", 
        lwd = 3, ...)
    par(xaxt = "s")
 
    for (i in 1:N) {
        lines(c(i, i), c(topBars[i], bottomBars[i]))
    }
    for (i in 1:ncol(lines)) {
        lines(as.numeric(labels), lines[, i], lwd = 0.5, lty = "dotted", 
            col = "gray")
    }
}


# Get annotations from raw files

getAnnotation <- function(accessions,  list.dat.files, name.annot = "Description", by.col = "Accession") {
  

  dat.all = data.frame(Accession=accessions, AnnotTemp=rep(NA, length(accessions)))


  for (i in 1:length(list.dat.files)) {
    
    dat.all[match(list.dat.files[[i]][,by.col], dat.all[,by.col]), "AnnotTemp"] = 
      list.dat.files[[i]][, name.annot]
    
  }
  
  colnames(dat.all)[grepl("AnnotTemp", colnames(dat.all))] = name.annot
  
  dat.all
}





# test cases

# Adam Walker


if(FALSE) {
	rm(list=ls())
	setwd( "\\\\APAF-HPV-FILE\\BioInfo_Project\\Projects\\TMT\\TMTPrePro_V2.2\\testset 1 Adam Walker spinal cord")
	source("..\\TMTPrepProV2_2.R")
	TMTPrepProV2_GP("-fRaw data.zip", "-dDesign_Runs1-3 2.xlsx", "-r1.5", "-p0.05", "-lYes")

}















