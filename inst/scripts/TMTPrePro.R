############################################################
# parameters needed
# 	metafile name, by default "meta.csv" in local folder
#	metafile must contains files and groups
# defaults if any
#  	in and out folder by default local folder
#############################################################

TMTPrePro <- function(...) 
{


library(TMTPrePro)
library(openxlsx)


args <- list(...)

# check type: overall or targeted
for(i in 1:length(args)) {
flag <- substring(args[[i]], 0, 2)
value <- substring(args[[i]], 3, nchar(args[[i]]))

if(flag=='-t') jobtype <- value;
if(flag=='-f') filename <- value;
if(flag=='-d') designfile <- value;
if(flag=='-r') FCCutoff <- as.numeric(value);
if(flag=='-z') ZScoreCutoff <- as.numeric(value);
if(flag=='-c') CountsCutoff <- as.numeric(value);
if(flag=='-p') PvalCutoff <- as.numeric(value);

} 


if (value == "OverallMultivar") runOverallMultivarJobTMT(...);
if (value == "Targeted") singleRunJob(pdFile=filename,comparisonInfoFile=designfile,
                      hirat=FCCutoff,lowrat=1/FCCutoff,zlim=ZScoreCutoff,plim=PvalCutoff);



}

