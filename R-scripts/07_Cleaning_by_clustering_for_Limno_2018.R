#Install and load required packages

list.of.packages <- c('scales', 'data.table','randomForest', 'plot3Drgl')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#This package has a different source and needs to be installed separately
source("http://bioconductor.org/biocLite.R")
biocLite("flowPeaks")

#Load packages

library(flowPeaks); library(scales);library(data.table);library(randomForest);library(plot3Drgl)
library(tidyverse)
#Set seed for repeatability
set.seed(42)

#Set your working directory, which contains both the data subset and the files to be cleaned (the latter can be in subfolders)

#Read in the data subset to be used for clustering
df <- read.csv('data/alldatasubset_100k_2018_Limno.csv')
df_raw <- read_csv("data-processed/anc4-particles.csv") 
df <- df_raw %>% 
	sample_n(size = 50000, replace = FALSE)
dat <- df_raw %>% 
	sample_n(size = 100000, replace = FALSE)

##Define new variables to help identify electonic noise
#These essentially quantify how different the first and last values are. 0.1 added because a) there
# are zeroes, b) it is just below the lowest non-zero measurement, c) we need to log transform later
# df$FL.Red.Gradient<- abs(df$FL.Red.First-df$FL.Red.Last)+0.1
# df$X2.FL.Red.Gradient<- abs(df$X2.FL.Red.First-df$X2.FL.Red.Last)+0.1

#Note: we do not really use SWS.Length in clustering because size  estimates at small sizes (<5 microns) are dodgy right now. Worth taking a look at, though, so it's used here for visualization only

#Select only the variables used for clustering and visualization
# dat.clust<-df[,c('FWS.Range','FL.Red.Range','X2.FL.Red.Range', 'FL.Orange.Range',
#                  'FL.Yellow.Range', 'X2.FL.Red.Gradient','FL.Red.Gradient',
#                  'SWS.Length','FWS.Fill.factor')]

dat.clust<-df[,c('fsc_a','ssc_a','fl1_a', 'fl2_a',
				 'fl3_a')]

#Subset data to remove entirely ridiculous cases that will mess with clustering.
dat.clust<-subset(dat.clust,dat.clust$ssc_a > 0)
dat.clust<-subset(dat.clust,dat.clust$fl1_a > 0)
dat.clust<-subset(dat.clust,dat.clust$fl2_a > 0)
dat.clust<-subset(dat.clust,dat.clust$fsc_a > 0)
dat.clust<-subset(dat.clust,dat.clust$fl3_a > 0)

#The next two are more of a judgement call. Take a look at plots and clustering before trying them.
# dat.clust<-subset(dat.clust,dat.clust$FWS.Range > 0.2)
# dat.clust<-subset(dat.clust,dat.clust$SWS.Length > 1.5)

#Take the log of all except for FWS.Fill.factor, which ranges from 0 to 1
# dat.clustlog<-data.frame(log10(dat.clust[,-c(5)]), dat.clust[,c(5)])
dat.clustlog<-data.frame(log10(dat.clust[,c(1:5)]))
# colnames(dat.clustlog)[9]<- 'FWS.Fill.factor'

#Cluster
subflow<-flowPeaks(dat.clustlog[,c(1:5)], tol=0.25, h0 = 0.05, h=2)
#1:5,7,9 works well

# examine details if interested. The 'weight' column indicatest the proportion of the data belonging
# to that cluster
summary(subflow)

# Examine scatterplot matrix coloured by cluster identity if desired. Can be useful, but each
# plot is small
# pairs(dat.clustlog,col=alpha(subflow$peaks.cluster, 0.5), pch='.',lower.panel=NULL,gap=0.2)

# Examine any pairwise plot of interest
# plot(dat.clustlog$FL.Red.Range~dat.clustlog$FWS.Range,col=alpha(subflow$peaks.cluster, 0.4))
# plot(dat.clustlog$FL.Red.Gradient~dat.clustlog$X2.FL.Red.Gradient,col=alpha(subflow$peaks.cluster, 0.7))

#Define colour palette for clusters in 3D plots
palette1<-c('black','red','green3','blue','cyan','magenta','yellow','gray','brown','purple','pink','orange','darkgreen')

#Examine any plots of interest
scatter3Drgl(dat.clustlog$fl1_a,dat.clustlog$fl2_a,dat.clustlog$fl3_a,colvar = subflow$peaks.cluster, col=palette1[1:max(subflow$peaks.cluster)],xlab='fl1', ylab='fl2', zlab='fl3',ticktype="detailed")

scatter3Drgl(dat.clustlog$FL.Red.Range,dat.clustlog$X2.FL.Red.Range, dat.clustlog$X2.FL.Red.Gradient,colvar = subflow$peaks.cluster, col=palette1[1:max(subflow$peaks.cluster)],xlab='FL.Red.Range', ylab='X2.FL.Red.Range', zlab='X2.FL.Red.Gradient')

scatter3Drgl(dat.clustlog$FL.Red.Range,dat.clustlog$X2.FL.Red.Range, dat.clustlog$FL.Red.Gradient,colvar = subflow$peaks.cluster, col=palette1[1:max(subflow$peaks.cluster)],xlab='FL.Red.Range', ylab='X2.FL.Red.Range', zlab='FL.Red.Gradient')

scatter3Drgl(dat.clustlog$FL.Red.Range,dat.clustlog$X2.FL.Red.Range, dat.clustlog$FL.Yellow.Range,colvar = subflow$peaks.cluster, col=palette1[1:max(subflow$peaks.cluster)], xlab='FL.Red.Range', ylab='X2.FL.Red.Range', zlab='FL.Yellow.Range')

scatter3Drgl(dat.clustlog$FL.Red.Range,dat.clustlog$X2.FL.Red.Range, dat.clustlog$FL.Orange.Range,colvar = subflow$peaks.cluster, col=palette1[1:max(subflow$peaks.cluster)], xlab='FL.Red.Range', ylab='X2.FL.Red.Range', zlab='FL.Orange.Range')

scatter3Drgl(dat.clustlog$FL.Red.Range,dat.clustlog$X2.FL.Red.Range, dat.clustlog$SWS.Length,colvar = subflow$peaks.cluster, col=palette1[1:max(subflow$peaks.cluster)], xlab='FL.Red.Range', ylab='X2.FL.Red.Range', zlab='SWS.Length')

scatter3Drgl(dat.clustlog$FWS.Range,dat.clustlog$FL.Red.Range, dat.clustlog$SWS.Length,colvar = subflow$peaks.cluster, col=palette1[1:max(subflow$peaks.cluster)], xlab='FWS.Range', ylab='FL.Red.Range', zlab='SWS.Length')

scatter3Drgl(dat.clustlog$FL.Red.Range,dat.clustlog$FL.Orange.Range, dat.clustlog$FL.Yellow.Range,colvar = subflow$peaks.cluster, col=palette1[1:max(subflow$peaks.cluster)], xlab='FL.Red.Range', ylab='FL.Orange.Range', zlab='FL.Yellow.Range')

#Visually identify junk clusters, which are characterized by low FWS.Average, high FWS.Fill.factor, and low fluorescence. High fluorescence on one channel but very low fluorescence on other channels can indicate one type of junk, electronic noise. Typically the first and second largest clusters (black and red colours) are junk, and some much smaller clusters as well. These are specific to each dataset and must be decided based on the plots. Once they are identified, use the cluster numbers to clean the dataset. In our case, clusters 1,2,4,5,6,8,9 are junk.

#Re-examine plot without the putative junk clusters (this may reveal smaller ones that were not
#visible earlier)
junk.clusterids<-c(2,3)
dat.clustlog2<-data.frame(dat.clustlog,subflow$peaks.cluster)
dat.good<- subset(dat.clustlog2, !(dat.clustlog2$subflow.peaks.cluster %in% junk.clusterids))

# Examine plots of the cleaned data to identify more junk clusters
# pairs(dat.good,col=alpha(dat.good$subflow.peaks.cluster, 0.6), pch='.',lower.panel=NULL,gap=0.2)
# plot(dat.good$FL.Red.First~dat.good$FL.Red.Last,col=alpha(subflow$peaks.cluster+2, 0.4))
# scatter3Drgl(dat.good$FL.Red.First,dat.good$FL.Red.Last,dat.good$FL.Red.Range,colvar = dat.good$subflow.peaks.cluster, col=palette1[1:max(subflow$peaks.cluster)], xlab='FL.Red.First', ylab='FL.Red.Last', zlab='FL.Red.Range')
scatter3Drgl(dat.good$fl1_a,dat.good$fl2_a,dat.good$fl3_a, colvar = dat.good$subflow.peaks.cluster, col=palette1[1:max(dat.good$subflow.peaks.cluster)],xlab='fl1', ylab='fl2', zlab='fl3',ticktype="detailed")

scatter3Drgl(dat.good$FL.Red.Range,dat.good$X2.FL.Red.Range, dat.good$X2.FL.Red.Gradient,colvar = dat.good$subflow.peaks.cluster, col=palette1[1:max(subflow$peaks.cluster)],xlab='FL.Red.Range', ylab='X2.FL.Red.Range', zlab='X2.FL.Red.Gradient')

scatter3Drgl(dat.good$FL.Red.Range,dat.good$X2.FL.Red.Range, dat.good$FL.Red.Gradient,colvar = dat.good$subflow.peaks.cluster, col=palette1[1:max(subflow$peaks.cluster)],xlab='FL.Red.Range', ylab='X2.FL.Red.Range', zlab='FL.Red.Gradient')

#Recluster cleaned data to see if additional junk clusters are now evident
# (uses a slightly different set of variables)
cleanflow<-flowPeaks(dat.good[,c(1:5)], tol=0.25,h0=0.05, h=2)

#examine details if interested
summary(cleanflow)

#Plot reclustered data, coloured by new cluster identity
# pairs(dat.good[1:9],col=alpha(cleanflow$peaks.cluster, 0.6), pch='.',lower.panel=NULL,gap=0.2)

scatter3Drgl(dat.good$FWS.Range,dat.good$FL.Red.Range,dat.good$X2.FL.Red.Range,colvar = cleanflow$peaks.cluster, col=palette1[1:max(subflow$peaks.cluster)], xlab='FWS.Range', ylab='FL.Red.Range', zlab='X2.FL.Red.Range')
## JB addition:
scatter3Drgl(dat.good$fl1_a,dat.good$fl2_a,dat.good$fl3_a, colvar = dat.good$subflow.peaks.cluster, col=palette1[1:max(dat.good$subflow.peaks.cluster)],xlab='fl1', ylab='fl2', zlab='fl3',ticktype="detailed")

#Specify any new junk clusters that emerge; leave blank if there are none
junk.clusterids.round2<-c()

#Add cluster IDs to dataset
dat.good3<-data.frame(dat.good,cleanflow$peaks.cluster)

#Train a random forest classifier to identify the clusters in both cleaning rounds
#First round of cleaning
rf <- randomForest(factor(subflow.peaks.cluster) ~ fl1_a + fl2_a + fl3_a +
				   fsc_a + ssc_a, data=dat.clustlog2,importance=TRUE, ntree=1001, na.action=na.omit)

#Second round of cleaning
rf2 <- randomForest(factor(cleanflow.peaks.cluster) ~ fl1_a + fl2_a + fl3_a +
						fsc_a + ssc_a,
                    data=dat.clustlog3,importance=TRUE, ntree=1001, na.action=na.omit)

#Load the list of raw data files
filenames <-list.files(indir,pattern=c('Allparameters_'),full.names=FALSE,recursive=FALSE)



dat.clust2<-dat[,c('fsc_a','ssc_a','fl1_a', 'fl2_a',
				 'fl3_a')]

#Subset data to remove entirely ridiculous cases that will mess with clustering.
dat.clust2<-subset(dat.clust2,dat.clust2$ssc_a > 0)
dat.clust2<-subset(dat.clust2,dat.clust2$fl1_a > 0)
dat.clust2<-subset(dat.clust2,dat.clust2$fl2_a > 0)
dat.clust2<-subset(dat.clust2,dat.clust2$fsc_a > 0)
dat.clust2<-subset(dat.clust2,dat.clust2$fl3_a > 0)
dat.clust2 <- dat.clust2 %>% 
	rownames_to_column(var = "id")


classifydat<-function(x){
  #Read in data file, exclude unhelpful columns unless needed
  # dat<-fread(paste0(indir,x))
  # colnames(dat)[which(colnames(dat) == "Particle.ID")] = "id"
  # dat = dat[,FL.Red.Gradient := abs(FL.Red.First-FL.Red.Last)+0.1]
  # dat = dat[,X2.FL.Red.Gradient := abs(X2.FL.Red.First-X2.FL.Red.Last)+0.1]
  # 
  # #Impose same conditions as at the beginning of this process, or clustering may go awry.
  # dat = subset(dat,FL.Red.Range > 0)
  # dat = subset(dat,X2.FL.Red.Range > 0)
  # dat = subset(dat,FL.Orange.Range > 0)
  # dat = subset(dat,FL.Yellow.Range > 0)

  # dat = as.data.table(dat.clust2)
  # dat = subset(dat,SWS.Length > 1)

  #Arrange by id. Done because the data.table package requires this in order to perform later
  #functions
  setkey(dat,id)

  #Select only the variables used for clustering
  dat2=dat[,list(log10(fl1_a), log10(fl2_a), log10(fl3_a),
                 log10(fsc_a), log10(ssc_a))]

  dat.clustlog2<-data.frame(dat[, 1], log10(dat[,c(2:6)]))
  dat2 <- as.data.table(dat.clustlog2)
  
  # Rename columns
  # setnames(dat2, 1:7, new=c('fl1_a', 'fl2_a', 'fl3_a', 'fcs_a',
  #                           'ssc_a'))

  #Predict the cluster identity of particles using the random forests classifier (first round)
  rf.class<-predict(rf, dat2)
  # dat.clustlog3 <- as.data.table(dat.clustlog2)
  #Add a new column with the predicted cluster identity
  newdat <- dat2[,classification1:=as.numeric(rf.class)]

  #Order by cluster. May not be necessary.
  setkey(newdat, classification1)

  #Remove the clusters that are believed to be junk
  newdat2 = newdat[!(classification1 %in% junk.clusterids)]

  # Subset data for second round of clustering (in case additional junk clusters have been detected)
  dat3 = newdat2[,list(log10(FWS.Range), log10(FL.Red.Range), log10(X2.FL.Red.Range),
                       log10(FL.Orange.Range), log10(X2.FL.Red.Gradient),log10(FL.Red.Gradient),
                       FWS.Fill.factor)]

  setnames(dat3, 1:7, new=c('FWS.Range', 'FL.Red.Range', 'X2.FL.Red.Range', 'FL.Orange.Range',
                            'X2.FL.Red.Gradient','FL.Red.Gradient', 'FWS.Fill.factor'))
  #Predict the cluster identity of particles using the random forests classifier (second round)
  rf2.class<-predict(rf2, dat3)

  #Add a new column with the predicted cluster identity
  newdat3 = newdat2[,classification2:=as.numeric(rf2.class)]

  #Remove the clusters that are believed to be junk
  newdat3.good = newdat3[!(classification2 %in% junk.clusterids.round2)]

  # Remove classification columns if desired.
  #newdatfinal=newdat3.good[,c('classification1','classification2'):=NULL]
  # else
  newdatfinal=newdat3.good
  
  ##########Developing reduced dataset for student use - will only output cluster assignments and calculated biovolume

  
  
  #Order by ID. May not be necessary.
  setkey(newdatfinal, id)

  #Write cleaned file to .csv, file name edited by removing 'listmode_' and appending '_cleaned'
  #You might need to alter the digits within the 'substr' commands depending on how your files
  #are named

  student_out<-newdatfinal
  student_out$biovol<-(student_out$X2.FL.Red.Average^(4/3))/(10^(127/75))
  student_out<-student_out[,c("biovol","classification1","classification2")]
  
  #Cleaned files are written to same subfolder as original file
  write.csv(newdatfinal, paste0(outdir,substr(x,1,nchar(x)-4),'_cleaned','.csv'),row.names = FALSE,quote = FALSE)
  write.csv(student_out, paste0(outdir,substr(x,1,nchar(x)-4),'_student','.csv'),row.names = FALSE,quote = FALSE)
}

#Run function on all files. Can be run on multiple processors except on Windows
mapply(classifydat, x=filenames)
