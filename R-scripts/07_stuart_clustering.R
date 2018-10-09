#Install and load required packages

#Load packages

library(flowPeaks); library(scales);library(data.table);library(randomForest);library(plot3Drgl)
library(cowplot)
library(tidyverse)
#Set seed for repeatability
set.seed(42)

#Set your working directory, which contains both the data subset and the files to be cleaned (the latter can be in subfolders)

#Read in the data subset to be used for clustering
# df <- read.csv('data/alldatasubset_100k_2018_Limno.csv')
df_raw <- read_csv("data-processed/anc4-particles.csv") 
df <- df_raw %>% 
	dplyr::filter(ssc_a > 0, fsc_a > 0, fl1_a > 0, fl2_a >0, fl3_a>0) %>% 
	select(plate, well, ssc_a, fsc_a, fl1_a, fl2_a, fl3_a) %>% 
	rownames_to_column(var = "id") %>% 
	mutate_at(vars(contains("_a")), ~ (log10(.)))


dat <- df_raw %>% 
	sample_n(size = 100000, replace = FALSE)




## pull out variables of interest
dat.clust<-df[,c('fsc_a','ssc_a','fl1_a', 'fl2_a',
				 'fl3_a')]

#Cluster
subflow <- flowPeaks(dat.clust[,c(1:5)], tol=0.25, h0 = 0.05, h=2)


# examine details if interested. The 'weight' column indicatest the proportion of the data belonging
# to that cluster
summary(subflow)

#Define colour palette for clusters in 3D plots
palette1 <- c('blue','cyan','magenta','yellow','gray','brown','purple','pink','orange','darkgreen')
dat.clust2<-data.frame(dat.clust,subflow$peaks.cluster)

df_class <- left_join(df, dat.clust2, by = c("fl1_a", "fl2_a", "fl3_a", "ssc_a", "fsc_a"))

### here we want cluster class 1, this is the algae

write_csv(df_class, "data-processed/anc4-classified.csv")

cell_counts <- df_class %>% 
group_by(plate, well, subflow.peaks.cluster) %>% 
	tally() %>% 
	dplyr::filter(subflow.peaks.cluster == "1")

write_csv(cell_counts, "data-processed/cell_counts_flow_peaks.csv")

#Examine any plots of interest
scatter3Drgl(dat.clust$fl1_a,dat.clust$fl2_a,dat.clust$fl3_a,colvar = subflow$peaks.cluster, col=palette1[1:max(subflow$peaks.cluster)],xlab='fl1', ylab='fl2', zlab='fl3',ticktype="detailed")


plot(dat.clustlog$fl1_a,dat.clustlog$fl2_a)

#Re-examine plot without the putative junk clusters (this may reveal smaller ones that were not
#visible earlier)
junk.clusterids<-c(2,3, 4)
dat.clustlog2<-data.frame(dat.clustlog,subflow$peaks.cluster)
dat.good<- subset(dat.clustlog2, !(dat.clustlog2$subflow.peaks.cluster %in% junk.clusterids))

# Examine plots of the cleaned data to identify more junk clusters
scatter3Drgl(dat.good$fl1_a,dat.good$fl2_a,dat.good$fl3_a, colvar = dat.good$subflow.peaks.cluster, col=palette1[1:max(dat.good$subflow.peaks.cluster)],xlab='fl1', ylab='fl2', zlab='fl3',ticktype="detailed")

#Recluster cleaned data to see if additional junk clusters are now evident
# (uses a slightly different set of variables)
cleanflow<-flowPeaks(dat.good[,c(1:5)], tol=0.25,h0=0.05, h=2)

#examine details if interested
summary(cleanflow)

#Plot reclustered data, coloured by new cluster identity
## JB addition:
scatter3Drgl(dat.good$fl1_a,dat.good$fl2_a,dat.good$fl3_a,colvar = cleanflow$peaks.cluster, col=palette1[1:max(subflow$peaks.cluster)], xlab='fl1', ylab='fl2', zlab='fl3')

#Specify any new junk clusters that emerge; leave blank if there are none
junk.clusterids.round2<-c(2,3)

#Add cluster IDs to dataset
dat.clustlog3 <- data.frame(dat.good,cleanflow$peaks.cluster)


# Random forest -----------------------------------------------------------

#Train a random forest classifier to identify the clusters in both cleaning rounds
#First round of cleaning
rf <- randomForest(factor(subflow.peaks.cluster) ~ fl1_a + fl2_a + fl3_a +
				   fsc_a + ssc_a, data = dat.clustlog2, importance=TRUE, ntree=1001, na.action=na.omit)

#Second round of cleaning
rf2 <- randomForest(factor(cleanflow.peaks.cluster) ~ fl1_a + fl2_a + fl3_a +
						fsc_a + ssc_a,
                    data=dat.clustlog3,importance=TRUE, ntree=1001, na.action=na.omit)

#Load the list of raw data files

dat <- df_raw %>% 
	sample_n(size = 100000, replace = FALSE)

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


dat.clustlog2<-data.frame(dat.clust2[, 1], log10(dat.clust2[,c(2:6)]))
dat2 <- as.data.table(dat.clustlog2)
  #Arrange by id. Done because the data.table package requires this in order to perform later
  #functions
  setkey(dat2,id)

  #Predict the cluster identity of particles using the random forests classifier (first round)
  rf.class<-predict(rf, dat2)
  # dat.clustlog3 <- as.data.table(dat.clustlog2)
  #Add a new column with the predicted cluster identity
  newdat <- dat2[,classification1:=as.numeric(rf.class)]
  
newdat %>% 
	ggplot(aes(x = fl1_a, y = fl3_a, color = factor(classification1))) + geom_point()

  #Order by cluster. May not be necessary.
  setkey(newdat, classification1)

  #Remove the clusters that are believed to be junk
  newdat2 = newdat[!(classification1 %in% junk.clusterids)]
  
  
 

  # Subset data for second round of clustering (in case additional junk clusters have been detected)
  dat3 = newdat2[,list(log10(fl1_a), log10(fl2_a), log10(fl3_a),
                       log10(fsc_a), log10(ssc_a))]

  setnames(dat3, 1:5, new=c('fl1_a', 'fl2_a', 'fl3_a', 'fsc_a',
                            'ssc_a'))
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
