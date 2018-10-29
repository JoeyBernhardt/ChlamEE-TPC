

library(mleTools)
library(growthTools)
library(tidyverse)



all_rfus4 <- read_csv("data-processed/globe-chlamy-acclimated-RFU-time.csv")

pop14 <- all_rfus3 %>% 
	filter(population == 14) %>% 
	mutate(ln.fluor = log(RFU))

res <- get.growth.rate(pop14$days, log(pop14$RFU))
res
res$best.model
res$best.slope
res$best.model.rsqr
res$best.se
res$best.model.slope.n
res$best.model.slope.r2
summary(res$models$gr.lag)
res$slopes

gdat <- pop14 %>%
	group_by(temperature) %>%
	do(grs=get.growth.rate(x=.$days,y=.$ln.fluor,id=.$temperature,plot.best.Q=T,fpath=NA))

pop14_growth <- gdat %>%
	summarise(temperature,mu=grs$best.slope,
			  best.model=grs$best.model,best.se=grs$best.se)



all4 %>% 
	filter(days < 0.5) %>% 
	ggplot(aes(x = days, y = RFU, color = factor(population), group = well_plate)) +
	geom_point(size = 2) +
	# scale_color_viridis_d(name = "Temperature") +
	xlab("Days") +
	facet_wrap( ~ temperature, scales = "free") 


all_rfus4 %>% 
	filter(population == 6, temperature == 40) %>% 
	ggplot(aes(x = days, y = RFU, color = factor(temperature), group = well_plate)) +
	geom_point(size = 2) +
	# scale_color_viridis_d(name = "Temperature") +
	xlab("Days") +
	facet_wrap( ~ population, scales = "free") +
	geom_line() 

pop8 <- all_rfus4 %>% 
	filter(population == 6, temperature != 22) %>% 
	mutate(ln.fluor = log(RFU))

gdat8 <- pop8 %>%
	group_by(temperature) %>%
	do(grs=get.growth.rate(x=.$days,y=.$ln.fluor,id=.$temperature,plot.best.Q=T,fpath=NA))

pop8_growth <- gdat8 %>%
	summarise(temperature,mu=grs$best.slope,
				  best.model=grs$best.model,best.se=grs$best.se)


nbcurve.traits8 <- get.nbcurve.tpc(pop8_growth$temperature,pop8_growth$mu,method='grid.mle2',
								 plotQ=T,conf.bandQ =T,fpath=NA)

