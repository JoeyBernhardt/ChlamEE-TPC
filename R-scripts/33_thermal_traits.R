
library(minpack.lm)
library(broom)
library(nlstools)
library(nls.multstart)
library(tidyverse)
library(cowplot)
library(rootSolve)
library(forcats)
library(growthTools)
library(readxl)
library(janitor)



rfu <- read_csv("data-processed/chlamee-exponential.csv")

rfu2 <- rfu %>% 
	rename(temp = temperature) %>% 
	filter(!is.na(RFU)) %>% 
	filter(population != "cc1629")

population_key <- read_excel("data-general/ChlamEE_Treatments_JB.xlsx") %>% 
	clean_names() %>% 
	mutate(treatment = ifelse(is.na(treatment), "ancestors", treatment)) %>% 
	filter(population != "cc1629")

rfu2 <- left_join(rfu2, population_key, by = c("population"))


fit_growth <- function(df){
	res <- try(nlsLM(RFU ~  mean(c(df$RFU[df$days < 0.2])) * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
					 data= df,  
					 start= list(z= 25,w= 25,a= 0.2, b= 0.07),
					 lower = c(z = 0, w= 0, a = 0, b = 0),
					 upper = c(z = 40, w= 60,a =  2, b = 2),
					 control = nls.control(maxiter=1024, minFactor=1/204800000)))
	if(class(res)!="try-error"){
		out1 <- tidy(res) %>%  
			select(estimate, term) %>%  
			spread(key = term, value = estimate)
		out2 <- glance(res)
	}
	all <- bind_cols(out1, out2)
	all
}


df_split <- rfu2 %>% 
	filter(!is.na(RFU)) %>% 
	split(.$population)


output <- df_split %>%
	map_df(fit_growth, .id = "population") 

og <- output %>% 
	mutate(growth_10 = a*exp(b*10)*(1-((10-z)/(w/2))^2)) %>% 
	mutate(growth_15 = a*exp(b*15)*(1-((15-z)/(w/2))^2)) %>% 
	mutate(growth_20 = a*exp(b*20)*(1-((20-z)/(w/2))^2)) %>% 
	mutate(growth_25 = a*exp(b*25)*(1-((25-z)/(w/2))^2)) %>% 
	mutate(growth_30 = a*exp(b*30)*(1-((30-z)/(w/2))^2)) %>% 
	mutate(growth_35 = a*exp(b*35)*(1-((35-z)/(w/2))^2)) %>% 
	mutate(growth_40 = a*exp(b*40)*(1-((40-z)/(w/2))^2))



# merge with P params -----------------------------------------------------


phosphate_params <- read.csv("data-processed/phosphate-params-monod-bootstrapped.csv") %>% 
	group_by(population) %>% 
	summarise_each(funs = mean, umax, ks)

m <- 0.56
all_traits <- left_join(og, phosphate_params) %>% 
	left_join(population_key) %>% 
	filter(population != "cc1629") %>% 
	mutate(rstar = ks*m/(umax-m)) %>% 
	mutate(growth_50P = umax*(50/(ks + 50)))

all_traits %>% 
	ggplot(aes(x = rstar, y = growth_30)) + geom_point(size = 1.5) +
	ylab("Growth at 30°C (/day)") + xlab("P* (uM P)") 
ggsave("figures/growth-30-pstar.png", width = 8, height = 6)

all_traits %>% 
	gather(contains("growth"), key = temperature, value = growth_rate) %>% 
	mutate(temperature = str_replace(temperature, "growth_", "")) %>% 
	mutate(temperature = as.factor(temperature)) %>% 
	ggplot(aes(x = rstar, y = growth_rate, color = temperature, fill = temperature)) + geom_point() +
	geom_smooth(method = "lm", aes(color = temperature, fill = temperature)) + ylim(0, 4) +
	ylab("Temperature dependent growth rate") + xlab("P*") +
	scale_color_viridis_d() + scale_fill_viridis_d()
ggsave("figures/growth-temps-pstar.png", width = 8, height = 6)


all_traits %>% 
	ggplot(aes(x = umax, y = growth_20)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("TPC growth rate at 20°C") + xlab("umax from P experiments") +
	ylim(0.6, 1.8) + xlim(0.6, 1.8)
ggsave("figures/tpc-growth-umax-P.png", width = 8, height = 6)

all_traits %>% 
	ggplot(aes(x = growth_50P, y = growth_20, color = treatment)) + geom_point() +
	geom_abline(intercept = 0, slope = 1) +
	ylab("TPC growth rate at 20°C") + xlab("P expt growth 50 uM P") +
	ylim(0.6, 1.8) + xlim(0.6, 1.8)
ggsave("figures/tpc-growth-P.png", width = 8, height = 6)




# plot the fits -----------------------------------------------------------

prediction_function <- function(df) {
	tpc <-function(x){
		res<-(df$a[[1]]*exp(df$b[[1]]*x)*(1-((x-df$z[[1]])/(df$w[[1]]/2))^2))
		res
	}
	
	pred <- function(x) {
		y <- tpc(x)
	}
	
	x <- seq(0, 50, by = 0.1)
	
	preds <- sapply(x, pred)
	preds <- data.frame(x, preds) %>% 
		rename(temperature = x, 
			   growth = preds)
}


bs_split <- output %>% 
	split(.$population)


all_preds <- bs_split %>% 
	map_df(prediction_function, .id = "population") %>% 
	left_join(population_key)


all_preds %>% 
	ggplot(aes(x = temperature, y = growth, color = factor(treatment))) + geom_line(size = 1) +
	ylim(0, 4.3) + xlim(0, 50) + geom_hline(yintercept = 0) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") + 
	scale_color_discrete(name = "Treatment") + facet_wrap( ~ population)





get_topt <- function(df){
	grfunc<-function(x){
		-nbcurve(x, z = df$z[[1]],w = df$w[[1]],a = df$a[[1]],b = df$b[[1]])
	}
	optinfo<-optim(c(x=df$z[[1]]),grfunc)
	opt <-c(optinfo$par[[1]])
	maxgrowth <- c(-optinfo$value)
	results <- data.frame(topt = opt, rmax = maxgrowth)
	return(results)
}

output_split <- output %>% 
	split(.$population)



tpc_params <- read_csv("data-processed/chlamee-acute-tpc-fits.csv") %>% 
	filter(population != "cc1629")

phosphate_params <- read.csv("data-processed/phosphate-params-monod-bootstrapped.csv") %>% 
	group_by(population) %>% 
	summarise_each(funs = mean, umax, ks)

all_traits <- left_join(tpc_params, phosphate_params)

all_traits %>% 
	ggplot(aes(x = tmax, y = ks, color = treatment)) + geom_point() 


### ok now let's make the umax R* trade-off plot at different temperatures

growth <- a*exp(b*temp)*(1-((temp-z)/(w/2))^2)


### ok this makes no sense...the TPC parameters don't seem to correct
all_traits %>% 
	rename(z = s) %>% 
	mutate(growth_10 = a*exp(b*10)*(1-((10-z)/(w/2))^2)) %>% 
	mutate(growth_20 = a*exp(b*20)*(1-((20-z)/(w/2))^2)) %>% View

### let's try refitting
rfu <- read_csv("data-processed/chlamee-exponential.csv")

rfu2 <- rfu %>% 
	rename(temp = temperature) %>% 
	filter(!is.na(RFU)) 


rfu2 %>% 
	ggplot(aes(x = days, y = RFU, group = population, color = temp)) + geom_point() +
	facet_wrap( ~ population)

library(nls.multstart)

output2 <- rfu2 %>% 
	group_by(population) %>%
	mutate(N0_mean = mean(N0)) %>% 
	filter(population == 1) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(RFU ~ N0_mean * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
												  data = .x,
												  iter = 500,
												  start_lower = c(z= 25,w= 25,a= 0.2, b= 0.1),
												  start_upper = c(z= 35,w= 35,a= 0.4, b= 0.2),
												  supp_errors = 'N',
												  na.action = na.omit,
												  lower = c(z = 0, w= 0, a = -0.2, b = 0),
												  upper = c(z = 40, w= 80,a =  2, b = 2),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))

info <- output2 %>%
	unnest(fit %>% map(glance))
params <- output2 %>%
	unnest(fit %>% map(tidy))

fit_growth <- function(df){
	res <- try(nlsLM(RFU ~  mean(c(df$RFU[df$days < 0.25])) * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
					 data= df,  
					 start= list(z= 25,w= 25,a= 0.2, b= 0.1),
					 lower = c(z = 0, w= 0, a = -0.2, b = 0),
					 upper = c(z = 40, w= 80,a =  2, b = 2),
					 control = nls.control(maxiter=1024, minFactor=1/204800000)))
	if(class(res)!="try-error"){
		out1 <- tidy(res) %>%  
			select(estimate, term) %>%  
			spread(key = term, value = estimate)
		out2 <- glance(res)
	}
	all <- bind_cols(out1, out2)
	all
}


df_split <- rfu2 %>% 
	filter(!is.na(RFU)) %>% 
	split(.$population)


output <- df_split %>%
	map_df(fit_growth, .id = "population")


output <- df_split %>%
	map_df(fit_growth, .id = "population") 

output_growth <- output %>% 
	mutate(growth_10 = a*exp(b*10)*(1-((10-z)/(w/2))^2)) %>% 
	mutate(growth_20 = a*exp(b*10)*(1-((10-z)/(w/2))^2)) %>% View


get_topt <- function(df){
	grfunc<-function(x){
		-nbcurve(x, z = df$z[[1]],w = df$w[[1]],a = df$a[[1]],b = df$b[[1]])
	}
	optinfo<-optim(c(x=df$z[[1]]),grfunc)
	opt <-c(optinfo$par[[1]])
	maxgrowth <- c(-optinfo$value)
	results <- data.frame(topt = opt, rmax = maxgrowth)
	return(results)
}

output_split <- output %>% 
	split(.$population)

library(rootSolve)
nbcurve<-function(temp,z,w,a,b){
	res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
	res
}

topts <- output_split %>% 
	map_df(get_topt, .id = "population") 

get_tmax <- function(df){
	uniroot.all(function(x) nbcurve(x, z = df$z[[1]],w = df$w[[1]],a = df$a[[1]], b = df$b[[1]]),c(30,150))
}

tmaxes <- output_split %>% 
	map(get_tmax) %>% 
	unlist() %>% 
	as_data_frame() %>% 
	rename(tmax = value) %>% 
	mutate(population = rownames(.)) 

get_tmin <- function(df){
	uniroot.all(function(x) nbcurve(x, z = df$z[[1]],w = df$w[[1]],a = df$a[[1]], b = df$b[[1]]),c(-40,25))
}

tmins <- output_split %>% 
	map(get_tmin) %>% 
	unlist() %>% 
	as_data_frame() %>% 
	rename(tmin = value) %>% 
	mutate(population = rownames(.))


limits <- left_join(tmins, tmaxes, by = "population")
limits_all <- left_join(limits, topts, by = "population")
