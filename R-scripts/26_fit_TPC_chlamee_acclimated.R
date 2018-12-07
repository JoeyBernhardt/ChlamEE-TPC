
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


nbcurve<-function(temp,z,w,a,b){
	res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
	res
}


rfu <- read_csv("data-processed/chlamee-exponential.csv")

# rfu <- exponential 
rfu2 <- rfu %>% 
	rename(temp = temperature) %>% 
	filter(!is.na(RFU)) 

population_key <- read_excel("data-general/ChlamEE_Treatments_JB.xlsx") %>% 
	clean_names()

rfu2 <- left_join(rfu2, population_key, by = c("population"))

unique(rfu2$temp)

df <- rfu2 %>% 
	filter(population == 3)

df %>% 
	ggplot(aes(x = days, y = RFU, color = factor(temp))) + geom_point()


fit_growth <- function(df){
	res <- try(nlsLM(RFU ~  mean(c(df$RFU[df$days < 0.2])) * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
					 data= df,  
					 start= list(z= 25,w= 25,a= 0.2, b= 0.07),
					 lower = c(z = 0, w= 0, a = 0, b = 0),
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

##using growthTools fit
output_split <- read_csv("data-processed/chlamee-acute-tpc-fits.csv") %>% 
	rename(z = s) %>% 
	split(.$population)
topts <- output_split %>% 
	map_df(get_topt, .id = "population") 


topts %>% 
	mutate(population = as.factor(population)) %>% 
	ggplot(aes(x = fct_reorder(population, topt), y = topt)) + geom_point() +
	xlab("Population") +ylab("Topt")

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
limits_all2 <- left_join(limits_all, population_key, by = "population")

limits_all2 %>% 
	mutate(population = as.factor(population)) %>% 
	ggplot(aes(x = fct_reorder(population, tmax), y = tmax, color = treatment)) + geom_point() +
	xlab("Population") +ylab("Tmax")

limits_all2 %>% 
	filter(!is.na(treatment)) %>% 
	ggplot(aes(x = treatment, y = topt)) + geom_boxplot() +
	ylab("Topt (°C)") + xlab("Selection treatment")
limits_all2 %>% 
	filter(!is.na(treatment)) %>% 
	ggplot(aes(x = treatment, y = rmax)) + geom_boxplot() +
	ylab("rmax") + xlab("Selection treatment")

limits_all2 %>% 
	filter(!is.na(treatment)) %>%
	ggplot(aes(x = treatment, y = rmax)) + geom_boxplot() +
	ylab("rmax") + xlab("Selection treatment")

limits_all2 %>% 
	filter(!is.na(treatment)) %>% 
	ggplot(aes(x = treatment, y = tmax)) + geom_boxplot() +
	ylab("Tmax (°C)") + xlab("Selection treatment")

limits_all2 %>% 
	filter(!is.na(treatment)) %>% 
	ggplot(aes(x = treatment, y = rmax, color = ancestor_id)) + geom_boxplot() +
	ylab("rmax") + xlab("Selection treatment") 

limits_all2 %>% 
	filter(!is.na(treatment)) %>% 
	ggplot(aes(x = ancestor_id, y = rmax, color = ancestor_id)) + geom_boxplot() +
	ylab("rmax") + xlab("Selection treatment") 

limits_all2 %>% 
	filter(tmax > 38) %>% 
	filter(!is.na(treatment)) %>% 
	lm(rmax ~ treatment, data = .) %>% summary()

limits_all2 %>% 
	filter(tmax > 38) %>% 
	filter(!is.na(treatment)) %>% 
	lm(rmax ~ ancestor_id, data = .) %>% summary()


limits_all2 %>% 
	# filter(population != 7) %>% 
	mutate(population = as.factor(population)) %>% 
	ggplot(aes(x = fct_reorder(population, tmin), y = tmin)) + geom_point() +
	xlab("Population") +ylab("Tmin")



limits_all %>% 
	# filter(population != 7) %>% 
	mutate(population = as.factor(population)) %>% 
	ggplot(aes(x = topt, y = rmax)) + geom_point() +
	xlab("Topt") + ylab("rmax") +
	geom_smooth(method = "lm")

limits_all %>% 
	filter(tmax > 38) %>% 
	mutate(population = as.factor(population)) %>% 
	ggplot(aes(x = tmax, y = rmax)) + geom_point() +
	xlab("Tmax") + ylab("rmax") +
	geom_smooth(method = "lm")

limits_all %>% 
	# filter(population != 7) %>% 
	mutate(population = as.factor(population)) %>% 
	ggplot(aes(x = topt)) + geom_histogram()


all_traits <- left_join(limits_all, output, by = "population")

all_traits %>% 
	filter(w < 60) %>% 
	ggplot(aes(x = rmax, y = w)) + geom_point() +
	geom_smooth(method = "lm", color = "purple") +ylab("Thermal breadth") + xlab("Maximum growth rate")


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

df <- output_split[[1]]
prediction_nb <- function(df) {
	tpc <-function(x){
		res<- nbcurve2(x, df$topt[[1]], df$w[[1]], df$a[[1]], df$b[[1]])
		res
	}
	
	pred <- function(x) {
		y <- tpc(x)
	}
	
	x <- seq(0, 50, by = 0.1)
	
	preds <- sapply(x, tpc)
	preds <- data.frame(x, preds) %>% 
		rename(temperature = x, 
			   growth = preds)
}



# bs_split <- output %>% 
	# split(.$population)

all_preds <- output_split %>% 
	map_df(prediction_nb, .id = "population")


# all_preds <- bs_split %>% 
	# map_df(prediction_function, .id = "population")


all_preds2 <- left_join(all_preds, population_key, by = "population") %>% 
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>% 
	filter(ancestor_id != "cc1629")

all_preds2 %>% 
	dplyr::filter(!is.na(population)) %>% 
	dplyr::filter(!is.na(treatment)) %>% 
	ggplot(aes(x = temperature, y = growth, color = treatment, group = population)) + geom_line(size = 1) +
	ylim(0, 4.3) + xlim(0, 50) + geom_hline(yintercept = 0) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") + 
	scale_color_discrete(name = "Treatment")  +
	facet_wrap( ~ ancestor_id)

all_preds2 %>% 
	dplyr::filter(!is.na(population)) %>% 
	dplyr::filter(!is.na(treatment)) %>% 
	ggplot(aes(x = temperature, y = growth, color = ancestor_id, group = population)) + geom_line(size = 1) +
	ylim(0, 4.3) + xlim(0, 50) + geom_hline(yintercept = 0) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") + 
	scale_color_discrete(name = "Ancestor")  +
	facet_wrap( ~ treatment)

rmax <- all_preds2 %>% 
	group_by(population, treatment, ancestor_id) %>% 
	summarise(rmax = max(growth)) %>% 
	ungroup() %>% 
	mutate(treatment = ifelse(is.na(treatment), "none", treatment))

controls <- rmax %>% 
	filter(ancestor_id %in% c("anc4", "anc5", "cc1690"))

mod <- controls %>% 
	mutate(evolved = ifelse(treatment %in% c("none", "C"), "no", "yes")) %>% 
	filter(treatment %in% c("BS", "none", "C")) %>% 
	lm(rmax ~ treatment, data = .)

summary(mod)

library(nlme)
growth_mod <-lme(rmax ~ treatment, random= ~ treatment| ancestor_id, data= rmax)
summary(growth_mod)
anova(growth_mod)
intervals(growth_mod)


ancestors <- all_preds2 %>% 
	filter(treatment %in% c("none"))
all_preds2 %>% 
	dplyr::filter(!is.na(population)) %>% 
	dplyr::filter(!is.na(treatment)) %>% 
	ggplot(aes(x = temperature, y = growth, color = treatment, group = population)) + geom_line(size = 1) +
	ylim(0, 4.3) + xlim(0, 50) + geom_hline(yintercept = 0) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") + 
	geom_line(aes(x = temperature, y = growth), data = ancestors, color = "black", size = 1.2)+
	scale_color_discrete(name = "Treatment") + facet_grid( ~ ancestor_id) 
ggsave("figures/ancestors_chlamee.pdf", width = 14, height = 3)


anc4_controls <- all_preds2 %>% 
	filter(ancestor_id == "anc4", treatment %in% c("none", "C"))
all_preds2 %>% 
	dplyr::filter(ancestor_id == "anc4") %>% 
	ggplot(aes(x = temperature, y = growth, color = treatment, group = population)) + geom_line(size = 1) +
	ylim(0, 4.3) + xlim(0, 50) + geom_hline(yintercept = 0) +
	geom_line(aes(x = temperature, y = growth), data = anc4_controls, color = "black", size = 2)+
	ylab("Exponential growth rate") + xlab("Temperature (°C)") + 
	scale_color_discrete(name = "Ancestor") 
ggsave("figures/ancestor4_chlamee.pdf", width = 8, height = 5)

anc5_controls <- all_preds2 %>% 
	filter(ancestor_id == "anc5", treatment %in% c("none", "C"))
all_preds2 %>% 
	dplyr::filter(ancestor_id == "anc5") %>% 
	ggplot(aes(x = temperature, y = growth, color = treatment, group = population)) + geom_line(size = 1) +
	ylim(0, 4.3) + xlim(0, 50) + geom_hline(yintercept = 0) +
	geom_line(aes(x = temperature, y = growth), data = anc5_controls, color = "black", size = 2)+
	ylab("Exponential growth rate") + xlab("Temperature (°C)") + 
	scale_color_discrete(name = "Ancestor") 
ggsave("figures/ancestor5_chlamee.pdf", width = 8, height = 5)

all_preds2 %>% 
	dplyr::filter(treatment == "none") %>% 
	ggplot(aes(x = temperature, y = growth, color = ancestor_id, group = population)) + geom_line(size = 1) +
	ylim(0, 4.3) + xlim(0, 50) + geom_hline(yintercept = 0) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") + 
	scale_color_discrete(name = "Ancestor") 
ggsave("figures/ancestors_only_chlamee.pdf", width = 8, height = 5)
