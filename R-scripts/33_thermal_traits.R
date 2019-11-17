
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
	filter(!is.na(RFU)) 

population_key <- read_excel("data-general/ChlamEE_Treatments_JB.xlsx") %>% 
	clean_names()

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



tpc_params <- read_csv("data-processed/chlamee-acute-tpc-fits.csv")
