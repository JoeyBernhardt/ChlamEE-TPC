
library(minpack.lm)
library(broom)
library(nlstools)
library(nls.multstart)
library(tidyverse)
library(cowplot)

rfu <- read_csv("data-processed/globe-chlamy-exponential-RFU-time.csv")

rfu <- exponential 
rfu2 <- rfu %>% 
	dplyr::rename(temp = temperature) %>% 
	filter(!is.na(RFU)) 

rfu2 %>% 
	filter(population == 4) %>% 
	ggplot(aes(x = days, y = RFU, color = factor(temp))) + geom_point()


p4 <- rfu2 %>% 
	filter(population == 2)

p4 %>% 
	ggplot(aes(x = days, y = RFU, color = factor(temp))) + geom_point()

##cell_density ~ 800 * exp(r*days)

df <- p4
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

library(forcats)
library(rootSolve)

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

limits_all %>% 
	mutate(population = as.factor(population)) %>% 
	ggplot(aes(x = fct_reorder(population, tmax), y = tmax)) + geom_point() +
	xlab("Population") +ylab("Tmax")

limits_all %>% 
	filter(population != 7) %>% 
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
	# filter(population != 7) %>% 
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
	filter(population != 7) %>% 
	# filter(population %in% c(11, 3, 5)) %>% 
	ggplot(aes(x = rmax, y = w)) + geom_point() +
	geom_smooth(method = "lm", color = "purple") +ylab("Thermal breadth") + xlab("Maximum growth rate")
ggsave("figures/generalist-specialist.pdf", width = 6, height = 4)


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
		dplyr::rename(temperature = x, 
			   growth = preds)
}


bs_split <- output %>% 
	split(.$population)


all_preds <- bs_split %>% 
	map_df(prediction_function, .id = "population")

all_preds %>% 
	# filter(population %in% c(11, 3, 5)) %>% 
	mutate(population = as.integer(population)) %>% 
	ggplot(aes(x = temperature, y = growth, color = factor(population))) + geom_line(size = 1) +
	ylim(0, 3.5) + xlim(0, 50) + geom_hline(yintercept = 0) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") + scale_color_discrete(name = "Population")
ggsave("figures/globe-chlamy-TPCs.pdf", width = 8, height = 4)
ggsave("figures/globe-chlamy-TPCs-acclimated.pdf", width = 12, height = 6)

all_preds_acclimated <- all_preds
write_csv(all_preds_acclimated, "data-processed/all_preds_acclimated.csv")
write_csv(all_preds, "data-processed/all_preds_not_acclimated.csv")

all_preds %>% 
	# filter(population == 14) %>% 
	mutate(population = as.integer(population)) %>% 
	ggplot(aes(x = temperature, y = growth, color = factor(population))) + geom_line(size = 1) +
	ylim(0, 3.5) + xlim(0, 50) +
	facet_wrap( ~ population) + geom_hline(yintercept = 0) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") + scale_color_discrete(name = "Population")
ggsave("figures/globe-chlamy-TPCs-facet.pdf", width = 12, height = 6)

# bootstrapping -----------------------------------------------------------

df <- rfu2 %>% 
	filter(population == 4, !is.na(RFU))

df$RFU[[2]]

fit1 <- nlsLM(RFU ~ mean(c(df$RFU[df$days ==0])) * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
				 data= df,  
			start= list(z= 25,w= 25,a= 0.2, b= 0.07),
			  lower = c(z = 0, w= 0, a = -0.2, b = 0),
			  upper = c(z = 40, w= 80,a =  2, b = 2),
				 control = nls.control(maxiter=1024, minFactor=1/204800000))


fit_bootstrap <- function(df){
bootnls <- df %>% 
	bootstrap(100) %>% 
	do(tidy(nlsLM(RFU ~ mean(c(.$RFU[1:28]))  * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
	  data= .,  
	  start= list(z= 25,w= 25,a= 0.2, b= 0.07),
	  lower = c(z = 0, w= 0, a = -0.2, b = 0),
	  upper = c(z = 40, w= 80,a =  2, b = 2),
	  control = nls.control(maxiter=1024, minFactor=1/204800000))))
}


df_split <- rfu2 %>% 
	filter(!is.na(RFU)) %>% 
	split(.$population)
output_bs2 <- df_split %>%
	map_df(fit_bootstrap, .id = "population") 


df3 <- df_split[[8]]
fit3 <- nlsLM(RFU ~ mean(c(df3$RFU[days == 0]))  * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
	  data= df3,  
	  start= list(z= 25,w= 25,a= 0.2, b= 0.07),
	  lower = c(z = 0, w= 0, a = -0.2, b = 0),
	  upper = c(z = 40, w= 80,a =  2, b = 2),
	  control = nls.control(maxiter=1024, minFactor=1/204800000))

boot8 <- nlsBoot(fit3)



nls_boot_coefs_1 <- as_data_frame(boot1$coefboot) %>% 
	mutate(population = 1)
nls_boot_coefs_2 <- as_data_frame(boot2$coefboot) %>% 
	mutate(population = 2)
nls_boot_coefs_3 <- as_data_frame(boot3$coefboot) %>% 
	mutate(population = 3)
nls_boot_coefs_4 <- as_data_frame(boot4$coefboot) %>% 
	mutate(population = 4)
nls_boot_coefs_5 <- as_data_frame(boot5$coefboot) %>% 
	mutate(population = 5)
nls_boot_coefs_6 <- as_data_frame(boot6$coefboot) %>% 
	mutate(population = 6)
nls_boot_coefs_7 <- as_data_frame(boot7$coefboot) %>% 
	mutate(population = 7)
nls_boot_coefs_8 <- as_data_frame(boot8$coefboot) %>% 
	mutate(population = 8)
nls_boot_coefs_9 <- as_data_frame(boot9$coefboot) %>% 
	mutate(population = 9)
nls_boot_coefs_10 <- as_data_frame(boot10$coefboot) %>% 
	mutate(population = 10)
nls_boot_coefs_11 <- as_data_frame(boot11$coefboot) %>% 
	mutate(population = 11)
nls_boot_coefs_12 <- as_data_frame(boot12$coefboot) %>% 
	mutate(population = 12)
nls_boot_coefs_13 <- as_data_frame(boot13$coefboot) %>% 
	mutate(population = 13)
nls_boot_coefs_14 <- as_data_frame(boot14$coefboot) %>% 
	mutate(population = 14)
nls_boot_coefs_15 <- as_data_frame(boot15$coefboot) %>% 
	mutate(population = 15)

all_boots <- bind_rows(nls_boot_coefs_1,nls_boot_coefs_2, nls_boot_coefs_3, nls_boot_coefs_4, nls_boot_coefs_5,
					   nls_boot_coefs_6, nls_boot_coefs_7, nls_boot_coefs_8, nls_boot_coefs_9, nls_boot_coefs_10,
					   nls_boot_coefs_11, nls_boot_coefs_12, nls_boot_coefs_13, nls_boot_coefs_14, nls_boot_coefs_15) %>% 
	rownames_to_column()


all_boots_split <- all_boots %>% 
	split(.$rowname)

all_preds_bs <- all_boots_split %>% 
	map_df(prediction_function, .id = "unique_id")

all_preds_bs <- left_join(all_preds_bs, all_boots, by = c("unique_id" = "rowname"))

limits_c <- all_preds_bs %>% 
	group_by(temperature, population) %>% 
	summarise(q2.5=quantile(growth, probs=0.025),
			  q97.5=quantile(growth, probs=0.975),
			  mean = mean(growth)) 



nls_boot_c <- nlsBoot(fit3, niter = 1000)
nls_boot_coefs_c <- as_data_frame(boot11$coefboot)
ctpc <- as_data_frame(best_fit_c) %>% 
	rownames_to_column(.) %>% 
	spread(key = rowname, value = value)




fit_boot <- function(df){
	res <- try(nlsLM(RFU ~ mean(c(df$RFU[df$days ==0])) * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
					 data= df,  
					 start= list(z= 25,w= 25,a= 0.2, b= 0.1),
					 lower = c(z = 0, w= 0, a = -0.2, b = 0),
					 upper = c(z = 40, w= 80,a =  2, b = 2),
					 control = nls.control(maxiter=1024, minFactor=1/204800000)))
	if(class(res)!="try-error"){
		nls_boot_c <- nlsBoot(res, niter = 1000)
		nls_boot_coefs_c <- as_data_frame(nls_boot_c$coefboot)
	}
	all <- bind_cols(nls_boot_coefs_c)
	all
}



df_split <- rfu2 %>% 
	filter(!is.na(RFU)) %>% 
	split(.$population)

output_boot <- df_split %>%
	map_df(fit_boot, .id = "population") 





# ok now let’s take the bs replicates and predict them! -------------------


nls_boot_c <- nlsBoot(fit1, niter = 1000)
summary(fit1)
tidy(fit1, conf.int = TRUE)
coef(fit1)
confint2(fit1)
best_fit_c <- coef(fit1)
nls_boot_c <- nlsBoot(fit1, niter = 1000)
nls_boot_coefs_c <- as_data_frame(nls_boot_c$coefboot)
ctpc <- as_data_frame(best_fit_c) %>% 
	rownames_to_column(.) %>% 
	spread(key = rowname, value = value)

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


bs_split <- output_bs %>%
	unite(unique_id, replicate, population, sep = "_", remove = FALSE) %>% 
	ungroup() %>% 
	select(unique_id, term, estimate) %>% 
	spread(key = term, value = estimate) %>% 
	split(.$unique_id)

bs_split[[1]]

all_preds_bs <- bs_split %>% 
	map_df(prediction_function, .id = "unique_id")

limits_c <- all_preds_bs %>% 
	separate(unique_id, into = c("replicate", "population")) %>% 
	group_by(temperature, population) %>% 
	summarise(q2.5=quantile(growth, probs=0.025),
			  q97.5=quantile(growth, probs=0.975),
			  mean = mean(growth)) 

# limits_c3 <- all_preds_bs %>% 
# 	group_by(temperature) %>% 
# 	summarise(q2.5=quantile(growth, probs=0.025),
# 			  q97.5=quantile(growth, probs=0.975),
# 			  mean = mean(growth)) 

unique(limits_c$population)

limits_c <- limits_c %>% 
	mutate(population = as.factor(as.numeric(population)))
all_preds <- all_preds %>% 
	mutate(population = as.factor(as.numeric(population)))

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p +
# geom_line(aes(x = temperature, y = growth), data = filter(all_preds, population == 1)) + 
	# geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA),
	# 			data = limits_c, fill = "orange", alpha = 0.5) +
	geom_line(aes(x = temperature, y = growth, color = factor(population)), data = all_preds) +
	geom_line(aes(x = temperature, y = mean, color = factor(population)), data = limits_c) +
	geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, fill = factor(population)),
				data = limits_c, alpha = 0.5) +
	coord_cartesian(ylim = c(0, 3), xlim = c(0, 50)) +
	# facet_wrap( ~ population) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") + scale_color_discrete(name = "Population")



p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + 
	# stat_function(fun = tpc_c, color = "green", size = 1) +
	geom_line(aes(x = temperature, y = growth, group = replicate), color = "cadetblue", data = all_preds_c, alpha = 0.2) +
	stat_function(fun = ctpc, color = "black", size = 1) +
	stat_function(fun = vtpc, color = "orange", size = 1) +
	# # geom_line(aes(x = temperature, y = growth, group = replicate), color = "orange", data = all_preds_v, alpha = 0.1) +
	geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_v, fill = "orange", alpha = 0.5) +
	geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_c, fill = "cadetblue", alpha = 0.5) +
	





# other stuff -------------------------------------------------------------


fitc1 <- nlsLM(RFU ~ 20 * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
			  data= filter(rfu2, population == 9),  
			  start=list(z= results$z[[1]],w= results$w[[1]],a= results$a[[1]], b= results$b[[1]]),
			  lower = c(z = 0, w= 0, a = -0.2, b = 0),
			  upper = c(z = 30, w= 80,a =  0.5, b = 0.15),
			  control = nls.control(maxiter=1024, minFactor=1/204800000))

best_fit_c <- coef(fitc1)
ctpc <- as_data_frame(best_fit_c) %>% 
	rownames_to_column(.) %>% 
	spread(key = rowname, value = value)


fc9 <- ctpc
tpc_c9<-function(x){
	res<-(fc9$a[[1]]*exp(fc9$b[[1]]*x)*(1-((x-fc9$z[[1]])/(fc9$w[[1]]/2))^2))
	res
}
fc7 <- ctpc
tpc_c7<-function(x){
	res<-(fc7$a[[1]]*exp(fc7$b[[1]]*x)*(1-((x-fc7$z[[1]])/(fc7$w[[1]]/2))^2))
	res
}

fc2 <- ctpc
tpc_c2<-function(x){
	res<-(fc2$a[[1]]*exp(fc2$b[[1]]*x)*(1-((x-fc2$z[[1]])/(fc2$w[[1]]/2))^2))
	res
}

fc12 <- ctpc
tpc_c12<-function(x){
	res<-(fc12$a[[1]]*exp(fc12$b[[1]]*x)*(1-((x-fc12$z[[1]])/(fc12$w[[1]]/2))^2))
	res
}

fc8 <- ctpc
tpc_c8<-function(x){
	res<-(fc8$a[[1]]*exp(fc8$b[[1]]*x)*(1-((x-fc8$z[[1]])/(fc8$w[[1]]/2))^2))
	res
}


fc4 <- ctpc
tpc_c4<-function(x){
	res<-(fc4$a[[1]]*exp(fc4$b[[1]]*x)*(1-((x-fc4$z[[1]])/(fc4$w[[1]]/2))^2))
	res
}


fc10 <- ctpc
tpc_c10 <-function(x){
	res<-(fc10$a[[1]]*exp(fc10$b[[1]]*x)*(1-((x-fc10$z[[1]])/(fc10$w[[1]]/2))^2))
	res
}

fc6 <- ctpc
tpc_c6 <-function(x){
	res<-(fc6$a[[1]]*exp(fc6$b[[1]]*x)*(1-((x-fc6$z[[1]])/(fc6$w[[1]]/2))^2))
	res
}

fc1 <- ctpc
tpc_c <-function(x){
	res<-(fc1$a[[1]]*exp(fc1$b[[1]]*x)*(1-((x-fc1$z[[1]])/(fc1$w[[1]]/2))^2))
	res
}

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p + 
	stat_function(fun = tpc_c, color = "black", size = 1) +
	stat_function(fun = tpc_c2, color = "cadetblue", size = 1) +
	stat_function(fun = tpc_c7, color = "darkgreen", size = 1) +
	stat_function(fun = tpc_c9, color = "red", size = 1) +
	stat_function(fun = tpc_c6, color = "blue", size = 1) +
	stat_function(fun = tpc_c4, color = "orange", size = 1) +
	stat_function(fun = tpc_c8, color = "pink", size = 1) +
	stat_function(fun = tpc_c12, color = "green", size = 1) +
	stat_function(fun = tpc_c10, color = "purple", size = 1) +
	xlim(5, 45) + 
	ylim(0, 2.5) + geom_hline(yintercept = 0) +
 ylab("Exponential growth rate") + xlab("Temperature (°C)") +
	# geom_line(aes(x = temperature, y = growth, group = replicate), color = "cadetblue", data = all_preds, alpha = 0.2) +
	# geom_line(aes(x = temperature, y = growth, group = replicate), color = "purple", data = all_preds_v, alpha = 0.2) +
	# geom_line(aes(x = temperature, y = growth, group = replicate), color = "orange", data = all_preds_NLA, alpha = 0.2) +
	# geom_ribbon(aes(x = temperature, ymin = prediction_lower, ymax = prediction_upper, linetype = NA), fill = "transparent", alpha = 0.01,
	# data = all_preds_average, linetype = "dashed", color = "red", size = 0.5) +
	# geom_ribbon(aes(x = temperature, ymin = q2.5, ymax = q97.5, linetype=NA), data = limits_v, fill = "purple", alpha = 0.5) +
	coord_cartesian()


fit_function <- nlsLM(RFU ~ 20 * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
					  data= p6,  
					  start=list(z= results$z[[1]],w= results$w[[1]],a= results$a[[1]], b= results$b[[1]]),
					  lower = c(z = 0, w= 0, a = -0.2, b = 0),
					  upper = c(z = 30, w= 80,a =  0.5, b = 0.15),
					  control = nls.control(maxiter=1024, minFactor=1/204800000))








# nls multstart -----------------------------------------------------------
# 
# # start= list(z= 25,w= 25,a= 0.2, b= 0.1),
# lower = c(z = 0, w= 0, a = -0.2, b = 0),
# upper = c(z = 40, w= 80,a =  2, b = 2),

library(nls.multstart)
output2 <- rfu2 %>% 
	group_by(population) %>% 
	nest() %>% 
	mutate(fit = purrr::map(data, ~ nls_multstart(RFU ~ 20 * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
												  data = .x,
												  iter = 500,
												  start_lower = c(z= 25,w= 25,a= 0.2, b= 0.1),
												  start_upper = c(z= 35,w= 35,a= 0.4, b= 0.2),
												  supp_errors = 'N',
												  na.action = na.omit,
												  lower = c(z = 0, w= 0, a = -0.2, b = 0),
												  upper = c(z = 40, w= 80,a =  2, b = 2),
												  control = nls.control(maxiter=1000, minFactor=1/204800000))))


data_sub <- df_split[[4]] %>% 
	filter(round == "single")
fit1 <- nls_multstart(RFU ~ 18 * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
			  data = data_sub,
			  iter = 500,
			  start_lower = c(z= 25,w= 25,a= 0.0, b= 0.1),
			  start_upper = c(z= 35,w= 35,a= 0.3, b= 0.2),
			  supp_errors = 'N',
			  na.action = na.omit,
			  lower = c(z = 0, w= 0, a = -0.2, b = -2),
			  upper = c(z = 40, w= 80,a =  1, b = 2),
			  control = nls.control(maxiter=1000, minFactor=1/204800000))

info <- output2 %>%
	unnest(fit %>% map(glance))
params <- output2 %>%
	unnest(fit %>% map(tidy))

warnings()

