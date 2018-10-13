
library(minpack.lm)
library(broom)
library(nlstools)
library(nls.multstart)

rfu <- read_csv("data-processed/globe-chlamy-exponential-RFU-time.csv")


fit <- nlsLM(RFU ~ 20 * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
	  data= filter(rfu2, population == 9),  
	  start= list(z= results$z[[1]],w= results$w[[1]],a= results$a[[1]], b= results$b[[1]]),
	  lower = c(z = 0, w= 0, a = -0.2, b = 0),
	  upper = c(z = 40, w= 80,a =  0.7, b = 0.15),
	  control = nls.control(maxiter=1024, minFactor=1/204800000))

d2 <- filter(rfu2, population == 10)

exponential_tpc <- function(a, b, z, w, temp, days) {
	RFU  <- 20 * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days))
	
}

fit <- nls_multstart(RFU ~ 20 * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
					 data = d2,
					 iter = 500,
					 start_lower = c(z = 0, w= 0, a = -0.2, b = 0),
					 start_upper = c(z = 30, w= 80, a =  0.5, b = 0.15),
					 supp_errors = 'N',
					 na.action = na.omit,
					 lower = c(z = 0, w= 0, a = -0.3, b = 0))



summary(fit)



p6 <- rfu %>% 
	filter(population == 6) %>% 
	rename(temp = temperature)

p1 <- rfu %>% 
	filter(population == 1) %>% 
	rename(temp = temperature)

p10 <- rfu %>% 
	filter(population == 10) %>% 
	rename(temp = temperature)

p4 <- rfu %>% 
	filter(population == 4) %>% 
	rename(temp = temperature)

p2 <- rfu %>% 
	filter(population == 2) 


growth_function <- function(temp, days)(20 * exp((results$a*exp(results$b*temp)*(1-((temp-results$z)/(results$w/2))^2))*(days)))

avals<-seq(-0.2,1.2,0.02)
bvals<-seq(-0.2,0.3,0.02)

avals<-seq(0,0.5,0.07)
bvals<-seq(0,0.16,0.08)
zvals<-seq(10,15,1)
wvals<-seq(34,37,1)


df <-  expand.grid(a = avals, b = bvals, z = zvals, w = wvals) %>% 
	mutate(unique_id = rownames(.)) %>% 
	mutate(sample_size = 1)


##cell_density ~ 800 * exp(r*days)

df <- df_split[[1]]
min(rfus, na.rm = TRUE)

fit_growth <- function(df){
	start_density <- df$RFU[[1]]
	res <- try(nlsLM(RFU ~ 20 * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
					 data= df,  
					 start= list(z= results$z[[1]],w= results$w[[1]],a= results$a[[1]], b= results$b[[1]]),
					 lower = c(z = 0, w= 0, a = -0.2, b = 0),
					 upper = c(z = 45, w= 70,a =  3, b = 1),
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
	split(.$population)

output <- df_split %>%
	map_df(fit_growth, .id = "population") 



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
	map_df(prediction_function, .id = "population")

all_preds %>% 
	mutate(population = as.integer(population)) %>% 
	ggplot(aes(x = temperature, y = growth, color = factor(population))) + geom_line(size = 1) +
	ylim(0, 3) + xlim(0, 45) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") + scale_color_discrete(name = "Population")



results <- output %>% 
	filter(df.residual > 5) %>%
	top_n(n = -1, wt = AIC)


rfu2 <- rfu %>% 
	rename(temp = temperature)


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