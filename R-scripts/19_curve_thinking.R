
library(vegan)
library(cowplot)
library(ggbiplot)
library(tidyverse)


curves_globe <-read_csv("data-processed/all_preds_not_acclimated.csv") %>% 
	mutate(history = "not acclimated")
curves_acc <- read_csv("data-processed/all_preds_acclimated.csv") %>% 
	mutate(history = "acclimated")

all_curves <- bind_rows(curves_globe, curves_acc)


all_curves %>% 
	filter(history == "not acclimated", population %in% c(10, 1, 12, 13, 14)) %>% 
	mutate(population = as.integer(population)) %>% 
	ggplot(aes(x = temperature, y = growth, color = factor(population))) + geom_line(size = 1) +
	ylim(0, 4.3) + xlim(0, 50) + geom_hline(yintercept = 0) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") 

all_curves %>% 
	filter(history == "not acclimated", population %in% c(1,14)) %>% 
	mutate(population = as.integer(population)) %>% 
	ggplot(aes(x = temperature, y = growth, color = factor(population))) + geom_line(size = 1) +
	ylim(0, 2.1) + xlim(0, 45) + geom_hline(yintercept = 0) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") +
	scale_color_viridis_d(name = "Population", begin = 0.8, end = 0.3) 
ggsave("figures/hot-cold-shift-example.pdf", width = 6, height = 3)


all_curves %>% 
	filter(history == "not acclimated", population %in% c(1,14, 11, 3)) %>% 
	mutate(population = as.integer(population)) %>% 
	ggplot(aes(x = temperature, y = growth, color = factor(population))) + geom_line(size = 1) +
	ylim(0, 3.5) + xlim(0, 45) + geom_hline(yintercept = 0) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") +
	scale_color_viridis_d(name = "Population", begin = 0.8, end = 0.3) 
all_curves %>% 
	filter(history == "not acclimated", population %in% c(11, 3, 5)) %>% 
	mutate(population = as.integer(population)) %>% 
	ggplot(aes(x = temperature, y = growth, color = factor(population))) + geom_line(size = 1) +
	ylim(0, 3.5) + xlim(0, 45) + geom_hline(yintercept = 0) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") +
	scale_color_viridis_d(name = "Population", begin = 0.8, end = 0.3) 


all_curves %>% 
	# filter(population == 14) %>% 
	mutate(population = as.integer(population)) %>% 
	ggplot(aes(x = temperature, y = growth, color = history)) + geom_line(size = 1) +
	ylim(0, 4.3) + xlim(0, 50) +
	facet_wrap( ~ population) + geom_hline(yintercept = 0) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") +
	scale_color_discrete(name = "Thermal history")

ggsave("figures/acclimation_effect.pdf", width = 10, height = 8)


library(tidyverse)

growth_16 <-curves_globe %>% 
	select(-history) %>% 
	filter(temperature %in% c(1:40)) %>% 
	# filter(temperature %in% c(10, 16, 22, 28, 34, 40)) %>% 
	spread(key = temperature, value = growth, 2:3) %>% 
	select(-population)

cov_mat <- cov(growth_16)
pca_res <- prcomp(cov_mat, center = TRUE,scale. = TRUE)
summary(pca_res)
ggbiplot(pca_res, labels= tpc_temps)
ggsave("figures/PCA-non-acclimated-biplot.pdf", width = 6, height = 6)

pca16 <- rda(growth_16)


loadings16 <- scores(pca_res,choices=c(1,2))
summary(eigenvals(pca_res))

tpc_temps <- c(1:40)
tpc_temps <- c(10, 16, 22, 28, 34, 40)
# tpc_temps <- unique(growth_16$temperature)

pcs16 <- as_data_frame((loadings16)) %>% 
	select(PC1)
pc1_16 <- pcs16 %>% 
	mutate(temperature = tpc_temps) 
pc1_16 %>% 
	ggplot(aes(x = temperature, y = PC1)) + geom_point() +
	xlim(0, 40) + 
	# geom_smooth() + 
	geom_hline(yintercept = 0) +
	xlab("Temperature (°C)") + geom_line() +
	ylab("PC1 loadings")

ggsave("figures/PCA-non-acclimated.pdf", width = 6, height = 4)



pcs216 <- as_data_frame((loadings16)) %>% 
	select(PC2)
pc2_16 <- pcs216 %>% 
	mutate(temperature = tpc_temps) 
pc2_16 %>% 
	ggplot(aes(x = temperature, y = PC2)) + geom_point() +
	xlim(0, 40) + 
	# geom_smooth() + 
	geom_hline(yintercept = 0) +
	xlab("Temperature (°C)") + geom_line() +
	ylab("PC2 loadings")
ggsave("figures/PCA2-non-acclimated.pdf", width = 6, height = 4)

# acclimated --------------------------------------------------------------
library(vegan)
growth_16 <- curves_acc %>% 
	select(-history) %>% 
	spread(key = temperature, value = growth, 2:3) %>% 
	select(-population)

pca16 <- rda(growth_16, scale=TRUE)

loadings16 <- scores(pca16,choices=c(1,2))
summary(eigenvals(pca16))


tpc_temps <- unique(curves_acc$temperature)
pcs16 <- as_data_frame((loadings16[[1]]))
pc1_16 <- pcs16 %>% 
	mutate(temperature = tpc_temps) 
pc1_16 %>% 
	ggplot(aes(x = temperature, y = PC1)) + geom_point() +
	xlim(0, 40) + 
	# geom_smooth() + 
	geom_hline(yintercept = 0) +
	xlab("Temperature (°C)") + geom_line()

ggsave("figures/PCA-acclimated.pdf", width = 6, height = 4)

