
library(vegan)


curves_globe <-read_csv("data-processed/all_preds_not_acclimated.csv") %>% 
	mutate(history = "not acclimated")
curves_acc <- read_csv("data-processed/all_preds_acclimated.csv") %>% 
	mutate(history = "acclimated")

all_curves <- bind_rows(curves_globe, curves_acc)


all_curves %>% 
	# filter(population == 14) %>% 
	mutate(population = as.integer(population)) %>% 
	ggplot(aes(x = temperature, y = growth, color = history)) + geom_line(size = 1) +
	ylim(0, 4.3) + xlim(0, 50) +
	facet_wrap( ~ population) + geom_hline(yintercept = 0) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") +
	scale_color_discrete(name = "Thermal history")

ggsave("figures/acclimation_effect.pdf", width = 10, height = 8)





growth_16 <- all_preds %>% 
	spread(key = temperature, value = growth, 2:3) %>% 
	select(-population)

pca16 <- rda(growth_16, scale=TRUE)

loadings16 <- scores(pca16,choices=c(1,2))
summary(eigenvals(pca16))


tpc_temps <- unique(all_preds$temperature)
pcs16 <- as_data_frame((loadings16[[1]]))
pc1_16 <- pcs16 %>% 
	mutate(temperature = tpc_temps) 
pc1_16 %>% 
	ggplot(aes(x = temperature, y = PC1)) + geom_point() +
	xlim(0, 40) + 
	# geom_smooth() + 
	geom_hline(yintercept = 0) +
	xlab("Temperature (°C)") + geom_line()

ggsave("figures/PCA-non-acclimated.pdf", width = 6, height = 4)
