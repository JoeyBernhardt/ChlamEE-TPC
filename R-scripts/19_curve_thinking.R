
library(vegan)

all_preds

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
