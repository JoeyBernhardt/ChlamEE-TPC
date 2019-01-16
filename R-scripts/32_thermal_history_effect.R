

library(tidyverse)



acclimated_curves <- read_csv("data-processed/chlamee-acclimated-tpc-curves.csv") %>% 
	mutate(thermal_history = "acclimated")


acute_curves <- read_csv("data-processed/chlamee-acute-tpc-curves.csv") %>% 
	mutate(thermal_history = "acute")


all_curves <- bind_rows(acclimated_curves, acute_curves)

str(all_curves)

all_curves %>% 
	dplyr::filter(!is.na(population)) %>% 
	dplyr::filter(!is.na(treatment)) %>% 
	ggplot(aes(x = temperature, y = growth, color = thermal_history, group = population)) + geom_point(size = 0.5) +
	ylim(0, 4.5) + xlim(0, 50) + geom_hline(yintercept = 0) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") + 
	scale_color_discrete(name = "Thermal history")  +
	facet_grid(treatment ~ ancestor_id)
ggsave("figures/chlamee-thermal-history-effect.pdf", width = 15, height = 10)
