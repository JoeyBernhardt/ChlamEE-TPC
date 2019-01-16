

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
	ungroup() %>% 
	ggplot(aes(x = temperature, y = growth, color = thermal_history)) + geom_line(size = 1) +
	ylim(0, 5) + xlim(0, 45) + geom_hline(yintercept = 0) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") + 
	scale_color_viridis_d(name = "Thermal history", end =0.7, begin = 0.2)  +
	facet_grid(treatment ~ ancestor_id)
ggsave("figures/chlamee-thermal-history-effect.pdf", width = 15, height = 10)


acute <- read_csv("data-processed/chlamee-acute-tpc-fits.csv") %>%
rename(z = s) %>% 
	mutate(thermal_history = "acute")

acclimated <- read_csv("data-processed/chlamee-acclimated-tpc-fits-2.csv") %>% 
	rename(z = s) %>% 
	mutate(thermal_history = "acclimated")


both <- bind_rows(acute, acclimated)


both %>%
	ggplot(aes(x = population, y = topt, color = thermal_history)) + geom_point() +
	facet_grid(ancestor_id ~ treatment)

both_wide <- left_join(acute, acclimated, by = "population")

both_wide %>%
	ggplot(aes(x = topt.x, y = topt.y)) + geom_point() +
	geom_abline(slope = 1, intercept = 0) + ylab("Topt acclimated") + xlab("Topt acute")

