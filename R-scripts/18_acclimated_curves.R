

all_preds_acclimated <- read_csv("data-processed/all_preds_acclimated.csv") %>% 
	mutate(type = "acclimated")
not_acclimated <- read_csv("data-processed/all_preds_not_acclimated.csv") %>% 
	mutate(type = "not acclimated")


all_curves <- bind_rows(all_preds_acclimated, not_acclimated)


all_curves %>% 
	# filter(population == 14) %>% 
	mutate(population = as.integer(population)) %>% 
	ggplot(aes(x = temperature, y = growth, color = type)) + geom_line(size = 1) +
	ylim(0, 3.5) + xlim(0, 50) +
	facet_wrap( ~ population) + geom_hline(yintercept = 0) +
	ylab("Exponential growth rate") + xlab("Temperature (°C)") +
	scale_color_discrete(name = "Thermal history")

ggsave("figures/acclimation_effect.pdf", width = 10, height = 8)



