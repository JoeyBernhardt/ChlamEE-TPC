library(tidyverse)
library(tidyr)
library(broom)
library(cowplot)

all4 <- read_csv("data-processed/chlamee-acclimated-cells-RFU-test.csv")


rfu_growth <- all4 %>% 
	dplyr::group_by(population) %>% 
	do(tidy(lm(log(RFU) ~ days, data = .), conf.int = TRUE)) %>% 
	filter(term == "days") %>% 
	mutate(method = "rfu")

count_growth <- all4 %>% 
	dplyr::group_by(population) %>% 
	do(tidy(lm(log(cells_per_ml) ~ days, data = .), conf.int = TRUE)) %>% 
	filter(term == "days") %>% 
	mutate(method = "cell_counts")


all_growth <-bind_rows(rfu_growth, count_growth)

all_growth %>% 
	ggplot(aes(x = population, y = estimate, color = method)) + geom_point() +
	geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.2) +
	ylab("Growth rate estimate")
ggsave("figures/chlamee-22C-growth-rate-comparison.pdf", width = 8, height = 6)

all4 %>% 
	filter(days < 1) %>% 
	ggplot(aes(x = days, y = RFU)) + geom_point() + 
	facet_wrap( ~ population, scales = "free") 

all4 %>% 
	filter(days < 1) %>% 
	ggplot(aes(x = n, y = RFU)) + geom_point() 
ggsave("figures/RFU-v-events-22C.pdf", width = 6, height = 4)

all4 %>% 
	# filter(days < 1) %>% 
	ggplot(aes(x = days, y = log(RFU))) + geom_point() + 
	facet_wrap( ~ population)  + geom_smooth(method = "lm")

all4 %>% 
	filter(days < 1) %>% 
	ggplot(aes(x = days, y = RFU)) + geom_point() + 
	facet_wrap( ~ population) 


all4 %>% 
	gather(key = "measurement_type", value = "abundance", RFU, cells_per_ml) %>% 
	group_by(measurement_type) %>% 
	ggplot(aes(x = days, y = log(abundance), color = measurement_type)) + geom_point() + 
	facet_wrap( ~ population)  + geom_smooth(method = "lm")
ggsave("figures/abundance_v_time_chlamee_acclimated_22C_test.pdf", width = 10, height = 8)
