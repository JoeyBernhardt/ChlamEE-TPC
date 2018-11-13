

library(tidyverse)
library(readxl)
library(janitor)
library(broom)

rfu <- read_csv("data-processed/chlamee-acclimated-rfu-time.csv")
population_key <- read_excel("data-general/ChlamEE_Treatments_JB.xlsx") %>% 
	clean_names()

rfu2 <- left_join(rfu, population_key, by = c("population"))

temp22 <- rfu2 %>% 
	filter(round == "repeat", temperature == 22)

temp22 %>% 
	ggplot(aes(x = days, y = RFU, color = treatment, group = well_plate)) + geom_point() +
	geom_line() +
	facet_wrap( ~ population)
ggsave("figures/exp_growth22.pdf", width = 15, height = 8)

growth_rates <- temp22 %>% 
	group_by(population, well_plate) %>% 
	do(tidy(nls(RFU ~ 10 * exp(r*days),
				data= .,  start=list(r=0.01),
				control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() 


growth2 <- left_join(growth_rates, population_key, by = "population")

growth2 %>% 
	filter(population != "cc1629") %>% 
	ggplot(aes(x = reorder(treatment, estimate), y = estimate, fill = treatment)) + geom_boxplot() +
	ylab("Exponential growth rate at 22C (per day)") + xlab("Selection treatment") +
	facet_wrap( ~ ancestor_id)
ggsave("figures/growth_at_22_pooled.pdf", width = 8, height = 6)
