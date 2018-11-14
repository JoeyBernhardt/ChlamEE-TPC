

library(tidyverse)
library(readxl)
library(janitor)
library(broom)

rfu <- read_csv("data-processed/chlamee-acclimated-rfu-time.csv")
population_key <- read_excel("data-general/ChlamEE_Treatments_JB.xlsx") %>% 
	clean_names()

rfu2 <- left_join(rfu, population_key, by = c("population"))

temp22 <- rfu2 %>% 
	filter(round == "repeat", temperature %in% c(22, 28, 34, 10, 40)) %>% 
	mutate(exponential = case_when(temperature == 28 & days < 1 ~ "yes",
								   temperature == 34 & days < 1 ~ "yes",
								   temperature == 22 & days < 2 ~  "yes",
								   temperature == 10 & days < 15 ~  "yes",
								   temperature == 40 & days < 15 & days > 1 ~  "yes",
								   TRUE  ~ "no")) %>% 
	group_by(temperature, population, well_plate) %>% 
	mutate(N0 = RFU[[1]]) %>% 
	filter(exponential == "yes")

write_csv(temp22, "data-processed/chlamee-exponential.csv")

unique(temp22$temperature)

temp22 %>% 
	# filter(temperature == 34, days < 1.2) %>%
	filter(exponential == "yes") %>% 
	ggplot(aes(x = days, y = RFU, color = temperature, group = well_plate)) + geom_point() +
	geom_line() +
	facet_wrap( ~ population, scales = "free")

temp22 %>% 
	ggplot(aes(x = days, y = RFU, color = treatment, group = well_plate)) + geom_point() +
	geom_line() +
	facet_wrap( ~ population)
ggsave("figures/exp_growth22.pdf", width = 15, height = 8)


temp22 %>% 
	ggplot(aes(x = treatment, y = N0, color = treatment, group = well_plate)) + geom_point() +
	geom_line() +
	facet_wrap( ~ population)


growth_rates <- temp22 %>%
	filter(exponential == "yes") %>% 
	group_by(temperature, population, well_plate) %>% 
	do(tidy(nls(RFU ~ N0 * exp(r*days),
				data= .,  start=list(r=0.01),
				control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() 


growth2 <- left_join(growth_rates, population_key, by = "population")

growth2 %>% 
	filter(population != "cc1629", temperature == 28) %>% 
	filter(!is.na(treatment)) %>% 
	ggplot(aes(x = reorder(treatment, estimate), y = estimate, fill = treatment)) + geom_boxplot() +
	ylab("Exponential growth rate (per day)") + xlab("Selection treatment") +
	facet_wrap( ~ ancestor_id + temperature, scales = "free")
ggsave("figures/growth_at_22_pooled.pdf", width = 8, height = 6)


growth2 %>% 
	filter(population != "cc1629", treatment != "C") %>% 
	filter(!is.na(treatment)) %>% 
	ggplot(aes(x = reorder(treatment, estimate), y = estimate, fill = treatment)) + geom_boxplot() +
	ylab("Exponential growth rate (per day)") + xlab("Selection treatment") +
	facet_wrap( ~ temperature + ancestor_id, scales = "free")

growth2 %>% 
	filter(population != "cc1629") %>% 
	filter(!is.na(treatment)) %>% 
	ggplot(aes(x = reorder(treatment, estimate), y = estimate, fill = treatment)) + geom_boxplot() +
	ylab("Exponential growth rate at 22C (per day)") + xlab("Selection treatment") +
	facet_wrap( ~ temperature, scales = "free")
