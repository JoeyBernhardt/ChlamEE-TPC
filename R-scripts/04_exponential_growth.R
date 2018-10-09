### Fitting exponential growth
library(lubridate)
library(broom)
library(tidyverse)


exp_days <- read_csv("data-processed/exponential_repeats_RFU_anc4.csv")
exp_cells <- read_csv("data-processed/all_temps_cells_RFU.csv")

cell_days_counts <- exp_cells %>% 
	mutate(start_time = min(date_time)) %>% 
	mutate(days = interval(start_time, date_time)/ddays(1)) %>% 
	unite(col = temp_treatment, temperature, Treatment, remove = FALSE)

cell_days_RFU <- exp_days %>% 
	mutate(start_time = min(date_time)) %>% 
	mutate(days = interval(start_time, date_time)/ddays(1)) %>% 
	unite(col = temp_treatment, temperature, Treatment, remove = FALSE)

cell_key <- cell_days_c %>% 
	unite(col = temp_well, temperature, well, remove = FALSE) %>% 
	select(temperature, Treatment, population_id, unique_id, Ancestor_ID, well) %>% 
	distinct(population_id, .keep_all = TRUE)

cell_days_counts %>% 
	filter(temperature == 35) %>% 
	ggplot(aes(x = days, y = cells_per_ml, color = temp_treatment, group = temp_treatment)) + geom_point() +
	facet_wrap( ~ temperature + Treatment, scales = "free") + theme(legend.position = "none") 


exp_cells %>% 
	group_by(population_id) %>% 
	summarise(min_RFU = min(RFU)) %>%
	ungroup() %>% 
	summarise_each(funs(min, max, mean), min_RFU) %>% View

### ok let's just fit a linear model to log transformed data, quick and dirty

counts_slopes <- cell_days_counts %>% 
	group_by(temperature, Treatment) %>% 
	do(tidy(lm(log(cells_per_ml + 1) ~ days, data = . ), conf.int = TRUE))

RFU_slopes <- cell_days_counts %>% 
	group_by(temperature, Treatment) %>% 
	do(tidy(lm(log(RFU + 1) ~ days, data = . ), conf.int = TRUE))


counts_slopes %>% 
	filter(term == "days") %>% 
	filter(!Treatment %in% c("COMBO", "P")) %>% 
	ggplot(aes(x = temperature, y = estimate, color = Treatment)) + geom_point() + 
	geom_smooth()


growth_rates_counts <- cell_days_counts %>% 
	# filter(temperature == 30) %>% 
	group_by(temperature, Treatment) %>% 
	do(tidy(nls(cells_per_ml ~ 214 * exp(r*days),
				data= .,  start=list(r=0.01),
				control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() 

growth_rates <- cell_days_c %>% 
group_by(population_id) %>% 
	do(tidy(nls(RFU ~ 8 * exp(r*days),
				data= .,  start=list(r=0.01),
				control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() 

growth2 <- left_join(growth_rates, cell_key, by = "population_id")

growth2 %>% 
	filter(!Treatment %in% c("COMBO", "P")) %>% 
	# filter(temperature > 10) %>% 
	ggplot(aes(x = temperature, y = estimate, color = Treatment, group = Treatment, fill = Treatment)) + geom_point() +
	ylab("Exponential growth rate (per day)") + xlab("Temperature (°C)") + geom_smooth() +
	geom_hline(yintercept = 0) +
	facet_wrap( ~ Treatment)
ggsave("figures/anc4-TPC.pdf", width = 6, height = 4)
ggsave("figures/anc4-TPC-together.pdf", width = 6, height = 4)


