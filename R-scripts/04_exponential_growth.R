### Fitting exponential growth
library(lubridate)
library(broom)


exp_days <- read_csv("data-processed/exponential_repeats_RFU_anc4.csv")

cell_days_c <- exp_days %>% 
	mutate(start_time = min(date_time)) %>% 
	mutate(days = interval(start_time, date_time)/ddays(1))

cell_key <- cell_days_c %>% 
	select(temperature, Treatment, population_id, unique_id, Ancestor_ID) %>% 
	distinct(population_id, .keep_all = TRUE)


exp_days %>% 
	group_by(population_id) %>% 
	summarise(min_RFU = min(RFU)) %>%
	ungroup() %>% 
	summarise_each(funs(min, max, mean), min_RFU) %>% View


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
	geom_hline(yintercept = 0)
ggsave("figures/anc4-TPC.pdf", width = 6, height = 4)



