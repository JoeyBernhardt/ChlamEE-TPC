

library(tidyverse)
library(readxl)
library(janitor)
library(broom)

theme_set(theme_cowplot())

rfu <- read_csv("data-processed/chlamee-acclimated-rfu-time.csv")
population_key <- read_excel("data-general/ChlamEE_Treatments_JB.xlsx") %>% 
	clean_names()

rfu2 <- left_join(rfu, population_key, by = c("population"))

temp22 <- rfu2 %>% 
	filter(round == "repeat", temperature %in% c(22, 28, 34, 10, 40, 16)) %>% 
	mutate(exponential = case_when(temperature == 28 & days < 1 ~ "yes",
								   temperature == 34 & days < 1 ~ "yes",
								   temperature == 22 & days < 2 ~  "yes",
								   temperature == 10 & days < 7 ~  "yes",
								   temperature == 16 & days < 3 ~  "yes",
								   temperature == 40 & days < 7 & days > 1 ~  "yes",
								   TRUE  ~ "no")) %>% 
	group_by(temperature, population, well_plate) %>% 
	mutate(N0 = RFU[[1]]) %>% 
	filter(exponential == "yes") %>% 
	mutate(treatment = ifelse(is.na(treatment), "none", treatment))

write_csv(temp22, "data-processed/chlamee-exponential.csv")

unique(temp22$temperature)

temp22 %>% 
	# filter(temperature == 34, days < 1.2) %>%
	filter(exponential == "yes") %>% 
	ggplot(aes(x = days, y = RFU, color = factor(temperature), group = well_plate)) + geom_point() +
	geom_line() +
	facet_grid(treatment ~ ancestor_id, scales = "free") + scale_color_viridis_d(name = "Temperature")
ggsave("figures/exp-growth-chlamee-acclimated.pdf", width = 15, height = 8)

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
	group_by(temperature, treatment, ancestor_id, population, well_plate) %>% 
	do(tidy(nls(RFU ~ N0 * exp(r*days),
				data= .,  start=list(r=0.01),
				control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() 


growth2 <- left_join(growth_rates, population_key, by = c("population","treatment","ancestor_id"))

growth2 %>% 
	filter(population != "cc1629") %>% 
	filter(!is.na(treatment)) %>% 
	ggplot(aes(x = reorder(treatment, estimate), y = estimate, fill = treatment)) + geom_boxplot() +
	ylab("Exponential growth rate (per day)") + xlab("Selection treatment") +
	facet_wrap( ~ ancestor_id + temperature, scales = "free")
ggsave("figures/growth_pooled.pdf", width = 25, height = 15)


growth2 %>% 
	filter(temperature == 22) %>% 
	lm(estimate ~ ancestor_id*treatment, data = .) %>% summary()
	tidy(., conf.int = TRUE) %>% View


growth2 %>% 
	filter(population != "cc1629") %>% 
	filter(!is.na(treatment)) %>% 
	ggplot(aes(x = reorder(treatment, estimate), y = estimate, fill = treatment)) + geom_boxplot() +
	ylab("Exponential growth rate (per day)") + xlab("Selection treatment") +
	facet_wrap( ~ temperature + treatment, scales = "free")

growth2 %>% 
	filter(population != "cc1629") %>% 
	filter(!is.na(treatment)) %>% 
	ggplot(aes(x = reorder(treatment, estimate), y = estimate, fill = treatment)) + geom_boxplot() +
	ylab("Exponential growth rate at 22C (per day)") + xlab("Selection treatment") +
	facet_wrap( ~ temperature, scales = "free")



# now estimate TPC --------------------------------------------------------

library(growthTools)


res <- growth2 %>% 
	# filter(population == 5) %>% 
	group_by(population, treatment, ancestor_id) %>% 
	do(tpcs=get.nbcurve.tpc(.$temperature,.$estimate,method='grid.mle2',plotQ=F,conf.bandQ=F))


clean.res <- res %>% 
	summarise(population, treatment, ancestor_id,topt=tpcs$o,tmin=tpcs$tmin,tmax=tpcs$tmax,rsqr=tpcs$rsqr, 
			  a = tpcs$a, b = tpcs$b, w = tpcs$w, s = tpcs$s)

write_csv(clean.res, "data-processed/chlamee-acclimated-tpc-fits-2.csv")


clean.res %>% 
	lm(topt ~ treatment, data = .) %>% summary()

clean.res %>% 
	ggplot(aes(x = treatment, y = topt)) + geom_boxplot()

library(nlme)
growth_mod <-lme(topt ~ treatment, random= ~1|ancestor_id, data= clean.res)
summary(growth_mod)
anova(growth_mod)
intervals(growth_mod)

temps <- seq(0,40, 1)
sapply(temps, nbcurve2(clean.res)
