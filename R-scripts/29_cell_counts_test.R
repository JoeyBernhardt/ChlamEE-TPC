library(tidyverse)
library(lubridate)
library(cowplot)


dilution_fcs <- read_csv("data-processed/dilution-test-low-density-particles.csv") 
dilution_high <- read_csv("data-processed/dilution-test-high-density-particles.csv") 

plate_info <- read_csv("data-processed/chlamee-acclimated-plate-info.csv")



all_fcs3_all <- dilution_high %>% 
	select(3:10, concentration) %>% 
	dplyr::filter(fl1_a > 5, fl3_a > 5) %>%
	mutate(well = str_to_upper(well))


sorted_all <- all_fcs3_all %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a > 85000, "algae", type)) %>% 
	mutate(type = ifelse(is.na(type), "background", type))

# all_sorted_all <- left_join(sorted_all, plate_info, by = c("well", "plate")) %>% 
# 	# dplyr::filter(!grepl("A", well)) %>% 
# 	mutate(plate = as.integer(plate))


all_algae <- sorted_all %>% 
	dplyr::filter(type == "algae")
# all_algae_25 <- all_sorted_all25 %>% 
# 	dplyr::filter(type == "algae")


sorted_all %>% 
	dplyr::filter(grepl("12", well)) %>% 
	ggplot(aes(x = fl1_a, y = fl3_a, color = type)) + geom_point() +
	facet_wrap( ~ well) + scale_y_log10() + scale_x_log10()


counts <- all_algae %>% 
	group_by(well) %>% 
	tally() %>% 
	mutate(cells_per_ml = n*40) %>% 
	dplyr::filter(!grepl("01", well))



# bring in rfus -----------------------------------------------------------

rfu_low_density <- read_csv("data-processed/dilution-test-low-density-rfu.csv") %>% 
	dplyr::filter(!grepl("01", well))
rfu_high_density <- read_csv("data-processed/dilution-test-high-density-rfu.csv") %>% 
	dplyr::filter(!grepl("01", well))


all_counts_rfu <- left_join(counts, rfu_high_density, by = "well") %>% 
	separate(col = well, into =c("row", "column"), remove = FALSE, sep = 1)


all_counts_rfu %>% 
	# dplyr::filter(RFU < 100) %>% 
	mutate(scale_RFU = scale(RFU)) %>% 
	mutate(scale_counts = scale(cells_per_ml)) %>% 
	# dplyr::filter(population == 25) %>% 
	ggplot(aes(x = RFU, y = cells_per_ml, color = column)) + geom_point() +
	geom_smooth(method = "lm", color = "black") + facet_grid( ~ population, scales = "free") 

mod1 <- lm(RFU ~ cells_per_ml, data = all_counts_rfu)
summary(mod1)


pop25_low <- all_counts_rfu %>% 
	dplyr::filter(RFU < 200, population == "25") 

mod1 <- lm(RFU ~ cells_per_ml, data = pop25_low)
summary(mod1)
