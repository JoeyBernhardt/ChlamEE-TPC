library(tidyverse)
library(lubridate)


all_fcs2_all <- read_csv("data-processed/chlamee-acclimated-particles.csv")
plate_info <- read_csv("data-processed/chlamee-acclimated-plate-info.csv")




plate2 <- plate_layout %>% 
	mutate(well = str_to_upper(well)) 

all_fcs3_all <- all_fcs2_all %>% 
	select(4:10) %>% 
	dplyr::filter(fl1_a > 5, fl3_a > 5) %>% 
	mutate(well = str_to_upper(well))

sorted_all <- all_fcs3_all %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a > 1250, "algae", type)) %>% 
	mutate(type = ifelse(is.na(type), "background", type))

all_sorted_all <- left_join(sorted_all, plate_info, by = c("well", "plate")) %>% 
	# dplyr::filter(!grepl("A", well)) %>% 
	mutate(plate = as.integer(plate))


all_sorted_all %>% 
	dplyr::filter(plate == 74) %>% 
	ggplot(aes(x = fl1_a, y = fl3_a, color = type)) + geom_point() +
	facet_wrap( ~ well) + scale_y_log10() + scale_x_log10()


# str(all_sorted_all)
# 
# all_sorted2 <- left_join(all_sorted_all, plate_key, by = "plate")

all_algae <- all_sorted_all %>% 
	dplyr::filter(type == "algae")


counts <- all_algae %>% 
	group_by(plate, well, temperature) %>% 
	tally() %>% 
	mutate(cells_per_ml = n*40*(1/0.3255))

counts2 <- counts %>% 
	mutate(sample = ifelse(is.na(temperature), "blank", "phyto"))


	# dplyr::filter(plate == 70) %>% 
	ggplot(aes(x = well, y = cells_per_ml, color = sample)) + geom_point() +
	facet_wrap( ~ plate, scales = "free")




### lets see how much carry over there is

counts %>% 
	separate(well, into = c("row", "column"), sep = 1, remove = FALSE) %>% 
	# dplyr::filter(row == "A") %>% 
	ggplot(aes(x = temperature, y = cells_per_ml, color = row)) + geom_point() +
	facet_wrap( ~ row, scales = "free")



# bring in RFU data -------------------------------------------------------

RFUS <- read_csv("data-processed/chlamee-acclimated-RFU-time.csv")

rfu2 <- RFUS %>% 
	dplyr::filter(date_time > ymd_hms("2018-10-11 18:00:00")) 

rfu2 %>% 
	ggplot(aes(x = date_time, y = RFU)) + geom_point() + facet_wrap( ~ plate)


all <- left_join(counts, rfu2, by = c("plate", "well", "temperature"))

all2 <- bind_rows(all, all_inocs2) 

unique(counts$plate)


all %>% 
	dplyr::filter(!is.na(population)) %>% 
	dplyr:: filter(plate == 70) %>% 
	# filter(population == 4) %>% 
	ggplot(aes(x = well, y = cells_per_ml, group = temperature, color = factor(temperature))) + geom_point() +
	facet_wrap( ~ population) 

all %>% 
	dplyr::filter(!is.na(population)) %>% 
	dplyr:: filter(plate == 74) %>% 
	dplyr::filter(date_time == ymd_hms("2018-11-13 13:31:49")) %>%
	# filter(population == 4) %>% 
	ggplot(aes(x = cells_per_ml, y = RFU, group = temperature, color = factor(temperature))) + geom_point() 


all2 %>% 
	ungroup() %>% 
	distinct(temperature, plate) %>% View

all2 %>% 
	dplyr::filter(!is.na(population)) %>% 
	dplyr::filter(temperature == 10) %>% 
	ggplot(aes(x = date_time, y = cells_per_ml, group = temperature, color = factor(plate))) + geom_point() +
	facet_wrap( ~ population, scales = "free") 
ggsave("figures/cells_time.pdf", width = 16, height = 16)


all2 %>% 
	ggplot(aes(x = cells_per_ml, y = RFU, color = factor(temperature))) + geom_point() +
	geom_smooth()
