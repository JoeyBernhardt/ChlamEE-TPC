library(tidyverse)
library(lubridate)


all_fcs2_all_25 <- read_csv("data-processed/chlamee-acclimated-particles.csv") %>% 
	mutate(sample_volume = 25)
all_fcs2_all_50 <- read_csv("data-processed/chlamee-acclimated-particles-50ul.csv") %>% 
	mutate(sample_volume = 50)
plate_info <- read_csv("data-processed/chlamee-acclimated-plate-info.csv")


all_fcs <- bind_rows(all_fcs2_all_25, all_fcs2_all_50) %>% 
	dplyr::filter(plate ==70)


all_fcs3_all <- all_fcs %>% 
	select(4:11, sample_volume) %>% 
	dplyr::filter(fl1_a > 5, fl3_a > 5) %>% 
	mutate(well = str_to_upper(well))
# all_fcs3_all25 <- all_fcs2_all %>% 
# 	select(4:11) %>% 
# 	dplyr::filter(fl1_a > 5, fl3_a > 5) %>% 
# 	mutate(well = str_to_upper(well))

# sorted_all25 <- all_fcs3_all25 %>% 
# 	mutate(type = NA) %>% 
# 	mutate(type = ifelse(fl3_a > 85000, "algae", type)) %>% 
# 	mutate(type = ifelse(is.na(type), "background", type))


sorted_all <- all_fcs3_all %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a > 85000, "algae", type)) %>% 
	mutate(type = ifelse(is.na(type), "background", type))

all_sorted_all <- left_join(sorted_all, plate_info, by = c("well", "plate")) %>% 
	# dplyr::filter(!grepl("A", well)) %>% 
	mutate(plate = as.integer(plate))

# all_sorted_all25 <- left_join(sorted_all25, plate_info, by = c("well", "plate")) %>% 
# 	# dplyr::filter(!grepl("A", well)) %>% 
# 	mutate(plate = as.integer(plate))


all_sorted_all25 %>% 
	dplyr::filter(plate == 74) %>% 
	dplyr::filter(grepl("G", well)) %>% 
	ggplot(aes(x = fl4_a, y = fl3_a, color = type)) + geom_point() +
	facet_wrap( ~ well) + scale_y_log10() + scale_x_log10() + geom_hline(yintercept = 85000)

all_sorted_all %>% 
	dplyr::filter(plate == 70) %>% 
	dplyr::filter(grepl("E", well)) %>% 
	ggplot(aes(x = fl1_a, y = fl3_a, color = type)) + geom_point() +
	facet_wrap( ~ well) + scale_y_log10() + scale_x_log10() + geom_hline(yintercept = 85000)


# str(all_sorted_all)
# 
# all_sorted2 <- left_join(all_sorted_all, plate_key, by = "plate")

all_algae <- all_sorted_all %>% 
	dplyr::filter(type == "algae")
# all_algae_25 <- all_sorted_all25 %>% 
# 	dplyr::filter(type == "algae")


counts <- all_algae %>% 
	group_by(plate, well, temperature, sample_volume) %>% 
	tally() %>% 
	mutate(cells_per_ml = ifelse(sample_volume == 25, n*40*(1/0.3255), n*20*(1/0.3255))) 


counts %>% 
	select(-n) %>% 
	spread(key = sample_volume, value = cells_per_ml) %>% 
	mutate(`50` = ifelse(is.na(`50`), 0, `50`)) %>% 
	mutate(`25` = ifelse(is.na(`25`), 0, `25`)) %>% 
	# dplyr::filter(`25` == 0) %>% 
	ggplot(aes(x = `50`, y = `25`)) + geom_point() + geom_abline(intercept = 0, slope = 1) +
	geom_abline(intercept = 0, slope = 2, color = "purple") +
	ylab("Cell count with sample volume 25") + xlab("Cell count with sample volume 50") 
ggsave("figures/sample_volume_comparison.pdf", width = 6, height = 6)
	



counts_25 <- all_algae_25 %>% 
	group_by(plate, well, temperature) %>% 
	tally() %>% 
	mutate(cells_per_ml = n*40*(1/0.3255))

counts3 <- counts_25 %>% 
	mutate(sample = ifelse(is.na(temperature), "blank", "phyto"))


counts3 %>% 
	ggplot(aes(x = well, y = cells_per_ml, color = sample)) + geom_point() +
	facet_wrap( ~ plate, scales = "free")

counts25 <- counts3 %>% 
	dplyr::filter(plate == 70) %>% 
	mutate(sample_volume = 25)
	
counts50 <- counts2 %>% 
	dplyr::filter(plate == 70) %>% 
	mutate(sample_volume = 50)


all_counts <- bind_rows(counts25, counts50)

counts_wide <- left_join(counts50, counts25, by = c("well", "plate")) %>% 
	rename(counts_50 = n.x,
		   counts_25 = n.y)


counts_wide %>% 
	ggplot(aes(x = counts_50, y = counts_25))  + geom_point() +
	geom_abline(intercept = 0, slope = 1)
	

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
	dplyr::filter(plate == 70, sample_volume == 50) %>% 
	# dplyr::filter(date_time == ymd_hms("2018-11-13 13:31:49")) %>%
	ggplot(aes(x = cells_per_ml, y = RFU, group = temperature, color = factor(temperature))) +
	geom_point() 

all %>% 
	dplyr::filter(!is.na(population)) %>% 
	dplyr::filter(plate == 70, sample_volume == 50) %>% 
	# dplyr::filter(date_time == ymd_hms("2018-11-13 13:31:49")) %>%
	ggplot(aes(x = cells_per_ml, y = RFU, group = sample_volume, color = factor(sample_volume))) +
	geom_point() 
ggsave("figures/sample_volume_rfu_v_counts.pdf", width = 10, height =6)

all2 %>% 
	ungroup() %>% 
	distinct(temperature, plate) %>% View

all4 <- all %>% 
	dplyr::filter(!is.na(population)) %>% 
	dplyr::filter(plate %in% c(74, 75)) %>% 
	mutate(keep = case_when(plate == 74 & date_time == ymd_hms("2018-11-13 13:31:49") ~ "yes",
							plate == 75 ~ "yes",
							TRUE ~ "no")) %>% 
	dplyr::filter(keep == "yes") 

all4%>% 
	ggplot(aes(x = date_time, y = cells_per_ml, group = temperature, color = factor(plate))) + geom_point() +
	facet_wrap( ~ population, scales = "free") 
ggsave("figures/cells_time.pdf", width = 16, height = 16)

library(broom)
library(tidyverse)

write_csv(all4, "data-processed/chlamee-acclimated-cells-RFU-test.csv")

all4 %>% 
	dplyr::filter(population == 27) %>% 
	lm(log(RFU) ~ days, data = .) %>% 
	tidy()

all4 %>% 
	dplyr::filter(population == 27) %>% 
	lm(log(cells_per_ml) ~ days, data = .) %>% 
	tidy()

all2 %>% 
	ggplot(aes(x = cells_per_ml, y = RFU, color = factor(temperature))) + geom_point() +
	geom_smooth()
