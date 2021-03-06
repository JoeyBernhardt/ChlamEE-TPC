



### RFU import for chlamee-acclimated

library(tidyverse)
library(cowplot)
library(readxl)
library(janitor)
library(lubridate)


# read in general data ----------------------------------------------------


plate_layout <- read_excel("data-general/chlamee-acclimated-plate-layout.xlsx") 

plate_key <- read_excel("data-general/chlamee-acclimated-plate-key.xlsx") 

plate_info <- left_join(plate_layout, plate_key, by = c("plate_key")) %>% 
	rename(column = colum) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>%
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	mutate(population = ifelse(population == "anc 2", "Anc 2", population))


write_csv(plate_info, "data-processed/chlamee-acclimated-plate-info.csv")


# read in RFU data --------------------------------------------------------
RFU_files <- c(list.files("data-raw/chlamee-RFU-acclimated", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xlsx", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "")

# RFU_files[grepl("104", RFU_files)]

# RFU_files <- RFU_files[!grepl("acc", RFU_files)]

all_plates <- map_df(RFU_files, read_excel, range = "B56:N64", .id = "file_name") %>%
	rename(row = X__1) %>% 
	filter(!grepl("dilution", file_name)) 

all_times <- map_df(RFU_files, read_excel, range = "A6:B8", .id = "file_name") %>% 
	clean_names() %>% 
	filter(!is.na(plate_1)) %>% 
	filter(!grepl("dilution", file_name)) %>% 
	spread(key = plate_number, value = plate_1) %>% 
	separate(Time, into = c("crap", "time"), sep = " ") %>% 
	select(-crap) %>% 
	separate(file_name, into = c("path", "plate"), sep = "plate", remove = FALSE)


all_plates2 <- dplyr::left_join(all_plates, all_times, by = "file_name")

all_temp_RFU <- all_plates2 %>% 
	gather(key = column, value = RFU, 3:14) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>%
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU)) %>% 
	mutate(plate = as.numeric(plate))


plate2 <- plate_layout %>%
	rename(column = colum) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") 




all_rfus_raw <- left_join(all_temp_RFU, plate_info, by = c("well", "plate")) %>% 
	mutate(plate = as.numeric(plate)) %>% 
	filter(!is.na(plate_key))



all_rfus2 <- all_rfus_raw %>%
	unite(col = date_time, Date, time, sep = " ") %>%
	mutate(date_time = ymd_hms(date_time))

all_rfus3 <- all_rfus2 %>% 
	mutate(round = ifelse(plate %in% c(90, 97, 104, 6, 13, 20, 27, 29, 36, 48, 55, 62, 69, 76, 83, 111, 118, 125), "repeat", "single")) %>% 
	group_by(temperature) %>% 
	mutate(start_time = min(date_time)) %>% 
	mutate(days = interval(start_time, date_time)/ddays(1)) %>% 
	unite(col = well_plate, well, plate, remove =  FALSE) 


write_csv(all_rfus3, "data-processed/chlamee-acclimated-rfu-time.csv")

unique(all_rfus3$temperature)

all_rfus3 %>%
	filter(round == "repeat") %>% 
	filter(temperature %in% c(16)) %>% 
	# filter(population == 30) %>% 
	ggplot(aes(x = days, y = RFU, color = factor(temperature), group = well_plate)) +
	geom_point(size = 2) +
	scale_color_viridis_d(name = "Temperature") +
	xlab("Days") +
	facet_wrap( ~ population, scales = "free_y") +
	geom_line() 


ggsave("figures/chlamee-acclimated-RFU-time-exponential.pdf", width = 16, height = 12)
ggsave("figures/chlamee-acclimated-RFU-time.pdf", width = 25, height = 12)
ggsave("figures/chlamee-acclimated-RFU-time-all.pdf", width = 25, height = 12)
ggsave("figures/chlamee-acclimated-RFU-time-repeats.pdf", width = 25, height = 12)
ggsave("figures/chlamee-acclimated-RFU-time-pop-facet.pdf", width = 25, height = 12)
ggsave("figures/chlamee-acclimated-RFU-time-34C.pdf", width = 14, height = 8)
ggsave("figures/chlamee-acclimated-RFU-time-28C.pdf", width = 14, height = 8)
ggsave("figures/chlamee-acclimated-RFU-time-22C.pdf", width = 14, height = 8)
ggsave("figures/chlamee-acclimated-RFU-time-40C.pdf", width = 14, height = 8)



# single repeats exponential ----------------------------------------------

single_exp <- all_rfus3 %>% 
	filter(round == "single", temperature %in% c(22, 28, 34, 10, 40, 16)) %>% 
	mutate(exponential = case_when(temperature == 28 & days < 1 ~ "yes",
								   temperature == 34 & days < 1 ~ "yes",
								   temperature == 22 & days < 2 ~  "yes",
								   temperature == 10 & days < 7 ~  "yes",
								   temperature == 16 & days < 5 ~  "yes",
								   temperature == 40 & days < 7 & days > 1 ~  "yes",
								   TRUE  ~ "no")) %>% 
	filter(exponential == "yes")

single_plate_list_acclimated <- single_exp %>%
	filter(round == "single") %>% 
	distinct(plate)
write_csv(single_plate_list_acclimated, "data-processed/chlamee-acclimated-single-plate-exp.csv")


single_exp %>%
	ggplot(aes(x = days, y = RFU, color = factor(temperature), group = well_plate)) +
	geom_point(size = 2) +
	scale_color_viridis_d(name = "Temperature") +
	xlab("Days") +
	facet_wrap( ~ temperature, scales = "free_y") 
ggsave("figures/single_exp_acclimated.pdf", width = 8, height = 6)	


# plot for comparison with R-star ---------------------------------------------
population_key <- read_excel("data-general/ChlamEE_Treatments_JB.xlsx") %>% 
	clean_names()

temp22 <- all_rfus3 %>%
	filter(round == "repeat") %>% 
	filter(temperature %in% c(22)) 
temp22_2 <- left_join(temp22, population_key, by = "population")


temp22_2 %>% 
	filter(population != "cc1629", RFU < 400) %>% 
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>% 
	ggplot(aes(x = days, y = RFU, color = factor(temperature), group = well_plate)) +
	geom_point(size = 2) +
	scale_color_viridis_d(name = "Temperature") +
	xlab("Days") +
	facet_wrap( ~ well_plate, scales = "free") +
	# facet_grid(rows = vars(treatment), cols = vars(ancestor_id), scales = "free_y") +
	geom_line()
ggsave("figures/temp22_time_series_chlamee_TPC_facet_low_RFU.pdf", height = 30, width = 30)


temp22_2 %>% 
	filter(population != "cc1629") %>% 
	mutate(treatment = ifelse(is.na(treatment), "none", treatment)) %>% 
	mutate(day_int = round(days, digits = 0)) %>% 
	group_by(well_plate, day_int) %>% 
	sample_n(size = 1) %>% 
	ggplot(aes(x = days, y = RFU, color = factor(well), group = well_plate)) +
	geom_point(size = 2) +
	scale_color_viridis_d(name = "Temperature") +
	xlab("Days") +
	facet_wrap( ~ well_plate, scales = "free") +
	# facet_grid(rows = vars(treatment), cols = vars(ancestor_id), scales = "free_y") +
	geom_line() + theme(legend.position = "none")
ggsave("figures/temp22_time_series_chlamee_TPC_one_time_per_day_facet.pdf", height = 30, width = 30)
