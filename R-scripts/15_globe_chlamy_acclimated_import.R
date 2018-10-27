# packages

library(tidyverse)
library(cowplot)
library(readxl)
library(janitor)
library(lubridate)


# read in general data ----------------------------------------------------


plate_layout <- read_excel("data-general/globe-acclimated-dilution-plate-1.xlsx") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	unite(col = well, row, column, remove = FALSE, sep = "") %>% 
	mutate(well = str_to_upper(well))
	
	
RFU_files <- c(list.files("data-raw/globe-chlamy-RFU-acclimated", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xlsx", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "")

all_plates <- map_df(RFU_files, read_excel, range = "B56:N64", .id = "file_name") %>% 
	rename(row = X__1) 
all_times <- map_df(RFU_files, read_excel, range = "A6:B8", .id = "file_name") %>% 
	clean_names() %>% 
	# filter(!is.na(plate_1)) %>% 
	# select(-plate_2) %>% 
	spread(key = plate_number, value = plate_1) %>% 
	separate(Time, into = c("crap", "time"), sep = " ") %>% 
	select(-crap) %>% 
	separate(file_name, into = c("path", "plate"), sep = "plate", remove = FALSE)

all_plates2 <- left_join(all_plates, all_times, by = "file_name")
all_temp_RFU <- all_plates2 %>% 
	gather(key = column, value = RFU, 3:14) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU)) %>% 
	mutate(plate = as.numeric(plate))

all_rfu <- left_join(all_temp_RFU, plate_layout, by = c("well", "plate"))  %>% 
	filter(!grepl("setup-", file_name)) %>% 
	filter(!is.na(temperature))

write_csv(all_rfu, "data-processed/RFU-globe-chlamy-setup-b-2018-10-10.csv")

dilutions <- all_rfu %>% 
	group_by(population, temperature) %>% 
	summarise(mean_rfu = mean(RFU)) %>% 
	mutate(dilution_full = mean_rfu/11) %>% 
	mutate(dilution2 = 1/dilution_full)

dilutions %>% 
	filter(temperature == 10) %>% View

high_conc <- dilutions %>% 
	# filter(dilution2 > 0.15) %>% 
	mutate(volume_needed = dilution2*15) %>% 
	mutate(combo_needed = 15-volume_needed) %>% 
	mutate(volume_needed = round(digits = 2, volume_needed)) %>% 
	mutate(combo_needed = round(digits = 2, combo_needed)) %>% 
	mutate(over_5 = ifelse(volume_needed > 5, "yes", "no")) %>% 
	mutate(over_10 = ifelse(volume_needed > 10, "yes", "no")) %>%
	mutate(over_5_divided_x2 = ifelse(over_5 == "yes", volume_needed/2, NA)) %>% 
	mutate(over_10_divided_x3 = ifelse(over_10 == "yes", volume_needed/5, NA)) %>% 
	mutate(combo_needed_divided_x2 = round(digits = 2, combo_needed/2)) 

high_conc %>% 
	filter(temperature == 34) %>% 
	select(-contains("over")) %>%
	select(population, volume_needed, contains("combo")) %>% View
	


low_conc <- dilutions %>% 
	filter(dilution2 < 0.15) %>% 
	mutate(one_in_two = mean_rfu/2) %>% 
	mutate(dilution_full2 = one_in_two/11) %>% 
	mutate(dilution3 = 1/dilution_full2) %>% 
	mutate(volume_needed = dilution3*30) %>% 
	mutate(combo_needed = 30-volume_needed) %>% 
	mutate(volume_needed = round(digits = 1, volume_needed)) %>% 
	mutate(combo_needed = round(digits = 1, combo_needed)) %>% 
	mutate(one_in_two_needed = "yes")
