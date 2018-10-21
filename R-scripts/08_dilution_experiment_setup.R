# packages

library(tidyverse)
library(cowplot)
library(readxl)
library(janitor)
library(lubridate)


# read in general data ----------------------------------------------------


plate_layout <- read_excel("data-general/globe-chlamy-dilution-plate.xlsx") %>% 
	mutate(well = str_to_upper(well)) 
RFU_files <- c(list.files("data-raw/globe-chlamy-RFU", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xlsx", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "")

all_plates <- map_df(RFU_files, read_excel, range = "B56:N64", .id = "file_name") %>% 
	rename(row = X__1)
all_times <- map_df(RFU_files, read_excel, range = "A6:B8", .id = "file_name") %>% 
	clean_names() %>% 
	filter(!is.na(plate_1)) %>% 
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
	filter(!is.na(RFU))

all_rfu <- left_join(all_temp_RFU, plate_layout, by = "well")  %>% 
	filter(!grepl("setup-", file_name))

write_csv(all_rfu, "data-processed/RFU-globe-chlamy-setup-b-2018-10-10.csv")

dilutions <- all_rfu %>% 
	group_by(population) %>% 
	summarise(mean_rfu = mean(RFU)) %>% 
	mutate(dilution_full = mean_rfu/11) %>% 
	mutate(dilution2 = 1/dilution_full)

high_conc <- dilutions %>% 
	# filter(dilution2 > 0.15) %>% 
	mutate(volume_needed = dilution2*30) %>% 
	mutate(combo_needed = 30-volume_needed) %>% 
	mutate(volume_needed = round(digits = 1, volume_needed)) %>% 
	mutate(combo_needed = round(digits = 1, combo_needed)) %>% 
	mutate(over_5 = ifelse(volume_needed > 5, "yes", "no")) %>% 
	mutate(over_10 = ifelse(volume_needed > 10, "yes", "no")) %>%
	mutate(over_5_divided_x2 = ifelse(over_5 == "yes", volume_needed/2, NA)) %>% 
	mutate(over_10_divided_x3 = ifelse(over_10 == "yes", volume_needed/5, NA)) %>% 
	mutate(combo_needed_divided_x3 = round(digits = 1, combo_needed/3)) 
	

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
	


write_csv(high_conc, "data-general/dilutions-chlamy-globe.csv")




# read in general data ----------------------------------------------------


plate_layout <- read_excel("data-general/globe-chlamy-acclimation-dilution-plate.xlsx") %>% 
	mutate(well = str_to_upper(well)) 
RFU_files <- c(list.files("data-raw/globe-chlamy-RFU", full.names = TRUE))
RFU_files <- "data-raw/globe-chlamy-RFU/globe-2018-10-18-acclimation-plate.xlsx"
RFU_files <- RFU_files[grepl(".xlsx", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "")

all_plates <- map_df(RFU_files, read_excel, range = "B56:N64", .id = "file_name") %>% 
	rename(row = X__1)
all_times <- map_df(RFU_files, read_excel, range = "A6:B8", .id = "file_name") %>% 
	clean_names() %>% 
	filter(!is.na(plate_1)) %>% 
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
	filter(!is.na(RFU))

all_rfu <- left_join(all_temp_RFU, plate_layout, by = "well")  %>% 
	filter(!grepl("setup-", file_name))

write_csv(all_rfu, "data-processed/RFU-globe-chlamy-setup-acclimation-2018-10-18.csv")

dilutions <- all_rfu %>% 
	group_by(population) %>% 
	summarise(mean_rfu = mean(RFU)) %>% 
	mutate(dilution_full = mean_rfu/10) %>% 
	mutate(dilution2 = 1/dilution_full) %>%
	mutate(volume_needed = dilution2*10*1000) 
	
	

high_conc <- dilutions %>% 
	# filter(dilution2 > 0.15) %>% 
	mutate(volume_needed = dilution2*30) %>% 
	mutate(combo_needed = 30-volume_needed) %>% 
	mutate(volume_needed = round(digits = 1, volume_needed)) %>% 
	mutate(combo_needed = round(digits = 1, combo_needed)) %>% 
	mutate(over_5 = ifelse(volume_needed > 5, "yes", "no")) %>% 
	mutate(over_10 = ifelse(volume_needed > 10, "yes", "no")) %>%
	mutate(over_5_divided_x2 = ifelse(over_5 == "yes", volume_needed/2, NA)) %>% 
	mutate(over_10_divided_x3 = ifelse(over_10 == "yes", volume_needed/5, NA)) %>% 
	mutate(combo_needed_divided_x3 = round(digits = 1, combo_needed/3)) 


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



write_csv(high_conc, "data-general/dilutions-chlamy-globe.csv")

6/13

