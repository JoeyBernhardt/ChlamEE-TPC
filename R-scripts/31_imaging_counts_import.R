####
### chlamee-imaging import

library(tidyverse)
library(cowplot)
library(readxl)
library(janitor)
library(lubridate)


# read in general data ----------------------------------------------------


plate_layout <- read_excel("data-general/chlamee-acclimated-plate-layout.xlsx") 

plate_key <- read_excel("data-general/chlamee-acute-plate-key.xlsx") 

plate_info <- left_join(plate_layout, plate_key, by = c("plate_key")) %>% 
	rename(column = colum) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>%
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	mutate(population = ifelse(population == "anc 2", "Anc 2", population)) %>% 
	rename(plate = plate_number)


write_csv(plate_info, "data-processed/chlamee-acute-plate-info.csv")


# read in RFU data --------------------------------------------------------
RFU_files <- c(list.files("data-raw/chlamee-acute-imaging-spreadsheets", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xls", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xls$", replacement = "")

# RFU_files[grepl("104", RFU_files)]

# RFU_files <- RFU_files[!grepl("acc", RFU_files)]

all_plates <- map_df(RFU_files, read_excel, range = "B50:N58", .id = "file_name") %>% 
	rename(row = X__1) %>% 
	filter(!grepl("dilution", file_name)) %>% 
	separate(file_name, into = c("path", "plate"), sep = "plate", remove = FALSE) %>% 
	separate(plate, into = c("plate", "other"), sep = "_") %>% 
	select(-other) 

all_times <- map_df(RFU_files, read_excel, range = "A6:B8", .id = "file_name") %>% 
	clean_names() %>% 
	filter(!is.na(plate_1)) %>% 
	filter(!grepl("dilution", file_name)) %>% 
	spread(key = plate_number, value = plate_1) %>% 
	separate(Time, into = c("crap", "time"), sep = " ") %>% 
	select(-crap) %>% 
	separate(file_name, into = c("path", "plate"), sep = "plate", remove = FALSE) %>% 
	separate(plate, into = c("plate", "other"), sep = "_") %>% 
	select(-other)


all_temp_RFU <- all_plates %>% 
	gather(key = column, value = RFU, 5:16) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>%
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU)) %>% 
	mutate(plate = as.numeric(plate))


all_temp_RFU %>% 
	select(-path) %>%
	mutate(round = ifelse(grepl("04d", file_name), "d", "b")) %>% 
	select(-file_name) %>% 
	spread(key = round, value = RFU) %>%
	# filter(d < 100) %>% 
	ggplot(aes(x = d, y = b)) + geom_point() +
	geom_abline(slope = 1, intercept = 0)



all_rfus_raw <- left_join(all_temp_RFU, plate_info, by = c("well", "plate")) %>% 
	mutate(plate = as.numeric(plate)) %>% 
	filter(!is.na(plate_key))



all_rfus2 <- all_rfus_raw %>%
	unite(col = date_time, Date, time, sep = " ") %>%
	mutate(date_time = ymd_hms(date_time))

all_rfus3 <- all_rfus2 %>% 
	mutate(round = ifelse(plate %in% c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 76, 81, 86), "repeat", "single")) %>% 
	group_by(temperature) %>% 
	mutate(start_time = min(date_time)) %>% 
	mutate(days = interval(start_time, date_time)/ddays(1)) %>% 
	unite(col = well_plate, well, plate, remove =  FALSE) 


