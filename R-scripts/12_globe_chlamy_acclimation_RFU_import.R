### RFU import for GlobeChlamy - Acclimation trial
library(tidyverse)
library(cowplot)
library(readxl)
library(janitor)
library(lubridate)


# read in general data ----------------------------------------------------


plate_layout <- read_excel("data-general/globe-chlamy-acclimation-plate-setup.xlsx") 

plate_key <- read_excel("data-general/globe-chlamy-acclimation-plate-key.xlsx") %>% 
	mutate(temperature = ifelse(temperature == 20, 22, temperature))

plate_info <- left_join(plate_layout, plate_key, by = "plate_key")



# read in RFU data --------------------------------------------------------
RFU_files <- c(list.files("data-raw/globe-chlamy-RFU-acclim", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xlsx", RFU_files)]
RFU_files <- RFU_files[!grepl("acclimation", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "")

# RFU_files <- RFU_files[!grepl("acc", RFU_files)]

all_plates <- map_df(RFU_files, read_excel, range = "B50:E52", .id = "file_name") %>% 
	rename(row = X__1) %>% 
	filter(!grepl("setup", file_name)) 

all_times <- map_df(RFU_files, read_excel, range = "A6:B8", .id = "file_name") %>% 
	clean_names() %>% 
	filter(!is.na(plate_1)) %>% 
	filter(!grepl("setup", file_name)) %>% 
	spread(key = plate_number, value = plate_1) %>% 
	separate(Time, into = c("crap", "time"), sep = " ") %>% 
	select(-crap) %>% 
	separate(file_name, into = c("path", "plate"), sep = "plate", remove = FALSE)



all_plates2 <- left_join(all_plates, all_times, by = "file_name")

all_temp_RFU <- all_plates2 %>% 
	gather(key = column, value = RFU, 3:5) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU)) %>% 
	mutate(plate = as.numeric(plate))


plate2 <- plate_layout %>% 
	mutate(well = str_to_upper(well)) 


all_rfus_raw <- left_join(all_temp_RFU, plate_info, by = c("well", "plate")) %>% 
	mutate(plate = as.numeric(plate))

# all_rfus <- left_join(all_rfus_raw, plate_key, by = "plate")

all_rfus <- all_rfus_raw

all_rfus2 <- all_rfus %>%
	unite(col = date_time, Date, time, sep = " ") %>%
	mutate(date_time = ymd_hms(date_time)) %>% 
	mutate(start_time = min(date_time)) %>% 
	mutate(days = interval(start_time, date_time)/ddays(1)) %>% 
	unite(col = well_plate, well, plate, remove =  FALSE) 


all_rfus2 %>% 
	separate(well, into =c("row", "column"), sep = 1, remove = FALSE) %>% 
	filter(temperature == 10) %>% 
	ggplot(aes(x = days, y = RFU, group = well_plate, color = row)) + geom_point() +
	facet_wrap( ~ population + temperature, scales = "free") +
	geom_line() + scale_color_viridis_d(name = "Row")
ggsave("figures/acclimation-time-well.pdf", width = 16, height = 16)
