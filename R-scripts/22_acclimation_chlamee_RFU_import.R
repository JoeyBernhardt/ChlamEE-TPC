### RFU import for ChlamEE
library(tidyverse)
library(cowplot)
library(readxl)
library(janitor)
library(lubridate)


# read in general data ----------------------------------------------------


plate_layout <- read_excel("data-general/ChlamEE-acclimation-plate_layout.xlsx") %>% 
	clean_names() %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "")

plate_key <- read_excel("data-general/chlamee-acclimation-plate-key.xlsx") %>% 
	mutate(plate = as.character(plate))

all_keys <- left_join(plate_key, plate_layout, by = "plate_key")




# read in RFU data --------------------------------------------------------
RFU_files <- c(list.files("data-raw/chlamee-acclimation", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xlsx", RFU_files)]

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
	filter(!is.na(RFU))


plate2 <- plate_layout %>% 
	mutate(well = str_to_upper(well)) 


all_rfus_raw <- left_join(all_temp_RFU, all_keys, by = c("plate", "well")) %>% 
	mutate(plate = as.numeric(plate))

# all_rfus <- left_join(all_rfus_raw, plate_key, by = "plate")

all_rfus <- all_rfus_raw

all_rfus2 <- all_rfus %>%
	unite(col = date_time, Date, time, sep = " ") %>%
	mutate(date_time = ymd_hms(date_time))


# plate2 <- plate_layout %>% 
# 	mutate(well = str_to_upper(well)) 
# 
# 
# all_rfus2 <- all_rfus %>%
# 	unite(col = date_time, Date, time, sep = " ") %>%
# 	mutate(date_time = ymd_hms(date_time))

all_rfus3 <- all_rfus2 %>% 
	# mutate(round = ifelse(plate %in% c(36, 30, 24, 18, 12, 6), "repeat", "single")) %>% 
	mutate(start_time = min(date_time)) %>% 
	mutate(days = interval(start_time, date_time)/ddays(1)) %>% 
	unite(col = well_plate, well, plate, remove =  FALSE) 

write_csv(all_rfus3, "data-processed/chlamee-acclimation-RFUs.csv")
all_rfus3 %>%
	# filter(!plate %in% c(37, 38, 39, 40)) %>% 
	# filter(population_id == "blank") %>% 
	# filter(temperature > 10, temperature < 40) %>% 
	# filter(population %in% c(1, 2, 5, 6)) %>% 
	filter(RFU < 1000) %>% 
	ggplot(aes(x = days, y = RFU, color = factor(temperature), group = well_plate)) +
	geom_point(size = 2) +
	scale_color_viridis_d(name = "Temperature") +
	xlab("Days") +
	facet_wrap( ~ population_id, scales = "free") +
	geom_line() 
ggsave("figures/acclimation-chlamee.pdf", width = 20, height = 20)


all_rfus3 %>%
filter(RFU < 1000) %>% 
	ggplot(aes(x = days, y = RFU, color = factor(temperature), group = well_plate)) +
	geom_point(size = 2) +
	scale_color_viridis_d(name = "Temperature") +
	xlab("Days") +
	facet_wrap( ~ temperature, scales = "free") +
	geom_line() 
