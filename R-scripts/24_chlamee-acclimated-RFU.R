



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





# read in RFU data --------------------------------------------------------
RFU_files <- c(list.files("data-raw/chlamee-RFU-acclimated", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xlsx", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "")

# RFU_files[grepl("44", RFU_files)]

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
	mutate(round = ifelse(plate %in% c(6, 13, 20, 27, 29, 36, 48, 55, 62, 69, 76, 83), "repeat", "single")) %>% 
	group_by(temperature) %>% 
	mutate(start_time = min(date_time)) %>% 
	mutate(days = interval(start_time, date_time)/ddays(1)) %>% 
	unite(col = well_plate, well, plate, remove =  FALSE) 


all_rfus3 %>%
	# filter(round == "repeat") %>% 
	filter(temperature %in% c(28)) %>% 
	ggplot(aes(x = days, y = RFU, color = factor(temperature), group = well_plate)) +
	geom_point(size = 2) +
	scale_color_viridis_d(name = "Temperature") +
	xlab("Days") +
	facet_wrap( ~ population, scales = "free_y") +
	geom_line(linetype = "dotted") 
ggsave("figures/chlamee-acclimated-RFU-time.pdf", width = 30, height = 20)
ggsave("figures/chlamee-acclimated-RFU-time-34C.pdf", width = 14, height = 8)
ggsave("figures/chlamee-acclimated-RFU-time-40C.pdf", width = 14, height = 8)
