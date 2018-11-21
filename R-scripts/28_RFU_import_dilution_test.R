### dilution RFU import


library(tidyverse)
library(cowplot)
library(readxl)
library(janitor)
library(lubridate)


# read in general data ----------------------------------------------------


plate_layout <- read_excel("data-general/rfu-count-dilution-plate-layout.xlsx") 

# plate_key <- read_excel("data-general/chlamee-acclimated-plate-key.xlsx") 

plate_info <- plate_layout
	rename(column = colum) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>%
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	mutate(population = ifelse(population == "anc 2", "Anc 2", population))



# read in RFU data --------------------------------------------------------
RFU_files <- c(list.files("data-raw/count-rfu-dilution-test-RFU/", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xlsx", RFU_files)]
RFU_files <- RFU_files[grepl("high", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "")

# RFU_files[grepl("104", RFU_files)]

# RFU_files <- RFU_files[!grepl("acc", RFU_files)]

all_plates <- map_df(RFU_files, read_excel, range = "B56:N64", .id = "file_name") %>%
	rename(row = X__1) 

all_times <- map_df(RFU_files, read_excel, range = "A6:B8", .id = "file_name") %>% 
	clean_names() %>% 
	filter(!is.na(plate_1)) %>% 
	# filter(!grepl("dilution", file_name)) %>% 
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
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") 




all_rfus_raw <- left_join(all_temp_RFU, plate2, by = c("well")) 



all_rfus2 <- all_rfus_raw %>%
	unite(col = date_time, Date, time, sep = " ") %>%
	mutate(date_time = ymd_hms(date_time))

write_csv(all_rfus2, "data-processed/dilution-test-low-density-rfu.csv")
write_csv(all_rfus2, "data-processed/dilution-test-high-density-rfu.csv")
