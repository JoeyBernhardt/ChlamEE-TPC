

### RFU import for GlobeChlamy
library(tidyverse)
library(cowplot)
library(readxl)
library(janitor)
library(lubridate)


# read in general data ----------------------------------------------------


plate_layout <- read_excel("data-general/plate-layout-globe-chlamy.xlsx") 

plate_key <- read_excel("data-general/globe-chlamy-acclimated-plate-key.xlsx") 





# read in RFU data --------------------------------------------------------
RFU_files <- c(list.files("data-raw/globe-chlamy-RFU-acclimated", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xlsx", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "")

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







all_plates2 <- left_join(all_plates, all_times, by = "file_name")

all_temp_RFU <- all_plates2 %>% 
	gather(key = column, value = RFU, 3:14) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU))


plate2 <- plate_layout %>% 
	mutate(well = str_to_upper(well)) 


all_rfus_raw <- left_join(all_temp_RFU, plate2, by = "well") %>% 
	mutate(plate = as.numeric(plate))

all_rfus <- left_join(all_rfus_raw, plate_key, by = "plate")


all_rfus2 <- all_rfus %>%
	unite(col = date_time, Date, time, sep = " ") %>%
	mutate(date_time = ymd_hms(date_time))

all_rfus3 <- all_rfus2 %>% 
	mutate(round = ifelse(plate %in% c(36, 30, 24, 18, 12, 6), "repeat", "single")) %>% 
	mutate(start_time = min(date_time)) %>% 
	mutate(days = interval(start_time, date_time)/ddays(1)) %>% 
	unite(col = well_plate, well, plate, remove =  FALSE) 


all_rfus3 %>%
	filter(!plate %in% c(37, 38, 39, 40)) %>% 
	# filter(round == "repeat") %>% 
	# filter(temperature == 22) %>% 
	# filter(population %in% c(1, 2, 5, 6)) %>% 
	# filter(days < 0.3) %>% 
	ggplot(aes(x = days, y = RFU, color = factor(temperature), group = well_plate)) +
	geom_point(size = 2) +
	scale_color_viridis_d(name = "Temperature") +
	xlab("Days") +
	facet_wrap( ~ population, scales = "free") +
	geom_line() 
ggsave("figures/globe-chlamy-acclimated-RFU-time.pdf", width = 10, height = 8)



all_rfus3 %>%
	filter(!plate %in% c(37, 38, 39, 40)) %>% 
	# filter(round == "repeat") %>% 
	filter(temperature ==10) %>% 
	# filter(population %in% c(1, 2, 5, 6)) %>% 
	# filter(days < 0.3) %>% 
	ggplot(aes(x = days, y = RFU, color = factor(temperature), group = well_plate)) +
	geom_point(size = 2) +
	scale_color_viridis_d(name = "Temperature") +
	xlab("Days") +
	facet_wrap( ~ population, scales = "free") +
	geom_line() 


all_rfus4 <- all_rfus3 %>%
	filter(!plate %in% c(37, 38, 39, 40)) 

write_csv(all_rfus4, "data-processed/globe-chlamy-acclimated-RFU-time.csv")


all4 <- all_rfus4 %>% 
	mutate(exponential = case_when(temperature == 34 & days < 1 ~ "yes",
								   temperature == 28 & days < 1 ~ "yes",
								   temperature == 22 & days < 2.5 ~ "yes",
								   temperature == 16 & days < 7 ~ "yes",
								   temperature == 10 & days < 20 & days > 1.5 ~ "yes",
								   temperature == 40 & days < 7 & days > 1.5 ~ "yes",
								   temperature == 22 & population %in% c(3, 4, 5, 6) ~ "no",
								   TRUE ~ "no"))

exponential <- all4 %>% 
	mutate(exponential = ifelse(temperature == 22 & population %in% c(3, 4, 5, 6), "no", exponential)) %>% 
	filter(exponential == "yes") 


write_csv(exponential, "data-processed/globe-acclimated-exponential.csv")
