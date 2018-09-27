### RFU import for the Anc4 pilot
library(tidyverse)
library(cowplot)
library(readxl)
library(janitor)

RFU_files <- c(list.files("data-raw/anc4-pilot-september2018/anc4-pilot-2018-09-26", full.names = TRUE))

RFU_files <- RFU_files[grepl(".xlsx", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "")

all_plates <- map_df(RFU_files, read_excel, range = "B56:N64", .id = "file_name") %>% 
	rename(row = X__1)
all_times <- map_df(RFU_files, read_excel, range = "A6:B8", .id = "file_name") %>% 
	clean_names() %>% 
	filter(!is.na(plate_1)) %>% 
	select(-plate_2) %>% 
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
	
plate_layout <- read_excel("data-general/plate-layout-anc4-2018-09-26.xlsx") %>% 
	select(-population)
color_key <- read_excel("data-general/colour-key-anc4-2018-09-26.xlsx")

plate_key <- read_csv("data-general/anc4-pilot-plate-key.csv") %>% 
	rename(plate = plate_id) %>% 
	mutate(plate = as.character(plate))

plate2 <- left_join(plate_layout, color_key, by = "colour") %>% 
	mutate(well = str_to_upper(well)) %>% 
	filter(!is.na(population)) %>% 
	rename(Population = population)

treatments <- read_excel("data-general/ChlamEE_Treatments.xlsx") %>% 
	mutate(Population = as.character(Population))


all_treatments <- left_join(plate2, treatments, by = "Population") %>% 
	unite(col = unique_id, well, Population, Ancestor_ID, Treatment, sep = "_", remove = FALSE) %>% 
	mutate(well = str_to_upper(well))

all_rfus <- left_join(all_temp_RFU, all_treatments) %>% 
	separate(col = file_name, into = c("path", "plate"), sep = "plate", remove = FALSE)

all_rfus %>% 
	mutate(Treatment = ifelse(is.na(Treatment), "COMBO", Treatment)) %>% 
	ggplot(aes(x = Treatment, y = RFU, color = Treatment)) + geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
	facet_wrap( ~ plate)
ggsave("figures/anc4-pilot-RFU-2018-09-26.pdf", width = 14, height =8)


library(lubridate)
all_rfus2 <- all_rfus %>%
	unite(col = date_time, Date, time, sep = " ") %>%
	mutate(date_time = ymd_hms(date_time))

all_rfus3 <- left_join(all_rfus2, plate_key, by = "plate") %>% 
	mutate(temperature = as.numeric(temperature))

all_rfus3 %>% 	
filter(plate %in% c("4", "8", "12", "16", "20", "24")) %>% 
	filter(temperature > 20) %>% 
	mutate(Treatment = ifelse(is.na(Treatment), "COMBO", Treatment)) %>% 
	ggplot(aes(x = date_time, y = RFU, color = Treatment, group = well)) + geom_point(alpha = 0.5, size = 2) +
	geom_line() +
	facet_wrap( ~ temperature)
ggsave("figures/anc4-pilot-RFU-time.pdf", width = 8, height = 6)

all_rfus3 %>% 
	mutate(population_id = paste(plate, well, sep = "_")) %>% 
	filter(plate %in% c("4", "8", "12", "16", "20", "24")) %>% 
	filter(temperature > 20) %>% 
	mutate(Treatment = ifelse(is.na(Treatment), "COMBO", Treatment)) %>% 
	ggplot(aes(x = date_time, y = RFU, color = factor(temperature), group = population_id)) + geom_point(alpha = 0.5, size = 2) +
	geom_line() +
	facet_wrap( ~ Treatment) + scale_color_viridis_d()

all_rfus3 %>% 	
	filter(plate %in% c(13, 16)) %>% 
	mutate(Treatment = ifelse(is.na(Treatment), "COMBO", Treatment)) %>% 
	ggplot(aes(x = date_time, y = RFU, color = Treatment, group = well)) + geom_point(alpha = 0.5, size = 2) +
	geom_line() +
	facet_wrap( ~ Treatment + plate)
ggsave("figures/anc4-pilot-RFU-plate16-35degrees.pdf", width = 12, height = 6)
