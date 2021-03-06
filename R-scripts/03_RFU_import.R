### RFU import for the Anc4 pilot
library(tidyverse)
library(cowplot)
library(readxl)
library(janitor)
library(lubridate)


# read in general data ----------------------------------------------------



plate_layout <- read_excel("data-general/plate-layout-anc4-2018-09-26.xlsx") %>% 
	select(-population)
color_key <- read_excel("data-general/colour-key-anc4-2018-09-26.xlsx")

plate_key <- read_csv("data-general/anc4-pilot-plate-key.csv") %>% 
	rename(plate = plate_id) %>% 
	mutate(plate = as.character(plate))

treatments <- read_excel("data-general/ChlamEE_Treatments.xlsx") %>% 
	mutate(Population = as.character(Population))

 


# read in RFU data --------------------------------------------------------
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
	

plate2 <- left_join(plate_layout, color_key, by = "colour") %>% 
	mutate(well = str_to_upper(well)) %>% 
	filter(!is.na(population)) %>% 
	rename(Population = population)


all_treatments <- left_join(plate2, treatments, by = "Population") %>% 
	unite(col = unique_id, well, Population, Ancestor_ID, Treatment, sep = "_", remove = FALSE) %>% 
	mutate(well = str_to_upper(well))

all_rfus <- left_join(all_temp_RFU, all_treatments) %>% 
	separate(col = file_name, into = c("path", "plate"), sep = "plate", remove = FALSE)

# all_rfus %>% 
# 	mutate(Treatment = ifelse(is.na(Treatment), "COMBO", Treatment)) %>% 
# 	ggplot(aes(x = Treatment, y = RFU, color = Treatment)) + geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
# 	facet_wrap( ~ plate)
# ggsave("figures/anc4-pilot-RFU-2018-09-26.pdf", width = 14, height =8)


all_rfus2 <- all_rfus %>%
	unite(col = date_time, Date, time, sep = " ") %>%
	mutate(date_time = ymd_hms(date_time))

all_rfus3 <- left_join(all_rfus2, plate_key, by = "plate") %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(Treatment = ifelse(colour == "blank", "COMBO", Treatment)) %>% 
	mutate(Treatment = ifelse(colour == "wx", "Anc4", Treatment)) %>% 
	mutate(population_id = paste(plate, well, sep = "_")) %>% 
	# filter(plate %in% c("4", "8", "12", "16", "20", "24")) %>% 
	mutate(round = ifelse(plate %in% c("4", "8", "12", "16", "20", "24"), "repeat", "single"))
	# filter(!is.na(temperature)) 



# Plot! -------------------------------------------------------------------

## RFU over time
all_rfus3 %>% 
	filter(plate %in% c("4", "8", "12", "16", "20", "24")) %>% 
	filter(!Treatment %in% c("COMBO", "P")) %>% 
	filter(temperature == 8) %>% 
	ggplot(aes(x = date_time, y = RFU, color = factor(temperature), group = population_id)) + geom_point(size = 2) +
	geom_line() +
	facet_wrap( ~ Treatment, scales = "free") + scale_color_viridis_d(name = "Temperature") + xlab("Date")
ggsave("figures/anc4-pilot-RFU-time.pdf", width = 12, height = 10)


single_plates <- all_rfus3 %>% 
	filter(!plate %in% c("4", "8", "12", "16", "20", "24")) %>% 
	distinct(plate, date_time, .keep_all = TRUE) 

write_csv(single_plates, "data-processed/single-plates.csv")

all_rfus3 %>% 
	# filter(plate %in% c("4", "8", "12", "16", "20", "24")) %>% 
	dplyr::filter(temperature == 30) %>% 
	# dplyr::filter(plate == 6) %>% View
	ggplot(aes(x = date_time, y = RFU, color = plate, group = population_id)) + geom_point(size = 2) +
	geom_line() +
	facet_wrap( ~ Treatment) + scale_color_viridis_d(name = "Temperature") + xlab("Date")



repeat_rfus <- all_rfus3 %>% 
	filter(plate %in% c("4", "8", "12", "16", "20", "24")) %>% 
	mutate(keep = NA) %>% 
	mutate(keep = case_when(temperature == 35 & date_time < ymd_hms("2018-09-28 18:05:06") ~ "yes",
							temperature == 30 & date_time < ymd_hms("2018-09-28 18:05:06") ~ "yes",
							temperature == 25 & date_time < ymd_hms("2018-09-28 18:05:06") ~ "yes",
							temperature == 20 & date_time < ymd_hms("2018-09-30 09:38:11") ~ "yes",
							temperature == 12 ~ "yes",
							temperature == 8 ~ "yes",
							TRUE ~ "no"))
exponential_repeats <- repeat_rfus %>% 
	filter(keep == "yes")

write_csv(exponential_repeats, "data-processed/exponential_repeats_RFU_anc4.csv")


all_repeat_rfus <- all_rfus3 %>% 
	# filter(plate %in% c("4", "8", "12", "16", "20", "24")) %>% 
	mutate(keep = NA) %>% 
	mutate(keep = case_when(temperature == 35 & date_time < ymd_hms("2018-09-28 18:05:06") ~ "yes",
							temperature == 30 & date_time < ymd_hms("2018-09-28 18:05:06") ~ "yes",
							temperature == 25 & date_time < ymd_hms("2018-09-28 18:05:06") ~ "yes",
							temperature == 20 & date_time < ymd_hms("2018-09-30 09:38:11") ~ "yes",
							temperature == 12 ~ "yes",
							temperature == 8 ~ "yes",
							is.na(temperature) ~ "yes",
							TRUE ~ "no"))
exponential_all <- all_repeat_rfus %>% 
	dplyr::filter(keep == "yes") %>% 
	dplyr::filter(plate == 6) %>% 
	mutate(inoculation = ifelse(date_time < ymd_hms("2018-09-27 04:05:06"), "yes", "no")) %>% View
	mutate(state = case_when(plate == 25 ~ "time0",
							 TRUE ~ "growth"))

write_csv(exponential_all, "data-processed/all_exponential_rfus.csv")
	
exponential_repeats %>%
	ggplot(aes(x = date_time, y = RFU, color = factor(temperature), group = population_id)) + geom_point(size = 2) +
	geom_line() +
facet_wrap( ~ Treatment) + scale_color_viridis_d(name = "Temperature") + xlab("Date")
	

all_rfus3 %>% 
	filter(temperature == 35) %>% 
	filter(date_time < ymd_hms("2018-09-28 18:05:06")) %>% 
	filter(plate %in% c("4", "8", "12", "16", "20", "24")) %>% 
	ggplot(aes(x = date_time, y = RFU, color = factor(temperature), group = population_id)) + geom_point(size = 2) +
	# geom_line() +
	facet_wrap( ~ temperature + Treatment) + scale_color_viridis_d(name = "Temperature") + xlab("Date")


all_rfus3 %>% 
	mutate(population_id = paste(plate, well, sep = "_")) %>% 
	# filter(!plate %in% c("4", "8", "12", "16", "20", "24")) %>% 
	mutate(round = ifelse(plate %in% c("4", "8", "12", "16", "20", "24"), "repeat", "single")) %>% 
	filter(!is.na(temperature)) %>% 
	# filter(round == "single") %>% 
	filter(temperature > 16) %>% 
	ggplot(aes(x = date_time, y = RFU, color = round, group = population_id)) + geom_point(size = 2) +
	geom_line() +
	facet_wrap( ~ temperature + Treatment, scales = "free") + scale_color_viridis_d(name = "Temperature") + xlab("Date")
ggsave("figures/anc4-pilot-RFU-time-single-repeat.pdf", width = 12, height = 10)

all_rfus3 %>% 	
filter(plate %in% c("4", "8", "12", "16", "20", "24")) %>% 
	filter(!is.na(temperature)) %>% 
	filter(temperature > 16) %>% 
	ggplot(aes(x = date_time, y = RFU, color = Treatment, group = well)) + geom_point(alpha = 0.5, size = 2) +
	geom_line() +
	facet_wrap( ~ temperature)
ggsave("figures/anc4-pilot-RFU-time.pdf", width = 10, height = 6)



all_rfus3 %>% 
	mutate(population_id = paste(plate, well, sep = "_")) %>% 
	# filter(plate %in% c("4", "8", "12", "16", "20", "24")) %>% 
	mutate(round = ifelse(plate %in% c("4", "8", "12", "16", "20", "24"), "repeat", "single")) %>% 
	filter(!is.na(temperature)) %>% 
	filter(temperature > 15) %>% 
	ggplot(aes(x = date_time, y = RFU, color = round, group = population_id)) + geom_point(size = 2) +
	geom_line() +
	facet_wrap( ~ Treatment + temperature) + scale_color_viridis_d(name = "Temperature") + xlab("Date")
ggsave("figures/anc4-pilot-RFU-repeat-single.pdf", width = 15, height = 15)


all_rfus3 %>% 	
	# filter(plate %in% c(13, 16)) %>% 
	mutate(Treatment = ifelse(is.na(Treatment), "COMBO", Treatment)) %>% 
	ggplot(aes(x = date_time, y = RFU, color = Treatment, group = well)) + geom_point(alpha = 0.5, size = 2) +
	geom_line() +
	facet_wrap( ~ Treatment + plate)
ggsave("figures/anc4-pilot-RFU-plate16-35degrees.pdf", width = 12, height = 6)
