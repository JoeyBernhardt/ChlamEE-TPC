### RFU import for GlobeChlamy
library(tidyverse)
library(cowplot)
library(readxl)
library(janitor)
library(lubridate)


# read in general data ----------------------------------------------------


plate_layout <- read_excel("data-general/plate-layout-globe-chlamy.xlsx") 

plate_key <- read_excel("data-general/globe-chlamy-plate-key.xlsx") 





# read in RFU data --------------------------------------------------------
RFU_files <- c(list.files("data-raw/globe-chlamy-RFU", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xlsx", RFU_files)]
	
names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "")
	
# RFU_files <- RFU_files[!grepl("acc", RFU_files)]

all_plates <- map_df(RFU_files, read_excel, range = "B56:N64", .id = "file_name") %>% 
	rename(row = X__1) %>% 
	filter(!grepl("setup", file_name)) %>%
	filter(!grepl("acc", file_name))
	
all_times <- map_df(RFU_files, read_excel, range = "A6:B8", .id = "file_name") %>% 
	clean_names() %>% 
	filter(!grepl("acc", file_name)) %>% 
	filter(!is.na(plate_1)) %>% 
	filter(!grepl("setup", file_name)) %>% 
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

### find innoculation densities

inoc_densities <- all_rfus2 %>% 
	filter(date_time < ymd_hms("2018-10-11 08:24:40")) %>% 
	unite(col = well_plate, well, plate, remove =  FALSE) %>% 
	group_by(population) %>% 
	summarise_each(funs(mean), RFU) %>% 
	mutate(date_time = ymd_hms("2018-10-10 20:33:46")) %>% 
	mutate(round ="repeat")

	inoc10 <- inoc_densities %>% 
	mutate(temperature = 10) %>% 
	mutate(plate = 4) %>% 
		mutate(days = 0)
inoc16 <- inoc_densities %>% 
	mutate(temperature = 16)%>% 
	mutate(plate = 8)%>% 
	mutate(days = 0)
inoc22 <- inoc_densities %>% 
	mutate(temperature = 22) %>% 
	mutate(plate = 12)%>% 
	mutate(days = 0)
inoc28 <- inoc_densities %>% 
	mutate(temperature = 28) %>% 
	mutate(plate = 16)%>% 
	mutate(days = 0)
inoc34 <- inoc_densities %>% 
	mutate(temperature = 34) %>% 
	mutate(plate = 20)%>% 
	mutate(days = 0)
inoc40 <- inoc_densities %>% 
	mutate(temperature = 40) %>% 
	mutate(plate = 24)%>% 
	mutate(days = 0)
all_inocs <- bind_rows(inoc10, inoc16, inoc22, inoc28, inoc34, inoc40)

all_rfus3 <- bind_rows(all_rfus2, all_inocs) %>% 
	mutate(round = ifelse(plate %in% c(16, 20, 24, 12, 8, 4), "repeat", "single")) %>% 
	mutate(start_time = min(date_time)) %>% 
	mutate(days = interval(start_time, date_time)/ddays(1)) %>% 
	unite(col = well_plate, well, plate, remove =  FALSE) 


all_rfus3 %>% 
	dplyr::filter(temperature == 10) %>% 
	# filter(days < 3) %>% 
	# filter(population == 1) %>% 
	ggplot(aes(x = days, y = RFU, color = factor(plate), group = well_plate)) +
	geom_point(size = 2) + scale_color_viridis_d(name = "Temperature") + xlab("Days") +
	facet_wrap( ~ population, scales = "free") +
	geom_line() 
ggsave("figures/globe-chlamy-RFU-time-10C.pdf", width = 12, height = 10)

### which parts are exponential


all_rfus3 %>% 
	filter(temperature !=20) %>% 
	# filter(round == "repeat") %>% 
	# filter(population == 1) %>% 
	ggplot(aes(x = days, y = RFU, color = round, group = well_plate)) +
	geom_point(size = 2) + scale_color_viridis_d(name = "Plate") + xlab("Days") +
	facet_wrap( ~ population + temperature, scales = "free") +
	# geom_line() +
	xlim(0, 5)
ggsave("figures/globe-chlamy-RFU-time-round.pdf", width = 20, height = 15)



all4 <- all_rfus3 %>% 
	mutate(exponential = case_when(temperature == 34 & days < 1.75 ~ "yes",
								   temperature == 28 & days < 1.75 ~ "yes",
								   temperature == 22 & days < 2.5 ~ "yes",
								   temperature == 16 & days < 3 ~ "yes",
								   temperature == 10 & days < 20 ~ "yes",
								   temperature == 40 & days < 7 ~ "yes",
								   TRUE ~ "no"))

exponential <- all4 %>% 
	filter(exponential == "yes")

write_csv(exponential, "data-processed/globe-chlamy-exponential-RFU-time.csv")


all_rfus3 %>% 
	filter(round == "repeat") %>% 
	ggplot(aes(x = days, y = RFU, color = factor(population), group = well_plate)) +
	geom_point(size = 3) + scale_color_discrete(name = "Population") + xlab("Days") +
	facet_wrap( ~ temperature) + geom_line()


all_rfus2 %>% 
	ggplot(aes(x = date_time, y = RFU, color = factor(population), group = population)) +
	geom_point(size = 2) + scale_color_viridis_d(name = "Population") + xlab("Date") +
	facet_wrap( ~ temperature)

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
