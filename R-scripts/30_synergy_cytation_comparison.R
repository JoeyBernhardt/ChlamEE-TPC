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


write_csv(plate_info, "data-processed/chlamee-acclimated-plate-info.csv")


# read in RFU data --------------------------------------------------------
RFU_files <- c(list.files("data-raw/synergy-cytation-comparisons/direct_comparisons", full.names = TRUE))
RFU_files <- RFU_files[grepl(".xls", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xls$", replacement = "")

# RFU_files[grepl("104", RFU_files)]

# RFU_files <- RFU_files[!grepl("acc", RFU_files)]

all_plates <- map_df(RFU_files, read_excel, range = "B56:N64", .id = "file_name") %>%
	rename(row = X__1) %>% 
	filter(!grepl("dilution", file_name)) %>% 
	gather(key = column, value = RFU, 3:14) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>%
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU)) %>% 
	separate(file_name, into = c("path", "machine"), sep = "direct_comparisons/", remove = FALSE) %>% 
	separate(machine, into = c("instrument", "run"), sep = "-", remove = FALSE) 


all_plates %>% 
	select(instrument, well, RFU) %>% 
	filter(RFU < 200) %>% 
	spread(key = instrument, value = RFU) %>% 
	ggplot(aes(x = cytation, y = synergy)) + geom_point() +
	geom_abline(slope = 1, intercept = 0) +
	geom_smooth(method = "lm")

	
all_plates %>% 
	select(instrument, well, RFU) %>% 
	filter(RFU < 200) %>% 
	spread(key = instrument, value = RFU) %>%
	lm(synergy ~ cytation, data =.) %>% summary()

all_plates$file_name[[1]]	
	
