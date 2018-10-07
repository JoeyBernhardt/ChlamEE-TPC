

library(tidyverse)

rfus <- read_csv("data-processed/all_exponential_rfus.csv")
cell_counts <- read_csv("data-processed/all_cell_counts_ANC4.csv") %>% 
	select(n, well, plate) %>% 
	unite(col = population_id, plate, well, remove = FALSE)

plate25 <- rfus %>% 
	dplyr::filter(plate == "25") 

singles <- rfus %>% 
	dplyr::filter(plate %in% c("5", "6", "9", "10", "13", "14")) %>% 
	dplyr::filter(date_time > ymd_hms("2018-09-27 07:45:41")) 

all_singles <- bind_rows(plate25, singles)

all_rfu_counts <- left_join(all_singles, cell_counts, by = c("population_id", "well", "plate"))

all_rfu_counts2 <- all_rfu_counts %>% 
	mutate(cell_count = ifelse(is.na(n), 0, n)) %>% 
	mutate(cells_per_ml = cell_count*4*10) %>% 
	dplyr::filter(date_time != ymd_hms("2018-09-28 11:22:42"))

all_rfu_counts2 %>% 
	# mutate(temperature = ifelse(is.na(temperature), 30, temperature)) %>% 
	ggplot(aes(x = date_time, y = cell_count, group = population_id, color = factor(temperature))) + geom_point() +
	facet_wrap( ~ Treatment) + 
	# scale_color_viridis_d(name = "Temperature") +
	xlab("Date") + geom_line()

all_rfu_counts2 %>% 
	# dplyr::filter(Treatment == "B") %>% 
	mutate(temperature = ifelse(is.na(temperature), 30, temperature)) %>% 
	ggplot(aes(x = cells_per_ml, y = RFU, color = factor(Treatment))) + geom_point() +
	geom_smooth(method = "lm")


