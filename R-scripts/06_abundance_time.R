

library(tidyverse)
library(lubridate)
library(cowplot)

rfus <- read_csv("data-processed/all_exponential_rfus.csv")
cell_counts <- read_csv("data-processed/all_cell_counts_ANC4.csv") %>% 
	select(n, well, plate) %>% 
	unite(col = population_id, plate, well, remove = FALSE)

plate25 <- rfus %>% 
	dplyr::filter(plate == "25") 

singles <- rfus %>% 
	dplyr::filter(plate %in% c("5", "6", "9", "10", "13", "14", "1", "2")) %>% 
	dplyr::filter(date_time > ymd_hms("2018-09-27 07:45:41")) 

all_singles <- bind_rows(plate25, singles)

all_rfu_counts <- left_join(all_singles, cell_counts, by = c("population_id", "well", "plate"))

all_rfu_counts2 <- all_rfu_counts %>% 
	mutate(cell_count = ifelse(is.na(n), 0, n)) %>% 
	mutate(cells_per_ml = cell_count*4*10) %>% 
	dplyr::filter(date_time != ymd_hms("2018-09-28 11:22:42"))

all_rfu_counts2 %>% 
	mutate(temperature = ifelse(is.na(temperature), 20, temperature)) %>% 
	dplyr::filter(temperature == 20) %>% 
	ggplot(aes(x = date_time, y = cells_per_ml, group = population_id, color = factor(temperature))) + geom_point() +
	facet_wrap( ~ Treatment) + 
	# scale_color_viridis_d(name = "Temperature") +
	xlab("Date") + geom_line()

all_20 <- all_rfu_counts2 %>% 
	mutate(temperature = ifelse(is.na(temperature), 20, temperature)) %>% 
	dplyr::filter(temperature == 20) %>% 
	ungroup() %>% 
	select(Treatment, well, RFU, cells_per_ml, date_time) %>% 
	gather(key = type, value = count, RFU, cells_per_ml) %>% 
	ggplot(aes(x = date_time, y = count, group = well, color = type)) + geom_point()+
	geom_line() + facet_wrap( ~ Treatment + type, scales = "free", ncol = 4)
ggsave("figures/cells_and_RFU_20C.pdf", width = 10, height = 12)

all_35 <- all_rfu_counts2 %>% 
	mutate(temperature = ifelse(is.na(temperature), 35, temperature)) %>% 
	dplyr::filter(temperature == 35) %>% 
	ungroup() %>% 
	select(Treatment, well, RFU, cells_per_ml, date_time) %>% 
	gather(key = type, value = count, RFU, cells_per_ml) %>% 
	ggplot(aes(x = date_time, y = count, group = well, color = type)) + geom_point()+
	geom_line() + facet_wrap( ~ Treatment + type, scales = "free", ncol = 4)
ggsave("figures/cells_and_RFU_35C.pdf", width = 10, height = 12)

all_30 <- all_rfu_counts2 %>% 
mutate(temperature = ifelse(is.na(temperature), 30, temperature)) %>% 
	dplyr::filter(temperature == 30) %>% 
	ungroup() %>% 
	select(Treatment, well, RFU, cells_per_ml, date_time) %>% 
	gather(key = type, value = count, RFU, cells_per_ml) %>% 
	ggplot(aes(x = date_time, y = count, group = well, color = type)) + geom_point()+
	geom_line() + facet_wrap( ~ Treatment + type, scales = "free", ncol = 4)
ggsave("figures/cells_and_RFU_30C.pdf", width = 10, height = 12)

all_25 <- all_rfu_counts2 %>% 
	mutate(temperature = ifelse(is.na(temperature), 25, temperature)) %>% 
	dplyr::filter(temperature == 25) %>% 
	ungroup() %>% 
	select(Treatment, well, RFU, cells_per_ml, date_time) %>% 
	gather(key = type, value = count, RFU, cells_per_ml) %>% 
	ggplot(aes(x = date_time, y = count, group = well, color = well)) + geom_point()+
	geom_line() + facet_wrap( ~ Treatment + type, scales = "free", ncol = 4)
ggsave("figures/cells_and_RFU_25C_wells.pdf", width = 10, height = 12)


all_rfu_counts2 %>% 
	mutate(temperature = ifelse(is.na(temperature), 25, temperature)) %>% 
	dplyr::filter(temperature == 25) %>% 
	ggplot(aes(x = cells_per_ml, y = RFU)) + geom_point()+
 facet_wrap( ~ Treatment, scales = "free", ncol = 4)


rfus %>% 
	mutate(temperature = ifelse(is.na(temperature), 25, temperature)) %>% 
	dplyr::filter(temperature == 25) %>% 
	dplyr::filter(!plate %in% c("4", "8", "12", "16", "20", "24")) %>% 
	ungroup() %>% 
	ggplot(aes(x = date_time, y = RFU, group = well)) + geom_point()+
	geom_line() + facet_wrap( ~ Treatment, scales = "free", ncol = 4)
ggsave("figures/cells_and_RFU_25C.pdf", width = 10, height = 12)


all_temps <- bind_rows(all_20, all_25, all_30, all_35)

all_rfu_counts2 %>% 
	# dplyr::filter(Treatment == "B") %>% 
	mutate(temperature = ifelse(is.na(temperature), 30, temperature)) %>% 
	ggplot(aes(x = cells_per_ml, y = RFU, color = factor(Treatment))) + geom_point() +
	geom_smooth(method = "lm")


all_temps %>% 
	ggplot(aes(x = date_time, y = cells_per_ml, group = well, color = Treatment)) + geom_point() +
	geom_line() + facet_wrap( ~ temperature, scales = "free")
ggsave("figures/cells_v_time.pdf", width = 10, height = 12)

all_temps %>% 
	ggplot(aes(x = date_time, y = RFU, group = well, color = Treatment)) + geom_point() +
	geom_line() + facet_wrap( ~ temperature, scales = "free")
ggsave("figures/RFU_v_time.pdf", width = 10, height = 12)

write_csv(all_temps, "data-processed/all_temps_cells_RFU.csv")

all_temps <- read_csv("data-processed/all_temps_cells_RFU.csv")

all_temps %>% 
	ggplot(aes(x = cells_per_ml, y = RFU, color = factor(temperature))) + geom_point() +
	facet_wrap( ~ Treatment)
ggsave("figures/RFU_v_cells.pdf", width = 10, height = 12)
