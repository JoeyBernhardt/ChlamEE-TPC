library(tidyverse)


all_fcs2_all <- read_csv("data-processed/globe-chlamy-particles.csv")
plate_key <- read_excel("data-general/globe-chlamy-plate-key.xlsx")
plate_layout <- read_excel("data-general/plate-layout-globe-chlamy.xlsx")

plate2 <- plate_layout %>% 
	mutate(well = str_to_upper(well)) 

all_fcs3_all <- all_fcs2_all %>% 
	select(4:10) %>% 
	dplyr::filter(fl1_a > 5, fl3_a > 5) %>% 
	mutate(well = str_to_upper(well))

sorted_all <- all_fcs3_all %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a > 1250, "algae", type)) %>% 
	mutate(type = ifelse(is.na(type), "background", type))

all_sorted_all <- left_join(sorted_all, plate2, by = "well") %>% 
	# dplyr::filter(!grepl("A", well)) %>% 
	mutate(plate = as.integer(plate))

str(all_sorted_all)

all_sorted2 <- left_join(all_sorted_all, plate_key, by = "plate")

all_algae <- all_sorted2 %>% 
	dplyr::filter(type == "algae")


counts <- all_algae %>% 
	group_by(plate, well, temperature) %>% 
	tally() %>% 
	mutate(cells_per_ml = n*40)

### lets see how much carry over there is

counts %>% 
	separate(well, into = c("row", "column"), sep = 1, remove = FALSE) %>% 
	# dplyr::filter(row == "A") %>% 
	ggplot(aes(x = temperature, y = cells_per_ml, color = row)) + geom_point() +
	facet_wrap( ~ row, scales = "free")



plate25 <- counts %>% 
	dplyr::filter(plate == 25) 

n <- c(1, 1)
well <- c("C08", "E02")

ones <- as.data.frame(n, well) %>% 
	rownames_to_column() %>% 
	rename(well = rowname)

plate25b <- bind_rows(plate25, ones) %>% 
	ungroup() %>% 
	mutate(plate = 25) %>% 
	mutate(temperature = 20) %>% 
	mutate(cells_per_ml = ifelse(is.na(cells_per_ml), 40, cells_per_ml))

inoc16 <- plate25b %>% 
	mutate(date_time = ymd_hms("2018-10-10 20:33:46")) %>% 
	mutate(temperature = 16)
inoc10 <- plate25b %>% 
	mutate(date_time = ymd_hms("2018-10-10 20:33:46")) %>% 
	mutate(temperature = 10)

inoc22 <- plate25b %>% 
	mutate(date_time = ymd_hms("2018-10-10 20:33:46")) %>% 
	mutate(temperature = 22)
inoc28 <- plate25b %>% 
	mutate(date_time = ymd_hms("2018-10-10 20:33:46")) %>% 
	mutate(temperature = 28)
inoc34 <- plate25b %>% 
	mutate(date_time = ymd_hms("2018-10-10 20:33:46")) %>% 
	mutate(temperature = 34)
inoc40 <- plate25b %>% 
	mutate(date_time = ymd_hms("2018-10-10 20:33:46")) %>% 
	mutate(temperature = 40)

all_inocs <- bind_rows(inoc22, inoc28, inoc34, inoc10, inoc16, inoc40) 

all_inocs2 <- left_join(all_inocs, plate2, by = "well") %>% 
	dplyr::filter(!grepl("A", well))


# bring in RFU data -------------------------------------------------------

RFUS <- read_csv("data-processed/globe-chlamy-exponential-RFU-time.csv")

rfu2 <- RFUS %>% 
	dplyr::filter(date_time > ymd_hms("2018-10-11 18:00:00")) 

rfu2 %>% 
	ggplot(aes(x = date_time, y = RFU)) + geom_point() + facet_wrap( ~ plate)


all <- left_join(counts, rfu2, by = c("plate", "well", "temperature"))

all2 <- bind_rows(all, all_inocs2) 

all2 %>% 
	dplyr::filter(plate == 25) %>% View

all2 %>% 
	dplyr::filter(!is.na(population)) %>% 
	# filter(population == 4) %>% 
	ggplot(aes(x = date_time, y = log(cells_per_ml), group = temperature, color = factor(temperature))) + geom_point() +
	facet_wrap( ~ population) + geom_smooth(method = "lm")

all2 %>% 
	ungroup() %>% 
	distinct(temperature, plate) %>% View

all2 %>% 
	dplyr::filter(!is.na(population)) %>% 
	dplyr::filter(temperature == 10) %>% 
	ggplot(aes(x = date_time, y = cells_per_ml, group = temperature, color = factor(plate))) + geom_point() +
	facet_wrap( ~ population, scales = "free") 
ggsave("figures/cells_time.pdf", width = 16, height = 16)


all2 %>% 
	ggplot(aes(x = cells_per_ml, y = RFU, color = factor(temperature))) + geom_point() +
	geom_smooth()
