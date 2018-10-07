


library(tidyverse)
library(cowplot)
library(janitor)
library(tidyverse)
library(readxl)
library(stringr)
library(flowCore)

#### read in plate layouts and treatments
plate_layout <- read_excel("data-general/plate-layout-anc4-2018-09-26.xlsx") %>% 
	select(-population)
color_key <- read_excel("data-general/colour-key-anc4-2018-09-26.xlsx")

plate2 <- left_join(plate_layout, color_key, by = "colour") %>% 
	mutate(well = str_to_upper(well)) %>% 
	dplyr::filter(!is.na(population)) %>% 
	rename(Population = population)

treatments <- read_excel("data-general/ChlamEE_Treatments.xlsx") %>% 
	mutate(Population = as.character(Population))


all_treatments <- left_join(plate2, treatments, by = "Population") %>% 
	unite(col = unique_id, well, Population, Ancestor_ID, Treatment, sep = "_", remove = FALSE) %>% 
	mutate(well = str_to_upper(well))



### let's read in all of the fcs files, plate 25
fcs_files_25 <- c(list.files("flow-cytometry-data/20181005_112556-ANC4-plate25", full.names = TRUE))

names(fcs_files_25) <- fcs_files_25 %>% 
	gsub(pattern = ".fcs$", replacement = "")


#### Step 3: read in all the files!


read_fcs <- function(file) {
	fcs_file <- as.data.frame((exprs(read.FCS(file)))) %>% 
		clean_names()
	
}
all_fcs_25 <- map_df(fcs_files_25, read_fcs, .id = "file_name")


all_fcs2_25 <- all_fcs_25 %>% 
	separate(file_name, into = c("file_path", "well"), sep = -3, remove = FALSE) %>% 
	mutate(plate = 25)


colfunc <- colorRampPalette(c("white", "lightblue", "green", "yellow", "red"))


all_fcs2_25 %>% 
	dplyr::filter(well == "B07", fl1_a > 0) %>%
	ggplot(aes(x=log(fsc_a), y=log(fl3_a))) +
	# ylim(0, 30000) +
	# xlim(0,50000) +
	stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
	scale_fill_gradientn(colours=colfunc(100)) + 
	geom_density2d(colour="black", bins=10) +
	facet_wrap( ~ well)

all_fcs3_25 <- all_fcs2_25 %>% 
	select(1:8) %>% 
	dplyr::filter(fl1_a > 5, fl3_a > 5) %>% 
	mutate(well = str_to_upper(well))

all_fcs3_25 %>% 
	dplyr::filter(grepl("B", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 160000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 160000, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 160000, color = "black")
ggsave("figures/anc4-fc-b-wells-fl1-3-gated-plate25.png", height = 8, width =12)

all_fcs3_25 %>% 
	dplyr::filter(grepl("C", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 160000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 160000, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 160000, color = "black")
ggsave("figures/anc4-fc-c-wells-fl1-3-gated-plate25.png", height = 8, width =12)

all_fcs3_25 %>% 
	dplyr::filter(grepl("D", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a > 1250, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 160000, color = "black")
ggsave("figures/anc4-fc-d-wells-fl1-3-gated-plate25.png", height = 8, width =12)


sorted_25 <- all_fcs3_25 %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a > 1250, "algae", type)) 

counts_25 <- sorted_25 %>% 
	mutate(type = ifelse(is.na(type), "background", type)) %>% 
	group_by(well, type) %>% 
	tally() %>% 
	dplyr::filter(type == "algae")



all_counts_25 <- left_join(counts_25, all_treatments) %>% 
	mutate(plate = 25)

### let's read in all of the fcs files, plate 5
fcs_files_5 <- c(list.files("flow-cytometry-data/20181005_122904-ANC4-plate5", full.names = TRUE))

names(fcs_files_5) <- fcs_files_5 %>% 
	gsub(pattern = ".fcs$", replacement = "")
all_fcs_5 <- map_df(fcs_files_5, read_fcs, .id = "file_name")


all_fcs2_5 <- all_fcs_5 %>% 
	separate(file_name, into = c("file_path", "well"), sep = -3, remove = FALSE) %>% 
	mutate(plate = 5)
all_fcs3_5 <- all_fcs2_5 %>% 
	select(1:8) %>% 
	dplyr::filter(fl1_a > 5, fl3_a > 5) %>% 
	mutate(well = str_to_upper(well))


all_fcs3_5 %>% 
	dplyr::filter(grepl("D", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	# mutate(type = ifelse(fl3_a < 160000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 1250, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) +
	# geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 1250, color = "black")

sorted_5 <- all_fcs3_5 %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a > 1250, "algae", type)) 

counts_5 <- sorted_5 %>% 
	mutate(type = ifelse(is.na(type), "background", type)) %>% 
	group_by(well, type) %>% 
	tally() %>% 
	dplyr::filter(type == "algae")

all_counts_5 <- left_join(counts_5, all_treatments) %>% 
	mutate(plate = 5)

### let's read in all of the fcs files, plate 6
fcs_files_6 <- c(list.files("flow-cytometry-data/20181005_155259-ANC4-plate6", full.names = TRUE))

names(fcs_files_6) <- fcs_files_6 %>% 
	gsub(pattern = ".fcs$", replacement = "")
all_fcs_6 <- map_df(fcs_files_6, read_fcs, .id = "file_name")


all_fcs2_6 <- all_fcs_6 %>% 
	separate(file_name, into = c("file_path", "well"), sep = -3, remove = FALSE) %>% 
	mutate(plate = 6)
all_fcs3_6 <- all_fcs2_6 %>% 
	select(1:8) %>% 
	dplyr::filter(fl1_a > 5, fl3_a > 5) %>% 
	mutate(well = str_to_upper(well))


all_fcs3_6 %>% 
	dplyr::filter(grepl("D", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	# mutate(type = ifelse(fl3_a < 160000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 1250, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) +
	# geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 1250, color = "black")

sorted_6 <- all_fcs3_6 %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a > 1250, "algae", type)) 

counts_6 <- sorted_6 %>% 
	mutate(type = ifelse(is.na(type), "background", type)) %>% 
	group_by(well, type) %>% 
	tally() %>% 
	dplyr::filter(type == "algae")

all_counts_6 <- left_join(counts_6, all_treatments) %>% 
	mutate(plate = 6)

all_plate_counts <- bind_rows(all_counts_25, all_counts_6, all_counts_5)


write_csv(all_plate_counts, "data-processed/all_plate_cell_counts_cytometry.csv")


fcs_files_all <- c(list.files("flow-cytometry-data/ANC4-experiment/",
							  full.names = TRUE, pattern = ".fcs$", recursive = TRUE))


names(fcs_files_all) <- fcs_files_all %>% 
	gsub(pattern = ".fcs$", replacement = "")

all_fcs_all <- map_df(fcs_files_all, read_fcs, .id = "file_name")

all_fcs_all$file_name[[1]]

all_fcs2_all <- all_fcs_all %>% 
	separate(file_name, into = c("file_path", "plate_well"), sep = c("/plate"), remove = FALSE) %>% 
	separate(plate_well, into = c("plate", "well"), sep = -3, remove = FALSE) %>% 
	mutate(plate = str_replace(plate, "_", "")) %>%
	mutate(plate = str_replace(plate, "_", ""))

all_fcs3_all <- all_fcs2_all%>% 
	select(4:10) %>% 
	dplyr::filter(fl1_a > 5, fl3_a > 5) %>% 
	mutate(well = str_to_upper(well))

sorted_all <- all_fcs3_all %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a > 1250, "algae", type)) 

counts_all <- sorted_all %>% 
	mutate(type = ifelse(is.na(type), "background", type)) %>% 
	group_by(plate, well, type) %>% 
	tally() %>% 
	dplyr::filter(type == "algae")

all_counts_all <- left_join(counts_all, all_treatments)
write_csv(all_counts_all, "data-processed/all_cell_counts_ANC4.csv")

