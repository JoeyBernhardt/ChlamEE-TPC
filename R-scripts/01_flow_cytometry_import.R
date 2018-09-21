
library(flowCore)
library(tidyverse)
library(cowplot)
library(janitor)

### let's read in all of the fcs files
fcs_files <- c(list.files("flow-cytometry-data/20180919_170718-chlamee-stain", full.names = TRUE))

names(fcs_files) <- fcs_files %>% 
	gsub(pattern = ".fcs$", replacement = "")


#### Step 3: read in all the files!


read_fcs <- function(file) {
	fcs_file <- as.data.frame((exprs(read.FCS(file)))) %>% 
		clean_names()
	
}
all_fcs <- map_df(fcs_files, read_fcs, .id = "file_name")


all_fcs2 <- all_fcs %>% 
	separate(file_name, into = c("file_path", "well"), sep = -3, remove = FALSE)


colfunc <- colorRampPalette(c("white", "lightblue", "green", "yellow", "red"))


all_fcs2 %>% 
	filter(well == "B07", fl1_a > 0) %>% 
	ggplot(aes(x=log(fsc_a), y=log(fl3_a))) +
	# ylim(0, 30000) +
	# xlim(0,50000) +
	stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
	scale_fill_gradientn(colours=colfunc(100)) + 
	geom_density2d(colour="black", bins=10) +
	facet_wrap( ~ well)


all_fcs2 %>% 
	filter(grepl("B", well)) %>% 
	filter(fl1_a > 0, fl3_a > 0, fl4_a>0) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well)
ggsave("figures/chlamee-fc-stain-b-wells-fl1-3.pdf", height = 8, width =12)

all_fcs2 %>% 
	filter(grepl("C", well)) %>% 
	mutate(fl1_a = ifelse(fl1_a == 0, 1, fl1_a)) %>% 
	mutate(fl3_a = ifelse(fl3_a == 0, 1, fl3_a)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well)
ggsave("figures/chlamee-fc-stain-c-wells-fl1-3.png", height = 8, width =12)

all_fcs2 %>% 
	filter(grepl("D", well)) %>% 
	mutate(fl1_a = ifelse(fl1_a == 0, 1, fl1_a)) %>% 
	mutate(fl3_a = ifelse(fl3_a == 0, 1, fl3_a)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well)
ggsave("figures/chlamee-fc-stain-d-wells-fl1-3.png", height = 8, width =12)


all_fcs3 <- all_fcs2 %>% 
	select(1:8) %>% 
	filter(fl1_a > 5, fl3_a > 5)
	
all_fcs3 %>% 
	filter(grepl("A", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 180000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 180000, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 180000, color = "black")
ggsave("figures/chlamee-fc-stain-a-wells-fl1-3-gated.png", height = 8, width =12)

all_fcs3 %>% 
	filter(grepl("B", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 180000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 180000, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 180000, color = "black")
ggsave("figures/chlamee-fc-stain-b-wells-fl1-3-gated.png", height = 8, width =12)


all_fcs3 %>% 
	filter(grepl("C", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 180000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 180000, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 180000, color = "black")
ggsave("figures/chlamee-fc-stain-c-wells-fl1-3-gated.png", height = 8, width =12)


all_fcs3 %>% 
	filter(grepl("D", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 180000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 180000, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 180000, color = "black")
ggsave("figures/chlamee-fc-stain-d-wells-fl1-3-gated.png", height = 8, width =12)

all_fcs3 %>% 
	filter(grepl("E", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 180000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 180000, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 180000, color = "black")
ggsave("figures/chlamee-fc-stain-e-wells-fl1-3-gated.png", height = 8, width =12)

all_fcs3 %>% 
	filter(grepl("F", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 180000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 180000, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 180000, color = "black")
ggsave("figures/chlamee-fc-stain-f-wells-fl1-3-gated.png", height = 8, width =12)

all_fcs3 %>% 
	filter(grepl("G", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 180000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 180000, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 180000, color = "black")
ggsave("figures/chlamee-fc-stain-g-wells-fl1-3-gated.png", height = 8, width =12)

	
### now let's do cell counts

sorted <- all_fcs3 %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 180000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 180000, "algae", type)) 

counts <- sorted %>% 
	mutate(type = ifelse(is.na(type), "background", type)) %>% 
	group_by(well, type) %>% 
	tally()

library(readxl)
library(stringr)
plate_layout <- read_csv("data-general/plate-layout.csv") %>% 
	rename(Population = replicate)
treatments <- read_excel("data-general/ChlamEE_Treatments.xlsx") %>% 
	mutate(Population = as.character(Population))


all_treatments <- left_join(plate_layout, treatments, by= "Population") %>% 
	unite(col = unique_id, well, Population, Ancestor_ID, Treatment, sep = "_", remove = FALSE) %>% 
	mutate(well = str_to_upper(well))

all_counts <- left_join(counts, all_treatments)

all_counts %>% 
	filter(type != "background") %>% 
	mutate(biovolume = NA) %>% 
	mutate(biovolume = ifelse(type == "bacteria", ((0.5^3)*pi*1.333333)*n, biovolume)) %>%
	mutate(biovolume = ifelse(type == "algae", ((5^3)*pi*1.333333)*n, biovolume)) %>%
	ggplot(aes(x = unique_id, y = biovolume, fill = type)) + geom_bar(stat  = "identity") +
	facet_wrap( ~ Treatment, scales = "free") + scale_fill_viridis_d() +
	theme(axis.text.x = element_text(angle = 75, hjust = 1))
ggsave("figures/algae_bacteria_biovolume.png", width = 12, height = 15)

all_counts %>% 
	filter(type != "background") %>% 
	mutate(biovolume = NA) %>% 
	mutate(biovolume = ifelse(type == "bacteria", ((0.5^3)*pi*1.333333)*n, biovolume)) %>%
	mutate(biovolume = ifelse(type == "algae", ((5^3)*pi*1.333333)*n, biovolume)) %>%
	ggplot(aes(x = unique_id, y = n, fill = type)) + geom_bar(stat  = "identity") +
	facet_wrap( ~ Treatment, scales = "free") + scale_fill_viridis_d() +
	theme(axis.text.x = element_text(angle = 75, hjust = 1))
ggsave("figures/algae_bacteria_cell_counts.png", width = 12, height = 15)


ggplot(fcm, aes(x=fsc_a, y=fl1_a)) +
	ylim(0, 30000) +
	# xlim(0,50000) +
	stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
	scale_fill_gradientn(colours=colfunc(400)) + # gives the colour plot
	geom_density2d(colour="black", bins=5) # draws the lines inside



# Non-stained chlamEE -----------------------------------------------------

### let's read in all of the fcs files
fcs_filesn <- c(list.files("flow-cytometry-data/20180919_181701-chlamee-no-stain", full.names = TRUE))

names(fcs_filesn) <- fcs_filesn %>% 
	gsub(pattern = ".fcs$", replacement = "")


#### Step 3: read in all the files!


read_fcs <- function(file) {
	fcs_file <- as.data.frame((exprs(read.FCS(file)))) %>% 
		clean_names()
	
}
all_fcsn <- map_df(fcs_filesn, read_fcs, .id = "file_name")


all_fcs2n <- all_fcsn %>% 
	separate(file_name, into = c("file_path", "well"), sep = -3, remove = FALSE)


colfunc <- colorRampPalette(c("white", "lightblue", "green", "yellow", "red"))


all_fcs2n %>% 
	filter(well == "B07", fl1_a > 0) %>% 
	ggplot(aes(x=log(fsc_a), y=log(fl3_a))) +
	# ylim(0, 30000) +
	# xlim(0,50000) +
	stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
	scale_fill_gradientn(colours=colfunc(100)) + 
	geom_density2d(colour="black", bins=10) 


all_fcs2n %>% 
	filter(grepl("B", well)) %>% 
	filter(fl1_a > 0, fl3_a > 0, fl4_a>0) %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl1_a < 200000 & fl1_a > 1000 & fsc_a < 2900000 & fsc_a > 180000,  "algae", type)) %>% 
	ggplot(aes(x=fsc_a, y=fl1_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) +
	geom_hline(yintercept = 200000, color = "black") +
	geom_hline(yintercept = 1000, color = "black") +
	geom_vline(xintercept = 180000, color = "black") +
	geom_vline(xintercept = 2900000, color = "black") 
ggsave("figures/chlamee-fc-no-stain-b-wells-fsc-fl1.png", height = 8, width =12)

all_fcs2n %>% 
	filter(grepl("A", well)) %>% 
	filter(fl1_a > 0, fl3_a > 0, fl4_a>0) %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl1_a < 200000 & fl1_a > 1000 & fsc_a < 2900000 & fsc_a > 180000,  "algae", type)) %>% 
	ggplot(aes(x=fsc_a, y=fl1_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) +
	geom_hline(yintercept = 200000, color = "black") +
	geom_hline(yintercept = 1000, color = "black") +
	geom_vline(xintercept = 180000, color = "black") +
	geom_vline(xintercept = 2900000, color = "black") 
ggsave("figures/chlamee-fc-no-stain-a-wells-fsc-fl1.png", height = 8, width =12)

all_fcs2n %>% 
	filter(grepl("C", well)) %>% 
	filter(fl1_a > 0, fl3_a > 0, fl4_a>0) %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl1_a < 200000 & fl1_a > 1000 & fsc_a < 2900000 & fsc_a > 180000,  "algae", type)) %>% 
	ggplot(aes(x=fsc_a, y=fl1_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) +
	geom_hline(yintercept = 200000, color = "black") +
	geom_hline(yintercept = 1000, color = "black") +
	geom_vline(xintercept = 180000, color = "black") +
	geom_vline(xintercept = 2900000, color = "black") 
ggsave("figures/chlamee-fc-no-stain-c-wells-fsc-fl1.png", height = 8, width =12)


all_fcs2n %>% 
	filter(grepl("D", well)) %>% 
	filter(fl1_a > 0, fl3_a > 0, fl4_a>0) %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl1_a < 200000 & fl1_a > 1000 & fsc_a < 2900000 & fsc_a > 180000,  "algae", type)) %>% 
	ggplot(aes(x=fsc_a, y=fl1_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) +
	geom_hline(yintercept = 200000, color = "black") +
	geom_hline(yintercept = 1000, color = "black") +
	geom_vline(xintercept = 180000, color = "black") +
	geom_vline(xintercept = 2900000, color = "black") 
ggsave("figures/chlamee-fc-no-stain-d-wells-fsc-fl1.png", height = 8, width =12)

all_fcs2n %>% 
	filter(grepl("E", well)) %>% 
	filter(fl1_a > 0, fl3_a > 0, fl4_a>0) %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl1_a < 200000 & fl1_a > 1000 & fsc_a < 2900000 & fsc_a > 180000,  "algae", type)) %>% 
	ggplot(aes(x=fsc_a, y=fl1_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) +
	geom_hline(yintercept = 200000, color = "black") +
	geom_hline(yintercept = 1000, color = "black") +
	geom_vline(xintercept = 180000, color = "black") +
	geom_vline(xintercept = 2900000, color = "black") 
ggsave("figures/chlamee-fc-no-stain-e-wells-fsc-fl1.png", height = 8, width =12)

all_fcs2n %>% 
	filter(grepl("F", well)) %>% 
	filter(fl1_a > 0, fl3_a > 0, fl4_a>0) %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl1_a < 200000 & fl1_a > 1000 & fsc_a < 2900000 & fsc_a > 180000,  "algae", type)) %>% 
	ggplot(aes(x=fsc_a, y=fl1_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) +
	geom_hline(yintercept = 200000, color = "black") +
	geom_hline(yintercept = 1000, color = "black") +
	geom_vline(xintercept = 180000, color = "black") +
	geom_vline(xintercept = 2900000, color = "black") 
ggsave("figures/chlamee-fc-no-stain-f-wells-fsc-fl1.png", height = 8, width =12)

all_fcs2n %>% 
	filter(grepl("G", well)) %>% 
	filter(fl1_a > 0, fl3_a > 0, fl4_a>0) %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl1_a < 200000 & fl1_a > 1000 & fsc_a < 2900000 & fsc_a > 180000,  "algae", type)) %>% 
	ggplot(aes(x=fsc_a, y=fl1_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) +
	geom_hline(yintercept = 200000, color = "black") +
	geom_hline(yintercept = 1000, color = "black") +
	geom_vline(xintercept = 180000, color = "black") +
	geom_vline(xintercept = 2900000, color = "black") 
ggsave("figures/chlamee-fc-no-stain-g-wells-fsc-fl1.png", height = 8, width =12)


all_fcs3 <- all_fcs2 %>% 
	select(1:8) %>% 
	filter(fl1_a > 5, fl3_a > 5)

all_fcs3 %>% 
	filter(grepl("A", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 180000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 180000, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 180000, color = "black")
ggsave("figures/chlamee-fc-stain-a-wells-fl1-3-gated.png", height = 8, width =12)

all_fcs3 %>% 
	filter(grepl("B", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 180000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 180000, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 180000, color = "black")
ggsave("figures/chlamee-fc-stain-b-wells-fl1-3-gated.png", height = 8, width =12)


all_fcs3 %>% 
	filter(grepl("C", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 180000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 180000, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 180000, color = "black")
ggsave("figures/chlamee-fc-stain-c-wells-fl1-3-gated.png", height = 8, width =12)


all_fcs3 %>% 
	filter(grepl("D", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 180000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 180000, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 180000, color = "black")
ggsave("figures/chlamee-fc-stain-d-wells-fl1-3-gated.png", height = 8, width =12)

all_fcs3 %>% 
	filter(grepl("E", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 180000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 180000, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 180000, color = "black")
ggsave("figures/chlamee-fc-stain-e-wells-fl1-3-gated.png", height = 8, width =12)

all_fcs3 %>% 
	filter(grepl("F", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 180000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 180000, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 180000, color = "black")
ggsave("figures/chlamee-fc-stain-f-wells-fl1-3-gated.png", height = 8, width =12)

all_fcs3 %>% 
	filter(grepl("G", well)) %>% 
	# filter(well == "C01") %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 180000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 180000, "algae", type)) %>% 
	ggplot(aes(x=fl1_a, y=fl3_a, color = type)) + geom_point(alpha = 0.1, size = 0.5) + scale_y_log10() + scale_x_log10() +
	facet_wrap( ~ well) + geom_vline(xintercept = 17000, color = "black") + 
	geom_hline(yintercept = 180000, color = "black")
ggsave("figures/chlamee-fc-stain-g-wells-fl1-3-gated.png", height = 8, width =12)


### now let's do cell counts

sorted <- all_fcs3 %>% 
	mutate(type = NA) %>% 
	mutate(type = ifelse(fl3_a < 180000 & fl1_a > 17000, "bacteria", type)) %>% 
	mutate(type = ifelse(fl3_a > 180000, "algae", type)) 

counts <- sorted %>% 
	mutate(type = ifelse(is.na(type), "background", type)) %>% 
	group_by(well, type) %>% 
	tally()

library(readxl)
library(stringr)
plate_layout <- read_csv("data-general/plate-layout.csv") %>% 
	rename(Population = replicate)
treatments <- read_excel("data-general/ChlamEE_Treatments.xlsx") %>% 
	mutate(Population = as.character(Population))


all_treatments <- left_join(plate_layout, treatments, by= "Population") %>% 
	unite(col = unique_id, well, Population, Ancestor_ID, Treatment, sep = "_", remove = FALSE) %>% 
	mutate(well = str_to_upper(well))

all_counts <- left_join(counts, all_treatments)

all_counts %>% 
	filter(type != "background") %>% 
	mutate(biovolume = NA) %>% 
	mutate(biovolume = ifelse(type == "bacteria", ((0.5^3)*pi*1.333333)*n, biovolume)) %>%
	mutate(biovolume = ifelse(type == "algae", ((5^3)*pi*1.333333)*n, biovolume)) %>%
	ggplot(aes(x = unique_id, y = biovolume, fill = type)) + geom_bar(stat  = "identity") +
	facet_wrap( ~ Treatment, scales = "free") + scale_fill_viridis_d() +
	theme(axis.text.x = element_text(angle = 75, hjust = 1))
ggsave("figures/algae_bacteria_biovolume.png", width = 12, height = 15)

all_counts %>% 
	filter(type != "background") %>% 
	mutate(biovolume = NA) %>% 
	mutate(biovolume = ifelse(type == "bacteria", ((0.5^3)*pi*1.333333)*n, biovolume)) %>%
	mutate(biovolume = ifelse(type == "algae", ((5^3)*pi*1.333333)*n, biovolume)) %>%
	ggplot(aes(x = unique_id, y = n, fill = type)) + geom_bar(stat  = "identity") +
	facet_wrap( ~ Treatment, scales = "free") + scale_fill_viridis_d() +
	theme(axis.text.x = element_text(angle = 75, hjust = 1))
ggsave("figures/algae_bacteria_cell_counts.png", width = 12, height = 15)


ggplot(fcm, aes(x=fsc_a, y=fl1_a)) +
	ylim(0, 30000) +
	# xlim(0,50000) +
	stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
	scale_fill_gradientn(colours=colfunc(400)) + # gives the colour plot
	geom_density2d(colour="black", bins=5) # draws the lines inside



