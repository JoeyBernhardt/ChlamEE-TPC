library(tidyverse)
library(cowplot)
library(janitor)
library(tidyverse)
library(readxl)
library(stringr)
library(flowCore)


read_fcs <- function(file) {
	fcs_file <- as.data.frame((exprs(read.FCS(file)))) %>% 
		clean_names()
	
}



fcs_files_all <- c(list.files("flow-cytometry-data/dilution-tests/",
							  full.names = TRUE, pattern = ".fcs$", recursive = TRUE))


names(fcs_files_all) <- fcs_files_all %>% 
	gsub(pattern = ".fcs$", replacement = "")

all_fcs_all <- map_df(fcs_files_all, read_fcs, .id = "file_name")

all_fcs_all$file_name[[1]]

all_fcs2_all <- all_fcs_all %>% 
	separate(file_name, into = c("file_path", "well"), sep = c("-dilution-test/"), remove = FALSE) %>% 
	mutate(concentration = "low_concentration")


write_csv(all_fcs2_all, "data-processed/dilution-test-low-density-particles.csv")
