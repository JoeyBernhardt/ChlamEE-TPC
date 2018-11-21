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



fcs_files_all <- c(list.files("flow-cytometry-data/dilution-tests/20181121_162441-high-density-dilution-count-test",
							  full.names = TRUE, pattern = ".fcs$", recursive = TRUE))


names(fcs_files_all) <- fcs_files_all %>% 
	gsub(pattern = ".fcs$", replacement = "")

all_fcs_all <- map_df(fcs_files_all, read_fcs, .id = "file_name")

all_fcs_all$file_name[[1]]

all_fcs2_all <- all_fcs_all %>% 
	separate(file_name, into = c("file_path", "well"), sep = c("-count-test/"), remove = FALSE) %>% 
	mutate(concentration = "high_concentration")


write_csv(all_fcs2_all, "data-processed/dilution-test-high-density-particles.csv")
