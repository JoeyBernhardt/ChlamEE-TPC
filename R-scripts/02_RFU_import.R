## compare RFUs to flow cytometry cell counts

library(tidyverse)
library(readxl)

rfu_raw <- read_excel("data-raw/rfu-chlamee-2018-09-19.xlsx", range = "B72:N80") %>% 
	rename(row = X__1) %>% 
	gather(key = column, value = RFU, 2:13) %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") 

chlamee_inoc_densities <- read_csv("data-processed/chlamEE-inoc-densities-sep-19-2018.csv")

all_inocs <- left_join(chlamee_inoc_densities, rfu_raw, by = "well")

all_inocs %>% 
	ggplot(aes(x = cell_concentration, y = RFU, color = Treatment)) +
	geom_point(size = 3) + scale_color_viridis_d() +
	xlab("Cells per ml")
ggsave("figures/RFU-cells.png", width = 6, height = 6)
