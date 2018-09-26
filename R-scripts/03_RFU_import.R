### RFU import for the Anc4 pilot


RFU_files <- c(list.files("data-raw/anc4-pilot-september2018/anc4-pilot-2018-09-26", full.names = TRUE))

RFU_files <- RFU_files[grepl(".xlsx", RFU_files)]

names(RFU_files) <- RFU_files %>% 
	gsub(pattern = ".xlsx$", replacement = "")

all_plates <- map_df(RFU_files, read_excel, range = "B56:N64", .id = "file_name") %>% 
	rename(row = X__1)


all_temp_RFU <- all_plates %>% 
	gather(key = column, value = RFU, 3:14) %>% 
	unite(row, column, col = "well", remove = FALSE, sep = "") %>% 
	mutate(column = formatC(column, width = 2, flag = 0)) %>% 
	mutate(column = str_replace(column, " ", "0")) %>% 
	unite(col = well, row, column, sep = "") %>% 
	filter(!is.na(RFU))
	
plate_layout <- read_excel("data-general/plate-layout-anc4-2018-09-26.xlsx") %>% 
	select(-population)
color_key <- read_excel("data-general/colour-key-anc4-2018-09-26.xlsx")

plate2 <- left_join(plate_layout, color_key, by = "colour") %>% 
	mutate(well = str_to_upper(well)) %>% 
	filter(!is.na(population)) %>% 
	rename(Population = population)

treatments <- read_excel("data-general/ChlamEE_Treatments.xlsx") %>% 
	mutate(Population = as.character(Population))


all_treatments <- left_join(plate2, treatments, by = "Population") %>% 
	unite(col = unique_id, well, Population, Ancestor_ID, Treatment, sep = "_", remove = FALSE) %>% 
	mutate(well = str_to_upper(well))

all_rfus <- left_join(all_temp_RFU, all_treatments) %>% 
	separate(col = file_name, into = c("path", "plate"), sep = "plate", remove = FALSE)

all_rfus %>% 
	mutate(Treatment = ifelse(is.na(Treatment), "COMBO", Treatment)) %>% 
	ggplot(aes(x = Treatment, y = RFU, color = Treatment)) + geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
	facet_wrap( ~ plate)
ggsave("figures/anc4-pilot-RFU-2018-09-26.pdf", width = 14, height =8)
