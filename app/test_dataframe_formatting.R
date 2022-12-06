library(dplyr)

df_macrophage_props <- read.csv("macrophage_props_1dpi_and_5dpi.csv")

df_macrophage_props

df_fishID_props  <- df_macrophage_props %>% group_by(short_name, time_point) %>% summarize(fish_id = fish_id[1], mean_area = mean(Area, na.rm = TRUE), time_in_h = min(time_in_h), dpi = dpi[1], macrophage_count = n())

df_fishID_props

df_macrophage_count  <- df_macrophage_props %>% group_by(short_name, time_point) %>% summarize(fish_id = fish_id[1], cancer_cells = cancer_cells[1],  time_in_h = min(time_in_h), dpi = dpi[1], macrophage_count = n(), mean_area = mean(Area, na.rm = TRUE))

df_macrophage_count