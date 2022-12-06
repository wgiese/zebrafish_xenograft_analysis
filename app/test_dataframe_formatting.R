library(dplyr)

df_macrophage_props <- read.csv("/home/wgiese/image-data-workflows/zebrafish-xenograph-model/zebrafish_xenograft_analysis/app/macrophage_props_1dpi_and_5dpi.csv")

df_macrophage_props

df_fishID_props  <- df_macrophage_props %>% group_by(fish_id,time_point)

df_fishID_props
