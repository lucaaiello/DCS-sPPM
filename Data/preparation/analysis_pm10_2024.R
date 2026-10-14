######## ---------- PREPROCESSING OF PM10 DATA 2024 ----------- ######## 

## ---- 0. Packages ---------------------------------------------
library(data.table)
library(tidyverse)
library(ggplot2)
library(terra)
library(ncdf4)
library(sf)
library(rnaturalearth)
library(arrow)

## ---- 1. Aggregate hourly data -----------------------------
# Stream and process directly from parquet files on disk
pm10_daily_2024 <- open_dataset("data/pm10_hourly", format = "parquet") %>%
  filter(
    Start >= as.POSIXct("2024-01-01"),
    Start < as.POSIXct("2025-01-01"),
    Value > 0,
    Validity == 1
  ) %>%
  mutate(Date = as.Date(Start)) %>%
  group_by(Samplingpoint, Date) %>%
  summarise(
    Value = mean(Value, na.rm = TRUE),
    valid_hours = n(),
    .groups = "drop"
  ) %>%
  collect() 

# Save daily summary

### Uncomment ###
# write_csv(pm10_daily_2024, "pm10_europe_hourly_aggr_2024.csv")


# Load metadata CSV containing coordinates 
meta <- read.csv("metadata.csv") 
meta = meta %>% 
  select(Sampling.Point.Id, Longitude, Latitude) %>% 
  distinct(Sampling.Point.Id, .keep_all = TRUE)


## ---- 2. Load the data (hourly aggregated and daily) -------------
pm10 <- fread("pm10_europe_hourly_aggr_2024.csv")
pm10_daily <- fread("pm10_europe_daily_2024.csv")


## ---- 3. Combine hourly aggregated and daily data with same metadata ------------
pm10 = pm10 %>% 
  mutate(
    # Remove country code prefix (e.g., "AT/", "DE/", "FR/")
    SamplingPoint_clean = str_remove(Samplingpoint, "^[A-Z]{2}/")
  )

pm10_daily = pm10_daily %>% 
  filter(
    Start >= as.POSIXct("2024-01-01"),
    Start < as.POSIXct("2025-01-01"),
    Value > 0,
    Validity == 1, 
    Value < 500
  ) %>% 
  mutate(
    # Remove country code prefix (e.g., "AT/", "DE/", "FR/")
    SamplingPoint_clean = str_remove(Samplingpoint, "^[A-Z]{2}/")
  ) %>% 
  mutate(Date = as.IDate(Start)) %>% 
  select(Samplingpoint, Date, Value, SamplingPoint_clean)

pm10_2024_spatial <- pm10 %>% 
  left_join(meta, by = c("SamplingPoint_clean" = "Sampling.Point.Id"))

pm10_daily_2024_spatial <- pm10_daily %>% 
  left_join(meta, by = c("SamplingPoint_clean" = "Sampling.Point.Id"))

shared_stations <- intersect(
  unique(pm10_2024_spatial$SamplingPoint_clean),
  unique(pm10_daily_2024_spatial$SamplingPoint_clean)
)


pm10_comparison <- pm10_2024_spatial %>%
  inner_join(
    pm10_daily_2024_spatial,
    by = c("SamplingPoint_clean", "Date"), 
    suffix = c("_ds1", "_ds2")
  ) %>%
  mutate(
    # Check exact equality or near-equality (accounting for small floating-point differences)
    is_equal = Value_ds1 == Value_ds2, 
    abs_diff = abs(Value_ds1 - Value_ds2)
  )

cat("Total matching station-date records:", nrow(pm10_comparison), "\n")
cat("Exactly identical records:", sum(pm10_comparison$is_equal, na.rm = TRUE), "\n")
cat("Records with differences:", sum(!pm10_comparison$is_equal, na.rm = TRUE), "\n")

mean(pm10_comparison$abs_diff)
max(pm10_comparison$abs_diff)
quantile(pm10_comparison$abs_diff,0.9)
ggplot(pm10_comparison, aes(x=abs_diff)) + geom_histogram(bins=50)

# Combine both datasets and keep the minimum PM10 value per station and date
pm10_combined_lowest <- bind_rows(pm10_2024_spatial, pm10_daily_2024_spatial) %>%
  group_by(SamplingPoint_clean, Date) %>% 
  slice_min(order_by = Value, n = 1, with_ties = FALSE) %>% 
  ungroup()

num_stations <- n_distinct(pm10_combined_lowest$Samplingpoint)


## ---- 4. Save final dataset in .rds format -----------------------------

### Uncomment ###
# saveRDS(pm10_combined_lowest, file = "pm10_combined_2024.rds")

###### Plots ######

# # Select the first station from your dataset
# target_station <- pm10_combined_lowest$Samplingpoint[1]
# 
# pm10_combined_lowest %>% 
#   filter(Samplingpoint == target_station) %>% 
#   ggplot(aes(x = Date, y = Value)) +
#   geom_line(color = "#2b5c8f", linewidth = 0.8) +
#   geom_point(color = "#1d3557", alpha = 0.5, size = 1) +
#   labs(
#     title = paste("Daily PM10 Concentration (2024)"),
#     subtitle = paste("Station:", target_station),
#     x = "Date",
#     y = expression(paste("PM"[10], " Concentration (", mu, "g/m"^3, ")"))
#   ) +
#   theme_minimal()
# 
# target_date <- as.IDate("2024-01-01")
# 
# # Filter for target date and remove missing coordinates/PM10 values
# single_day_spatial <- pm10_combined_lowest %>% 
#   filter(
#     Date == target_date, 
#     !is.na(Latitude), 
#     !is.na(Longitude),
#     !is.na(Value)
#   ) %>% 
#   st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
# 
# # Get European map background
# europe_map <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")
# 
# ggplot() +
#   geom_sf(data = europe_map, fill = "#f8f9fa", color = "#d3d3d3") +
#   geom_sf(data = single_day_spatial, aes(color = Value), size = 1, alpha = 0.8) +
#   coord_sf(xlim = c(-10, 35), ylim = c(35, 70), expand = FALSE) +
#   scale_color_viridis_c(
#     option = "magma", 
#     name = expression(paste("PM"[10], " (", mu, "g/m"^3, ")"))
#   ) +
#   labs(
#     title = paste("PM10 Concentrations Across Europe"),
#     subtitle = paste("Date:", target_date),
#     x = "Longitude",
#     y = "Latitude"
#   ) +
#   theme_minimal()


