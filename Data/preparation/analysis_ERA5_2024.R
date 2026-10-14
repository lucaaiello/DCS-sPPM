######## ---------- PREPROCESSING OF ERA5 DATA 2024 ----------- ######## 

## ---- 0. Packages ---------------------------------------------
library(ggplot2)
library(ncdf4)

## ---- 1. Load each of the four 3-month NetCDF datasets (separated into 4 datasets, one per variable) ----

### Choose the corresponding files and uncomment ###

# nc_file_temp <- nc_open("ERA5_2024/01_02_03_2024_euro/2m_temperature_0_daily-mean.nc")
# nc_file_prec <- nc_open("ERA5_2024/01_02_03_2024_euro/total_precipitation_0_daily-mean.nc")
# nc_file_wind_u <- nc_open("ERA5_2024/01_02_03_2024_euro/10m_u_component_of_wind_stream-oper_daily-mean.nc")
# nc_file_wind_v <- nc_open("ERA5_2024/01_02_03_2024_euro/10m_v_component_of_wind_0_daily-mean.nc")

# nc_file_temp <- nc_open("ERA5_2024/04_05_06_2024_euro/2m_temperature_0_daily-mean.nc")
# nc_file_prec <- nc_open("ERA5_2024/04_05_06_2024_euro/total_precipitation_0_daily-mean.nc")
# nc_file_wind_u <- nc_open("ERA5_2024/04_05_06_2024_euro/10m_u_component_of_wind_stream-oper_daily-mean.nc")
# nc_file_wind_v <- nc_open("ERA5_2024/04_05_06_2024_euro/10m_v_component_of_wind_0_daily-mean.nc")

# nc_file_temp <- nc_open("ERA5_2024/07_08_09_2024_euro/2m_temperature_0_daily-mean.nc")
# nc_file_prec <- nc_open("ERA5_2024/07_08_09_2024_euro/total_precipitation_0_daily-mean.nc")
# nc_file_wind_u <- nc_open("ERA5_2024/07_08_09_2024_euro/10m_u_component_of_wind_stream-oper_daily-mean.nc")
# nc_file_wind_v <- nc_open("ERA5_2024/07_08_09_2024_euro/10m_v_component_of_wind_0_daily-mean.nc")

nc_file_temp <- nc_open("ERA5_2024/10_11_12_2024_euro/2m_temperature_0_daily-mean.nc")
nc_file_prec <- nc_open("ERA5_2024/10_11_12_2024_euro/total_precipitation_0_daily-mean.nc")
nc_file_wind_u <- nc_open("ERA5_2024/10_11_12_2024_euro/10m_u_component_of_wind_stream-oper_daily-mean.nc")
nc_file_wind_v <- nc_open("ERA5_2024/10_11_12_2024_euro/10m_v_component_of_wind_0_daily-mean.nc")


# Print file metadata and variable names
# names(nc_file_temp$var)
# names(nc_file_prec$var)
# names(nc_file_wind_u$var)
# names(nc_file_wind_v$var)


## ---- 2. Extract and create dataset from 4 variables -----------------------------

latitude <- ncvar_get(nc_file_temp, "latitude")
longitude <- ncvar_get(nc_file_temp, "longitude")
time_raw <- ncvar_get(nc_file_temp, "valid_time") 

temp_array <- ncvar_get(nc_file_temp, "t2m") - 273.15 # to celsius
prec_array <- ncvar_get(nc_file_prec, "tp")
wind_u_array <- ncvar_get(nc_file_wind_u, "u10")
wind_v_array <- ncvar_get(nc_file_wind_v, "v10")
wind_array <- sqrt(wind_u_array^2 + wind_v_array^2) # wind intensity

# Always close the file connection when done
nc_close(nc_file_temp)
nc_close(nc_file_prec)
nc_close(nc_file_wind_u)
nc_close(nc_file_wind_v)

# Spatio-temporal grid for new dataset
df_combined <- expand.grid(
  Longitude = longitude,
  Latitude  = latitude,
  Date      = time_raw
)

df_combined$Temperature   <- as.vector(temp_array)
df_combined$Precipitation <- as.vector(prec_array)
df_combined$Wind           <- as.vector(wind_array)

## ---- 3. Test the temporal and spatial plots -----------------------------

### Temporal plot ###

# Extraction 
temp_station <- temp_array[1,1, ] 
prec_station <- prec_array[1,1, ] 
wind_station <- wind_array[1,1, ] 

df_station <- data.frame(
  Date = time_raw,
  Temperature = temp_station,
  Precipitation = prec_station,
  Wind = wind_station
)

ggplot(df_station, aes(x = Date, y = Temperature)) +
  geom_line(color = "firebrick", linewidth = 1) +
  theme_minimal() +
  labs(
    x = "Date",
    y = "Température (°C)"
  )

ggplot(df_station, aes(x = Date, y = Precipitation)) +
  geom_line(color = "firebrick", linewidth = 1) +
  theme_minimal() +
  labs(
    x = "Date",
    y = "Precipitation"
  )

ggplot(df_station, aes(x = Date, y = Wind)) +
  geom_line(color = "firebrick", linewidth = 1) +
  theme_minimal() +
  labs(
    x = "Date",
    y = "Wind intensity"
  )

### Spatial plot ###

day_idx <- 1
selected_date <- time_raw[day_idx]

grid_data <- expand.grid(Longitude = longitude, Latitude = latitude)
grid_data$Temperature <- as.vector(temp_array[, , day_idx]) 
grid_data$Precipitation <- as.vector(prec_array[, , day_idx]) 
grid_data$Wind <- as.vector(wind_array[, , day_idx]) 

ggplot(grid_data, aes(x = Longitude, y = Latitude, fill = Temperature)) +
  geom_raster(interpolate = FALSE) +
  scale_fill_viridis_c(option = "magma", name = "Temp (°C)") +
  coord_quickmap() +
  theme_minimal() +
  labs(
    title = paste("Temperature -", selected_date),
    x = "Longitude",
    y = "Latitude"
  )

ggplot(grid_data, aes(x = Longitude, y = Latitude, fill = Precipitation)) +
  geom_raster(interpolate = FALSE) +
  scale_fill_viridis_c(option = "magma", name = "Temp (°C)") +
  coord_quickmap() +
  theme_minimal() +
  labs(
    title = paste("Precipitation -", selected_date),
    x = "Longitude",
    y = "Latitude"
  )

ggplot(grid_data, aes(x = Longitude, y = Latitude, fill = Wind)) +
  geom_raster(interpolate = FALSE) +
  scale_fill_viridis_c(option = "magma", name = "Temp (°C)") +
  coord_quickmap() +
  theme_minimal() +
  labs(
    title = paste("Wind -", selected_date),
    x = "Longitude",
    y = "Latitude"
  )


## ---- 3. Save each of the four 3-month datasets -----------------------------

### Choose the corresponding file and uncomment ###
# saveRDS(df_combined, file = "era5_combined_01_02_03_2024.rds")
# saveRDS(df_combined, file = "era5_combined_04_05_06_2024.rds")
# saveRDS(df_combined, file = "era5_combined_07_08_09_2024.rds")
# saveRDS(df_combined, file = "era5_combined_10_11_12_2024.rds")

## ---- 4. Combine the four 3-month datasets -----------------------------

df_combined_123 <- readRDS("era5_combined_01_02_03_2024.rds")
df_combined_456 <- readRDS("era5_combined_04_05_06_2024.rds")
df_combined_789 <- readRDS("era5_combined_07_08_09_2024.rds")
df_combined_101112 <- readRDS("era5_combined_10_11_12_2024.rds")

df_final <- bind_rows(df_combined_123,df_combined_456,df_combined_789,df_combined_101112)

df_final <- df_final[order(df_final$Longitude, df_final$Latitude, df_final$Date), ]

df_site <- df_final %>%
  filter(Longitude == -11, Latitude == 72) %>%
  arrange(Date)

ggplot(df_site, aes(x = Date, y = Temperature)) +
  geom_line(color = "firebrick", linewidth = 0.8) +
  labs(
    title = "Temperature Profile at (Long: -11, Lat: 72)",
    x = "Date",
    y = "Temperature"
  ) +
  theme_minimal()


## ---- 5. Save final dataset in .rds format -----------------------------

### Uncomment ###
# saveRDS(df_final, file = "era5_combined_2024.rds")
