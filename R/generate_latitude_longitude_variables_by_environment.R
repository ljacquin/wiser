# function to generate latitude and longitude variables by environment
generate_latitude_longitude_variables_by_environment <- function(df) {
  for (site in unique(df$Site)) {
    for (year in unique(df$Year)) {
      for (block_num in unique(df$Block)) {
        # variables new names
        latitude_col_name <- paste0(
          substr(site, 1, 2), "_", year,
          "_block_", block_num, "_latitude"
        )
        longitude_col_name <- paste0(
          substr(site, 1, 2), "_", year,
          "_block_", block_num, "_longitude"
        )
        # add new columns
        df[[latitude_col_name]] <- ifelse(
          df$Site == site & df$Year == year &
            df$Block == block_num,
          df$Latitude, 0
        )
        df[[longitude_col_name]] <- ifelse(
          df$Site == site & df$Year == year &
            df$Block == block_num,
          df$Longitude, 0
        )
      }
    }
  }
  return(df)
}
