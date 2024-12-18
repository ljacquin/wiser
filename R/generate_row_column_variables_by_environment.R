# function to generate row and column variables by environment
generate_row_column_variables_by_environment <- function(df) {
  for (site in unique(df$Site)) {
    for (year in unique(df$year)) {
      for (management in unique(df$Management)) {
        for (block_num in unique(df$block)) {
          # variables new names
          row_col_name <- paste0(
            substr(site, 1, 3), year, management,
            "_block_", block_num, "_row"
          )
          pos_col_name <- paste0(
            substr(site, 1, 3), year, management,
            "_block_", block_num, "_column"
          )
          # add new columns
          df[[row_col_name]] <- ifelse(
            df$Site == site & df$year == year &
              df$Management == management & df$block == block_num,
            df$Row, 0
          )
          df[[pos_col_name]] <- ifelse(
            df$Site == site & df$year == year &
              df$Management == management & df$block == block_num,
            df$Column, 0
          )
        }
      }
    }
  }
  return(df)
}