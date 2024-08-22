find_columns_with_multiple_unique_values <- function(df_) {
  cols_with_multiple_values <- colnames(df_)[apply(
    df_, 2, function(col) length(unique(col)) > 1
  )]
  return(cols_with_multiple_values)
}