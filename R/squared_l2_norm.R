# function to compute distance between observed and simulated phenotypes
squared_l2_norm <- function(y, y_sim) {
  return(sum((y - y_sim)^2))
}