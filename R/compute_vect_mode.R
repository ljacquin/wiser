compute_vect_mode <- function(vect) {
  freq <- table(vect)
  mode <- names(freq)[freq == max(freq)][1]
  return(mode)
}