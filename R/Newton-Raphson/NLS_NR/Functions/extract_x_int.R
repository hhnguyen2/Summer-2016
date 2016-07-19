# EXTRACT X intercept
# Input: simulation or card data with form (1,x.1,x.2,...)
# Output: n-by-p matrix: 1,x.1,x.2,...,x.n
extract_x_int <- function(data){
  x_indices <- grep("x.",colnames(data))
  abcd <- as.matrix(data[,c(1,x_indices)])
}