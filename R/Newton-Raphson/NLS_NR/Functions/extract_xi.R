# EXTRACT Xi
# Input: simulation or card data with form (1,x.1,x.2,...)
# Output: n-by-p matrix: 1,x.1,x.2,...,x.n
extract_xi <- function(data){
  as.matrix(data[,c(1,
                        grep("x.",colnames(data)))])
}