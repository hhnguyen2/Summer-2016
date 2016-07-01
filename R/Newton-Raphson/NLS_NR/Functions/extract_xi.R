# EXTRACT Xi
# Input: simulation data
# Output: n-by-p matrix: 1,x1,x2,...,xn
extract_xi <- function(sim.data){
  as.matrix(sim.data[,c(1,
                        grep("x.",colnames(sim.data)))])
}