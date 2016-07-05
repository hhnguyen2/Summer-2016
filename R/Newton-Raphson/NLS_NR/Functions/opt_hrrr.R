# Define function that evaluates 
# squared sum from i=1,...,n of (1,Zi)*(Ri-B0-BaZi)(1/(f.hat.zx))
opt_hrrr <- function(B){
  # b == Beta, where b[1] = Beta_0, and b[2] = Beta_a
  B <- as.numeric(B)
  
  # Evaluate (1,Zi)*(Ri-B0-BaZi)(1/(f.hat.zx))
  ans <- t(cbind(1,sim.data$zi)) %*% {{sim.data$R - B[1] - B[2]*sim.data$zi} * {1/sim.data$f.hat.zx}}
  t(ans) %*% ans
}