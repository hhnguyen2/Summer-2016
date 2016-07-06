# Define function that evaluates 
# squared sum from i=1,...,n of (1,Zi)*(Ri-B0-BaZi)(1/(f.hat.zx))
opt_hrrr <- function(B,data){
  # b == Beta, where b[1] = Beta_0, and b[2] = Beta_a
  B <- as.numeric(B)
  R <- data$R
  zi <- data$zi
  f.hat.zx <- data$f.hat.zx
  
  # Evaluate (1,Zi)*(Ri-B0-BaZi)(1/(f.hat.zx))
  ans <- t(cbind(1,zi)) %*% {{R - B[1] - B[2]*zi} * {1/f.hat.zx}}
  t(ans) %*% ans
}