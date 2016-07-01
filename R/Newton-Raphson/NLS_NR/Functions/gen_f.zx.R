#Generate f(z|x) from zi's and x's
gen_f.zx <- function(alpha_hat,sim.data){
  # Extract xi and zi
  xi <- extract_xi(sim.data)
  zi <- sim.data$zi
  # Allocate empty output matrix
  f.zx <- rep(NA,1000)
  
  # Vectorized: f(z|x) = expit(alpha_hat'x)     if z==1, 
  #                    = 1 - expit(alpha_hat'x) otherwise
  f.zx[which(zi==1)] <- expit(xi[which(zi==1),] %*% alpha_hat)
  f.zx[which(zi==0)] <- 1 - expit(xi[which(zi==0),] %*% alpha_hat)
  f.zx
}