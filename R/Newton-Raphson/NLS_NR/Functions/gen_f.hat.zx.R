#Generate f(z|x) from zi's and x's
gen_f.hat.zx <- function(gamma_hat,sim.data){
  # Extract xi and zi
  xi <- extract_xi(sim.data)
  zi <- sim.data$zi
  # Allocate empty output matrix
  f.hat.zx <- rep(NA,1000)
  
  # Vectorized: f(z|x) = expit(gamma_hat'x)     if z==1, 
  #                    = 1 - expit(gamma_hat'x) otherwise
  f.hat.zx[which(zi==1)] <- expit(xi[which(zi==1),] %*% gamma_hat)
  f.hat.zx[which(zi==0)] <- 1 - expit(xi[which(zi==0),] %*% gamma_hat)
  f.hat.zx
}