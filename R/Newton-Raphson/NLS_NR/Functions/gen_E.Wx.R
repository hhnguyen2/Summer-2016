gen_E.Wx <- function(alpha_hat,sim.data){
  # Extract xi
  xi <- extract_xi(sim.data)
  
  # Generate E.Wx
  #tanh(xi %*% alpha_hat)         # with covariates
  tanh(sim.data$ones * alpha_hat) # intercept model
}