gen_E.Wx <- function(alpha_hat,sim.data){
  # Extract xi
  xi <- extract_xi(sim.data)
  
  # Generate E.Wx
  tanh(xi %*% alpha_hat)
}