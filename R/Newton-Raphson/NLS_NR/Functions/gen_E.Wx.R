gen_E.Wx <- function(omega_hat,sim.data){
  # Extract xi
  xi <- extract_xi(sim.data)
  
  # Generate E.Wx
  expitm1(xi %*% omega_hat)
}