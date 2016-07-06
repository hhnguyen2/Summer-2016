# Define function fr that evalutes f(gamma)
opt_fr <- function(gamma,data){
  # INPUTS: gamma, with form (gamma_0, gamma_1, ..., gamma_n)
  #         sim.data, in form (1,x.1,x.2,...)
  # OUTPUTS: f(gamma) as column vector
  # ----------------
  # Extract xi
  gamma <- as.numeric(gamma) # coerce data type to ensure operation works
  xi <- extract_xi(data)
  zi <- data$zi
  
  # Evaluate (X)(z - expit(xi %*% gamma))
  col.minus.expit <- as.numeric(zi - expit(xi %*% gamma))  # zi - expit(gamma'x)
  ans <- as.matrix(colSums(xi * col.minus.expit))
  
  # output answer
  t(ans) %*% ans
}