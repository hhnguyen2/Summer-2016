# Define function fr that evalutes f(gamma)
# NOTE: Assumes that data frame sim.data exists. Generate using gen_sim.data()
opt_fr <- function(gamma){
  # INPUTS: gamma, with form (gamma_0, gamma_1, ..., gamma_n)
  # OUTPUTS: f(gamma) as column vector
  # ----------------
  # Extract xi
  gamma <- as.numeric(gamma) # coerce data type to ensure operation works
  
  # Evaluate (X)(z - expit(xi %*% gamma))
  col.minus.expit <- as.numeric(sim.data$zi - expit(xi %*% gamma))  # zi - expit(gamma'x)
  ans <- as.matrix(colSums(xi * col.minus.expit))
  
  # output answer
  t(ans) %*% ans
}