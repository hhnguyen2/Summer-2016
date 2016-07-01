# Define function fr that evalutes f(eta)
test_fr <- function(eta){
  # INPUTS: eta, with form (eta_0, eta_1, ..., eta_n)
  # OUTPUTS: f(eta) as column vector
  # ----------------
  # Extract xi
  eta <- as.numeric(eta) # coerce data type to ensure operation works
  
  # Evaluate (X)(z - expit(xi %*% eta))
  col.minus.expit <- as.numeric(sim.data$zi - expit(xi %*% eta))  # zi - expit(eta'x)
  tgt <- as.matrix(colSums(xi * col.minus.expit))
  t(tgt) %*% tgt
}