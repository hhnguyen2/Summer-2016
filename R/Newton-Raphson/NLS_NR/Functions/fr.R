# Define function fr that evalutes f(eta)
fr <- function(eta,sim.data,usem1 = FALSE){
  # INPUTS: eta, with form (eta_0, eta_1, ..., eta_n)
  # OUTPUTS: f(eta) as column vector
  # ----------------
  # Extract xi
  xi <- extract_xi(sim.data)
  eta <- as.numeric(eta) # coerce data type to ensure operation works
  
  # Evaluate f(eta); output as column vector
  if(usem1){ # If usem1 == TRUE, then use expitm1, since this is A|Z,x
    col.minus.expit <- as.numeric(sim.data$W - tanh(xi %*% eta)) #  W - expitm1(eta'x)
  } else{  # Otherwise just use expit, since this is Z|x
    col.minus.expit <- as.numeric(sim.data$zi - expit(xi %*% eta))  # zi - expit(eta'x)
  }
  
  # Output answer as column vector;
  as.matrix(colSums(xi * col.minus.expit))
}