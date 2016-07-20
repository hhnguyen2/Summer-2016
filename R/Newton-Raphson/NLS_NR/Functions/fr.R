# Define function fr that evalutes f(my.gamma)
fr <- function(my.gamma,data,usem1 = FALSE){
  # INPUTS: my.gamma, with form (my.gamma_0, my.gamma_1, ..., my.gamma_n)
  # OUTPUTS: f(my.gamma) as column vector
  # ----------------
  # Extract xi
  xi <- extract_xi(data)
  my.gamma <- as.numeric(my.gamma) # coerce data type to ensure operation works
  
  # Evaluate f(my.gamma); output as column vector
  if(usem1){ # If usem1 == TRUE, then use expitm1, since this is A|Z,x; my.gamma is actually alpha
    #col.minus.expit <- as.numeric(data$W - tanh(xi %*% my.gamma))   # W - tanh(my.gamma'x)
    col.minus.expit <- as.numeric(data$W - tanh(my.gamma[1]))        # intercept-model
  } else{  # Otherwise just use expit, since this is Z|x
    col.minus.expit <- as.numeric(data$zi - expit(xi %*% my.gamma))  # zi - expit(my.gamma'x)
  }
  
  # Output answer as column vector;
  as.matrix(colSums(xi * col.minus.expit))
}