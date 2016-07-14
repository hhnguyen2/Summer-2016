# Define function grr that evalutes f(alpha)
opt_grr <- function(alpha,data){
  # INPUTS: alpha, with form (eta_0, eta_1, ..., eta_n)
  # OUTPUTS: f(alpha) as column vector
  # ----------------
  # Extract xi
  alpha <- as.numeric(alpha) # coerce data type to ensure operation works
  W <- data$W
  xi <- extract_xi(data)
  
  # Evaluate (X)(z - expit(xi %*% alpha))
  col.minus.expit <- as.numeric(W - tanh(xi %*% alpha))  # zi - expit(alpha'x)
  ans <- as.matrix(colSums(xi * col.minus.expit))
  t(ans) %*% ans
}