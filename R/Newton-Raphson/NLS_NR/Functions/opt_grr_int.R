# Define function grr that evalutes f(alpha)
opt_grr_int <- function(alpha,data){
  # INPUTS: alpha, with form (alpha_0, alpha_1, ..., alpha_n)
  # OUTPUTS: tanh(2*alpha_0)^2
  # ----------------
  # Extract xi
  W <- data$W
  xi <- extract_xi(data)
  
  # Evaluate (X)(z - expit(xi %*% alpha))
  col.minus.expit <- as.numeric(W - tanh(alpha))  # zi - expit(alpha'x)
  ans <- as.matrix(colSums(xi * col.minus.expit))
  t(ans) %*% ans
}