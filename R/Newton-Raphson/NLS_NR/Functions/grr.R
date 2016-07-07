# Define function grr that evalutes f'(alph)
grr <- function(alph,data,usem1 = FALSE){
  # INPUTS: alph, with form (alph_0, alph_1, ..., alph_n)
  # OUTPUTS: Inverse Jacobian of f(alph) 
  # ----------------
  # Extract xi
  xi <- extract_xi(data)
  # Allocate empty answer matrix (accumulation matrix)
  ans <- matrix(0,length(alph),length(alph))

  ## Evaluate inverse Jacobian. Return NA matrix if singular. 
  if(usem1){ # If expitm1, we want gradient of tanh(x) = -4 * dexpit(2*x)
    # Iteratively generate Jacobian matrix and add to allocated matrix
    for (i in 1:nrow(data)){
      # Analytic coeff is outer product of xi*t(xi), so take adv of this. 
      add.me <- -4 * dexpit(2 * sum(xi[i,] * alph)) * {xi[i,] %*% t(xi[i,])}
      ans <- ans + add.me
    }
  } else{    # Otherwise, we want gradient of expit(x) = -1 * dexpit(x)
    # Iteratively generate Jacobian matrix and add to allocated matrix
    for (i in 1:nrow(data)){
      add.me <- -1 * dexpit(sum(xi[i,] * alph)) * {xi[i,] %*% t(xi[i,])}
      ans <- ans + add.me
    }
  }
  
  # Inverse accumulation matrix 'ans' and output answer
  solve(ans)
}