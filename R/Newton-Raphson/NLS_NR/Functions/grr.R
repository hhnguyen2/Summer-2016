# Define function grr that evalutes f'(alph)
grr <- function(alph,data,usem1 = FALSE){
  # INPUTS: alph, with form (alph_0, alph_1, ..., alph_n)
  # OUTPUTS: Inverse Jacobian of f(alph) 
  # ----------------
  # Extract xi
  xi <- extract_xi(data)
  # Allocate empty answer matrix
  ans <- matrix(0,length(alph),length(alph))
  
  # Iteratively generate Jacobian matrix and add to allocated matrix
  for (i in 1:nrow(data)){
    # Analytic coeff is outer product of xi*t(xi), so take adv of this. 
    add.me <- -1 * dexpit(sum(xi[i,] * alph)) * {xi[i,] %*% t(xi[i,])}
    ans <- ans + add.me
  }
  
  # Evaluate inverse Jacobian. Return NA matrix if singular. 
  if(usem1){ # If using expitm1, dexpit should be scaled by 2
    solve(ans*4) 
  } else{    # Otherwise, just solve as is 
    solve(ans)
  }
  
}