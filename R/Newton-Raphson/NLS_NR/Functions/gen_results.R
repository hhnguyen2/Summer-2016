# gen_results
# INPUT: a vector of Ba_hat estimates
# OUTPUT: 
gen_results <- function(Ba_hat,mu){
  ## Generate statistics
  true_mu <- mu[1]
  Ba_hat.avg <- mean(Ba_hat)
  Ba_hat.bias <- Ba_hat.avg - true_mu
  Ba_hat.pcbias <- (Ba_hat.bias / true_mu) * 100
  Ba_hat.var <- var(Ba_hat)
  Ba_hat.SE <- sd(Ba_hat)
  Ba_hat.lcl <- Ba_hat.avg - 1.96*Ba_hat.SE
  Ba_hat.ucl <- Ba_hat.avg + 1.96*Ba_hat.SE
  Ba_hat.mean2 <- Ba_hat.bias^2 + Ba_hat.var
  
  ## Tabulated results
  data.frame(mean=Ba_hat.avg,
             bias=Ba_hat.bias,
             pc_bias=Ba_hat.pcbias,
             var=Ba_hat.var,
             mean2_err=Ba_hat.mean2,
             lwr_95_ci=Ba_hat.lcl,
             upr_95_ci=Ba_hat.ucl)
}