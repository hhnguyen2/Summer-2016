# dexpit: expit'(x) = expit(x)(1-expit(x))
dexpit <- function(x){
  expit(x)*{1 - expit(x)}
}