# The following function estimates the relative risk RR and the confidence bounds
# on it for a population group g, an array bp of hospitalization counts in different
# population groups (including g) before the peak, the corresponding array ap after
# the peak, and the sample size N (we generally use N = 100,000) according to the 
# method described in the Epidemics paper 
# http//www.sciencedirect.com/science/article/pii/S1755436515000511

# Input:   g = group g
#         bp = vector of before peak case counts, with length g; bp[g] = before peak count for group g
#         ap = vector of after peak case counts, with length g;  ap[g] = after peak count for group g
#          N = population size

est_ba = function(g,bp,ap,N){
  l = length(bp)
  zz=rep(0,N)
  for(i in 1:N){
    x=rep(0,l)
    y=rep(0,l)
    
    for(j in 1:l){
      x[j]=rgamma(1,(bp[j]+1),1)
      y[j]=rgamma(1,(ap[j]+1),1)
    }
    
    f1=x[g]/sum(x)
    f2=y[g]/sum(y)
    zz[i]=f1/f2
  }
  v = rep(0,3)
  v[1]=mean(zz)
  v[2:3]=quantile(zz,c(0.025,0.975))
  v
}