## Read in hospital data
hosd <- read.csv("Ed_1988_2012.csv")
RR_est.RS <- read.csv("RR_RR_est.csv")

## Function to obtain RR RR_estimates and confidence bounds
source("est_ba.R")

## rs: Array of weekly RSV hosptalization counts in the 10 age groups for the 47 states
rs=array(0,dim=c(47,10,11,52)) ## 47 states : 10 age groups : 11 years : 52 weeks (1 yr)
sc=1
for(i in 2:nrow(hosd)){
  # state counter updated once new state is encountered
  if(hosd[i,1] != hosd[(i-1),1]){sc=sc+1}
  if(hosd[i,2]>=706 & hosd[i,2]<=1278 & hosd[i,2]!=1018){
    # week number assigned for study period
    if (hosd[i,2]<1018){wk=hosd[i,2]-705}
    if(hosd[i,2]>1018){wk=hosd[i,2]-706}
    ag=hosd[i,4]+1
    if(ag>=10){ag=10}
    yr=ceiling(wk/52)
    w=wk-((yr-1)*52)
    rs[sc,ag,yr,w]=rs[sc,ag,yr,w]+hosd$RSV[i]
  }
}

## sy: Indicator matrix for whether states reported for each season 
sy=matrix(0,nrow=47,ncol=11)
for(i in 1:47){
  for(j in 1:11){
    if(sum(rs[i,1:10,j,1:25])>0&sum(rs[i,1:10,j,28:52])>0){sy[i,j]=1}
  }}

## pws: Peak week for each state during each season for which it reported
pws=matrix(0,nrow=47,ncol=11)
for(st in 1:47){
  for (yr in 1:11){
    if(sy[st,yr]==1){
      sm=0
      for(wk in 1:52){
        x=sum(rs[st,1:10,yr,wk])
        if(x>sm){
          sm=x
          pws[st,yr]=wk
        }}
    }}}

## Pre-allocate before-peak (bp) and after-peak (ap) matrices
bp=ap=RR_est=RR_est.lwb=RR_est.upb=matrix(0,nrow=11,ncol=10)
rownames(bp) = rownames(ap) = rownames(RR_est) = rownames(RR_est.lwb) = rownames(RR_est.upb)
  c("2001-02","2002-03","2003-04","2004-05","2005-06","2006-07","2007-08","2008-09","2009-10","2010-11","2011-12")
colnames(bp) = colnames(ap) = colnames(RR_est) = colnames(RR_est.lwb) = colnames(RR_est.upb)
  c("<1y","y,1","y,2","y,3","y,4","y,5-11","y,12-17","y,18-49","y,50-64","y,65+")

for(yr in 1:11){
  for(ag in 1:10){
    x=0
    y=0
    for(st in 1:47){
      mw=pws[st,yr]
      if(mw>0){
        x=x+sum(rs[st,ag,yr,1:(mw-2)])
        y=y+sum(rs[st,ag,yr,(mw+2):52])
      }
    }
    bp[yr,ag]=x
    ap[yr,ag]=y
  }
}

## Use est_ba to find point RR_estimate and confidence interval
N=100000
set.seed(65711)

for(g in 1:10){
  V <- est_ba(bp[g,],ap[g,],N)
  RR_est[g,] <- V[,1]
  RR_est.lwb[g,] <- V[,2]
  RR_est.upb[g,] <- V[,3]
}

## Organize results in easy to read fashion
RR_est <- matrix(V[,1],11,10)
RR_est.lwb <- matrix(V[,2],11,10)
RR_est.upb <- matrix(V[,3],11,10)

RR_est
RR_est.lwb
RR_est.upb