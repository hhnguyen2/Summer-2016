## Working with hospital data
hosd <- read.csv("Ed_1988_2012.csv")
source("est_ba.R")

## Array of weekly RSV hosptalization counts in the 10 age groups for the 47 states
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

## Determine states that reported for each season 
sy=matrix(0,nrow=47,ncol=11)
for(i in 1:47){
  for(j in 1:11){
    if(sum(rs[i,1:10,j,1:25])>0&sum(rs[i,1:10,j,28:52])>0){sy[i,j]=1}
  }}

## Determine peak week for each state during each season for which it reported
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

## Find RR Estimate per state
Rs=matrix(0,nrow=11,ncol=10)
for(yr in 1:11){
  for(ag in 1:10){
    x=0
    y=0
    bp=0
    ap=0
    for(st in 1:47){
      mw=pws[st,yr]
      if(mw>0){
        bp=bp+sum(rs[st,1:10,yr,1:(mw-2)])
        ap=ap+sum(rs[st,1:10,yr,(mw+2):52])
        x=x+sum(rs[st,ag,yr,1:(mw-2)])
        y=y+sum(rs[st,ag,yr,(mw+2):52])
      }
    }
    Rs[yr,ag]=(x/bp)/(y/ap)
  }}

## Use est_ba to find point estimate and confidence interval
est_ba(g,bp,ap,N)