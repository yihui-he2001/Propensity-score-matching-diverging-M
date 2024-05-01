# Design 5 in Abadie&Imbens 2016
library(tidyverse)
create_dataset3=function(N){
  
  truncated_rexp <- function(n, rate, lower, upper) {
    x <- rexp(n, rate)
    x <- x[x >= lower & x <= upper]
    return(x)
  }
  
  X_1o=runif(N,0,1)
  X_2o=runif(N,0,1)
  X_1=map_dbl(X_1o,function(x) if (x<0.7){runif(1,-0.5,0)}else{runif(1,0,0.5)})
  X_2=map_dbl(X_2o,function(x) if (x<0.6){-truncated_rexp(35,1,0,0.5)[1]}else{truncated_rexp(35,1,0,0.5)[1]})
  U_0=rnorm(N,0,1)
  U_1=rnorm(N,0,1)
  Y_0=3*X_1-3*X_2+U_0
  Y_1=5+5*X_1+X_2+U_1
  P=exp(2*X_2)/(1+exp(2*X_2))
  W=map_dbl(P,function(data) rbinom(1,1,data))
  Y=pmap_dbl(tibble(Y_0,Y_1,W),function(Y_0,Y_1,W) Y_0*(1-W)+Y_1*W)
  ATE=4.616
  ATT=4.900
  #bate and batt are the root squares of the asymptotic variances of ATT and ATE in this dataset
  bate=2.426
  batt=2.760
  save(ATE,ATT,Y,W,X_1,X_2,P,file="working_data/dataset3.Rdata")
}

