### Design 1 in Abadie&Imbens 2016
library(tidyverse)
create_dataset1=function(N){
X_1=runif(N,-0.5,0.5)
X_2=runif(N,-0.5,0.5)
U_0=rnorm(N,0,1)
U_1=rnorm(N,0,1)
Y_0=3*X_1-3*X_2+U_0
Y_1=5+5*X_1+X_2+U_1
P=exp(X_1+2*X_2)/(1+exp(X_1+2*X_2))
W=map_dbl(P,function(data) rbinom(1,1,data))
Y=pmap_dbl(tibble(Y_0,Y_1,W),function(Y_0,Y_1,W) Y_0*(1-W)+Y_1*W)
ATE=5
save(ATE,Y,W,X_1,X_2,P,file="working_data/dataset1.Rdata")
}