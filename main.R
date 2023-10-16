library(dplyr)
library(tidyverse)
library(FNN)
library("optparse")
library(feather)
main=function(q,T=2000,N){
# --------------------------------------------------
# source the functions
# --------------------------------------------------
source("create_dataset1.R")
source("create_dataset2.R")
source("produce_table.R")
  


# --------------------------------------------------
# setup the session
# --------------------------------------------------
set.seed(123)
option_list = list(
    make_option(c("--output_feather_path"), type="character", default=paste0("working_data/out",q,".feather"), 
                help="output filename", metavar="character"),
    make_option(c("--resume"), type="logical", default=TRUE, 
                help="Resume from previous intermediate results?", metavar="logical"),
    make_option(c("--N_runs"), type="integer", default=2000L, 
                help="Number of runs per estimator", metavar="integer"),
    make_option(c("--N_chunk"), type="integer", default=200L, 
                help="After N_chunk runs, intermediate results will be saved", metavar="integer"),
    make_option(c("--N_workers"), type="integer", default=10, 
                help="If n_workers > 1, runs will be parallelized using 'future.apply'", metavar="integer")
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (opt$N_workers > 1) {
    library("future.apply")
    plan(multisession, workers = opt$N_workers)
    options(future.globals.maxSize = 3000 * 1024 ^ 2)
    replicate = future_replicate
}



# --------------------------------------------------
# start iterations
# --------------------------------------------------
for (k in 1:T){
print(paste0("Iteration:",k))
# --------------------------------------------------
# generate and load data
# --------------------------------------------------
if (q==1){
create_dataset1(N)
load("working_data/dataset1.Rdata")
}else{
create_dataset2(N)
load("working_data/dataset2.Rdata")
}



# --------------------------------------------------
# estimate propensity scores
# --------------------------------------------------
data=tibble(W=W,X_1=X_1,X_2=X_2)
model=glm(W~X_1+X_2,data,family=binomial(link="logit"))
hat_P=predict(model,data,type="response")



# --------------------------------------------------
# match with estimated propensity scores
# --------------------------------------------------
M_max=2^(floor(log2(N)/2))
Q=round(N^{1/3})
L=4
vec0=which(W==0)
vec1=which(W==1)
Index01=get.knnx(hat_P[vec0], hat_P[vec1], k = M_max, algorithm = "kd_tree")$nn.index
Index10=get.knnx(hat_P[vec1], hat_P[vec0], k = M_max, algorithm = "kd_tree")$nn.index
Index00=get.knnx(hat_P[vec0], hat_P[vec0], k = Q, algorithm = "kd_tree")$nn.index
Index11=get.knnx(hat_P[vec1], hat_P[vec1], k = Q, algorithm = "kd_tree")$nn.index
for (m in c(2^{0:log2(M_max)},round(N^{1/3}))){
K=rep(0,N)
K[vec0[1:length(tabulate(c(Index01[,1:m])))]]=tabulate(c(Index01[,1:m]))
K[vec1[1:length(tabulate(c(Index10[,1:m])))]]=tabulate(c(Index10[,1:m]))
tau_local=Y
tau_local[vec0]=tau_local[vec0]-apply(as.matrix(Index10[,1:m]),1,function(x) mean(Y[vec1[x]]))
tau_local[vec1]=tau_local[vec1]-apply(as.matrix(Index01[,1:m]),1,function(x) mean(Y[vec0[x]]))
s2_local=rep(0,N)
s2_local[vec0]=Q/(Q-1)*(apply(Index00,1,function(x) mean(Y[vec0[x]]^2))-apply(Index00,1,function(x) mean(Y[vec0[x]]))^2)
s2_local[vec1]=Q/(Q-1)*(apply(Index11,1,function(x) mean(Y[vec1[x]]^2))-apply(Index11,1,function(x) mean(Y[vec1[x]]))^2)
hat_tau=mean(pmap_dbl(tibble(Y,K,W),function(Y,K,W) Y*(K/m+1)*(2*W-1)))
hat_s2=mean(pmap_dbl(tibble(K,tau_local,s2_local),function(K,tau_local,s2_local) tau_local^2+((K/m)^2+(2*m-1)*K/(m^2))*s2_local))-hat_tau^2
c0=matrix(0,nrow=N,ncol=2)
c0[vec0,]=L/(L-1)*t(apply(Index00[,1:L],1,function(x) c(mean(Y[vec0[x]]*X_1[vec0[x]]),mean(Y[vec0[x]]*X_2[vec0[x]])))-apply(Index00[,1:L],1,function(x) c(mean(Y[vec0[x]])*mean(X_1[vec0[x]]),mean(Y[vec0[x]])*mean(X_2[vec0[x]]))))
c0[vec1,]=L/(L-1)*t(apply(Index01[,1:L],1,function(x) c(mean(Y[vec0[x]]*X_1[vec0[x]]),mean(Y[vec0[x]]*X_2[vec0[x]])))-apply(Index01[,1:L],1,function(x) c(mean(Y[vec0[x]])*mean(X_1[vec0[x]]),mean(Y[vec0[x]])*mean(X_2[vec0[x]]))))
c1=matrix(0,nrow=N,ncol=2)
c1[vec0,]=L/(L-1)*t(apply(Index10[,1:L],1,function(x) c(mean(Y[vec1[x]]*X_1[vec1[x]]),mean(Y[vec1[x]]*X_2[vec1[x]])))-apply(Index10[,1:L],1,function(x) c(mean(Y[vec1[x]])*mean(X_1[vec1[x]]),mean(Y[vec1[x]])*mean(X_2[vec1[x]]))))
c1[vec1,]=L/(L-1)*t(apply(Index11[,1:L],1,function(x) c(mean(Y[vec1[x]]*X_1[vec1[x]]),mean(Y[vec1[x]]*X_2[vec1[x]])))-apply(Index11[,1:L],1,function(x) c(mean(Y[vec1[x]])*mean(X_1[vec1[x]]),mean(Y[vec1[x]])*mean(X_2[vec1[x]]))))
c_mat=apply(cbind(c0,c1,hat_P),1,function(x) x[1:2]*x[5]+x[3:4]*(1-x[5]))
hat_c=apply(c_mat,1,mean)
X=map2(X_1,X_2,function(a,b) c(a,b))
I_list=map2(hat_P,X,function(p,x) p*(1-p)*matrix(x %*% t(x), nrow = length(x), ncol = length(x)))
hat_I=Reduce(`+`, I_list) / length(I_list)
hat_var=hat_s2-t(hat_c) %*%  solve(hat_I) %*% hat_c
hat_s2=sqrt(hat_s2)
hat_var=sqrt(hat_var)



# --------------------------------------------------
# record the results in the feather file
# --------------------------------------------------
if (opt$resume & file.exists(opt$output_feather_path)) {runs = as.data.frame(feather::read_feather(opt$output_feather_path))} else {runs = NULL}
runs =rbind(runs, data.frame(N=N, M=m, est=hat_tau, se=hat_var))
feather::write_feather(runs, opt$output_feather_path)
} 
}
}


for (q in 1:1){
  file.remove(paste0("working_data/out",q,".feather"))
  for (N in 2^{9:13}){
  print(paste0("Design ",q,":N=",N))
  main(q,T=2000,N)
  }
  produce_table(q)
}
