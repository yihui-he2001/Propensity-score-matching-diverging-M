library(gt)
produce_table=function(q){
suppressMessages(library("dplyr"))
suppressMessages(library("xtable"))
runfile = paste0("working_data/out",q,".feather")
outfile1 = paste0("tables/ATE:T",q,".txt")
outfile2 = paste0("tables/ATT:T",q,".txt")
runs = feather::read_feather(runfile) %>% as.tibble %>%
  mutate(error = (est1 - ate), se=se1) %>%
  mutate(covered95 = abs(error)<qnorm(0.975)*se1/sqrt(N)) %>%
  mutate(covered90 = abs(error)<qnorm(0.95)*se1/sqrt(N)) 
results1 = runs %>%
  group_by(N, M) %>%
  summarise(RMSE = sqrt(mean(error^2)), MAE = mean(abs(error)),
            Coverage95 = mean(covered95), Coverage90 = mean(covered90),
            NSD=sqrt(mean((est1)^2)-mean(est1)^2)) %>%
  ungroup() %>%
  mutate_at(vars(matches("NSD")), ~ .*sqrt(N)) %>%
  mutate_at("N", ~ round(.)) %>% mutate_at("M", ~ round(.)) %>%
  arrange(N, M) %>%
  as.tibble
results1$M=as.integer(results1$M)
results1$N=as.integer(results1$N)
print(xtable(results1, digits = 3, type = "latex"), file = outfile1, include.rownames=FALSE)

runs = feather::read_feather(runfile) %>% as.tibble %>%
  mutate(error = (est2 - att), se=se2) %>%
  mutate(covered95 = abs(error)<qnorm(0.975)*se2/sqrt(N)) %>%
  mutate(covered90 = abs(error)<qnorm(0.95)*se2/sqrt(N)) 
results2 = runs %>%
  group_by(N, M) %>%
  summarise(RMSE = sqrt(mean(error^2)), MAE = mean(abs(error)),
            Coverage95 = mean(covered95), Coverage90 = mean(covered90),
            NSD=sqrt(mean((est2)^2)-mean(est2)^2)) %>%
  ungroup() %>%
  mutate_at(vars(matches("NSD")), ~ .*sqrt(N)) %>%
  mutate_at("N", ~ round(.)) %>% mutate_at("M", ~ round(.)) %>%
  arrange(N, M) %>%
  as.tibble
results2$M=as.integer(results2$M)
results2$N=as.integer(results2$N)
print(xtable(results2, digits = 3, type = "latex"), file = outfile2, include.rownames=FALSE)
}