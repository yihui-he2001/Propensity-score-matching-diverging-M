library(gt)
produce_table=function(q){
suppressMessages(library("dplyr"))
suppressMessages(library("xtable"))
runfile = paste0("working_data/out",q,".feather")
outfile = paste0("tables/table",q,".txt")
ate=5
runs = feather::read_feather(runfile) %>% as.tibble %>%
  mutate(error = (est - ate), se=se) %>%
  mutate(covered95 = abs(error)<qnorm(0.975)*se/sqrt(N)) %>%
  mutate(covered90 = abs(error)<qnorm(0.95)*se/sqrt(N)) 
results = runs %>%
  group_by(N, M) %>%
  summarise(RMSE = sqrt(mean(error^2)), MAE = mean(abs(error)),
            Coverage95 = mean(covered95), Coverage90 = mean(covered90),
            NSE=mean(se), NSD=sqrt(mean((est)^2)-mean(est)^2)) %>%
  ungroup() %>%
  mutate_at(vars(matches("NSD")), ~ .*sqrt(N)) %>%
  mutate_at("N", ~ round(.)) %>% mutate_at("M", ~ round(.)) %>%
  arrange(N, M) %>%
  as.tibble
print(xtable(results, digits = 3, type = "latex"), file = outfile, include.rownames=FALSE)
}