setwd("Y:/Script/diversity_in_e")
source("functions.r")
records<-NULL

for (g in c("Amphibians", "Birds", "Reptiles", "Mammals")){
  df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", g))
  item<-data.frame(group=g, n=nrow(df_list))
  records<-bind(records, item)
}
write.csv(records, "../../Objects/N_species_per_group.csv")
sum(records$n)
