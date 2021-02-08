library(dplyr)

ttt=0
threshold<-seq(0.2, 0.5, by=0.1)
th=0.2
result<-NULL
for (ttt in c(0, 1, 2)){
  sp_mean<-read.csv(sprintf("../../Figures_Full_species/when_where_extinction_all/when_extinct_%d.csv", ttt))
  sp_mean_all<-sp_mean%>%dplyr::group_by(SSP, extinct_year, dispersal, all_sp, label, exposure, da)%>%
    dplyr::summarise(extinct_ratio=mean(mean_extinct_ratio))
  sp_mean_all$Group<-"ALL"
  for (th in threshold){
    sp_mean_item<-sp_mean%>%dplyr::filter(mean_extinct_ratio>th)
    sp_mean_item<-sp_mean_item %>% 
      group_by(Group, SSP, exposure, da)%>% 
      slice(which.min(extinct_year))%>%
      dplyr::select(Group, SSP, exposure, da, extinct_year)
    
    sp_mean_all_item<-sp_mean_all%>%dplyr::filter(extinct_ratio>th)
    sp_mean_all_item<-sp_mean_all_item %>% 
      group_by(Group, SSP, exposure, da)%>% 
      slice(which.min(extinct_year))%>%
      dplyr::select(Group, SSP, exposure, da, extinct_year)
    
    item<-bind_rows(sp_mean_item, sp_mean_all_item)
    item$cell_threshold<-ttt
    item$extinct_threshold<-th
    if (is.null(result)){
      result<-item
    }else{
      result<-bind_rows(result, item)
    }
  }
}
write.table(result, "../../Figures_Full_species/when_where_extinction_all/proportion_by_year.csv", row.names=F, sep=",")
