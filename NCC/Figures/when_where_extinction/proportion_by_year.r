library(dplyr)
threshold<-seq(0.1, 0.5, by=0.1)
th=0.3
result<-NULL

sp_mean<-read.csv(sprintf("../../Figures/when_where_extinction_all/when_extinct.csv"))
sp_mean[which(sp_mean$exposure==" no exposure"), ]$exposure<-" no climate resilience"
sp_mean[which(sp_mean$exposure=="5-year exposure"), ]$exposure<-"climate resilience"

sp_mean_all<-sp_mean%>%dplyr::group_by(SSP, extinct_year, dispersal, label, exposure, da)%>%
  dplyr::summarise(extinct_ratio=sum(mean_extinct_ratio*all_sp)/sum(all_sp),
                   extinct_ratio2=mean(mean_extinct_ratio))
sp_mean_all$Group<-"ALL"
if (F){
  ggplot(sp_mean_all)+geom_line(aes(x=extinct_year, y=extinct_ratio2, color=factor(SSP), linetype=factor(da)))+
    facet_wrap(~exposure)
}
for (th in threshold){
  sp_mean_item<-sp_mean%>%dplyr::filter(mean_extinct_ratio>th)
  sp_mean_item<-sp_mean_item %>% 
    group_by(Group, SSP, exposure, da)%>% 
    slice(which.min(extinct_year))%>%
    dplyr::select(Group, SSP, exposure, da, extinct_year, mean_extinct_ratio)
  colnames(sp_mean_item)[6]<-"extinct_ratio"
  sp_mean_all_item<-sp_mean_all%>%dplyr::filter(extinct_ratio>th)
  sp_mean_all_item<-sp_mean_all_item %>% 
    group_by(Group, SSP, exposure, da)%>% 
    slice(which.min(extinct_year))%>%
    dplyr::select(Group, SSP, exposure, da, extinct_year, extinct_ratio)
  
  item<-bind_rows(sp_mean_item, sp_mean_all_item)
  item$extinct_threshold<-th
  if (is.null(result)){
    result<-item
  }else{
    result<-bind_rows(result, item)
  }
}

write.table(result, "../../Figures/when_where_extinction_all/proportion_by_year.csv", row.names=F, sep=",")

