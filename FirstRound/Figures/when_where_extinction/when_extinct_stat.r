library(dplyr)
library(Rmisc)
library(data.table)

#When removing 1 cell: What is average extinction by 2100 (+/- std dev) for each emission scenario and group? 
#What about for each emission scenario and across all groups?


ttt=0
result<-list()
for (ttt in c(0, 1, 2)){
  sp_mean<-read.csv(sprintf("../../Figures_Full_species/when_where_extinction_all/when_extinct_%d.csv", ttt))
  sp_mean<-sp_mean%>%dplyr::filter(mean_extinct_ratio>0)
  start_year<-sp_mean%>%dplyr::group_by(Group, SSP, exposure, da, all_sp)%>%
    dplyr::summarise(st_year=min(extinct_year))
  table(start_year$st_year)
  year_2021<-sp_mean%>%dplyr::filter(extinct_year==2021)
  year_2021_SSP<-year_2021%>%dplyr::group_by(SSP, Group)%>%
    dplyr::summarise(mean_extinction=mean(mean_extinct_ratio),
                     sd_extinction=sd(mean_extinct_ratio),
                     CI_extinction=CI(mean_extinct_ratio)[1]-CI(mean_extinct_ratio)[2],
                     year=2021,
                     threshold=ttt)
  year_2021_SSP_all_group<-year_2021%>%dplyr::group_by(SSP)%>%
    dplyr::summarise(mean_extinction=mean(mean_extinct_ratio),
                     sd_extinction=sd(mean_extinct_ratio),
                     CI_extinction=CI(mean_extinct_ratio)[1]-CI(mean_extinct_ratio)[2],
                     year=2021,
                     threshold=ttt,
                     Group="ALL")
  year_2021_SSP_all_group<-year_2021_SSP_all_group%>%
    dplyr::select(SSP, Group, mean_extinction, sd_extinction, CI_extinction, year, threshold)
  
  year_2100<-sp_mean%>%dplyr::filter(extinct_year==2100)
  year_2100_SSP<-year_2100%>%dplyr::group_by(SSP, Group)%>%
    dplyr::summarise(mean_extinction=mean(mean_extinct_ratio),
                     sd_extinction=sd(mean_extinct_ratio),
                     CI_extinction=CI(mean_extinct_ratio)[1]-CI(mean_extinct_ratio)[2],
                     year=2100,
                     threshold=ttt)
  year_2100_SSP_all_group<-year_2100%>%dplyr::group_by(SSP)%>%
    dplyr::summarise(mean_extinction=mean(mean_extinct_ratio),
                     sd_extinction=sd(mean_extinct_ratio),
                     CI_extinction=CI(mean_extinct_ratio)[1]-CI(mean_extinct_ratio)[2],
                     year=2100,
                     threshold=ttt,
                     Group="ALL")
  year_2100_SSP_all_group<-year_2100_SSP_all_group%>%
    dplyr::select(SSP, Group, mean_extinction, sd_extinction, CI_extinction, year, threshold)
  
  result[[sprintf("t_%d", ttt)]]<-rbindlist(list(year_2021_SSP, year_2021_SSP_all_group,
                                                 year_2100_SSP, year_2100_SSP_all_group))
}
result<-rbindlist(result)
write.table(result, "../../Objects_Full_species/when_extinction_stat/when_extinction_stat.csv",
            row.names = F, sep=",")
