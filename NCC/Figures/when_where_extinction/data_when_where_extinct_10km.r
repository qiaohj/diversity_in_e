library(dplyr)
library(data.table)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
exposure<-5
dispersal<-1
if (F){
  for (dispersal in c(0, 1)){
    for (exposure in c(0, 5)){
      print(paste(dispersal, exposure))
      sp_dis_all<-readRDS(sprintf("../../Figures/N_Extinction/sp_dis_all_%d_10km_2_100km.rda", exposure))
      extinct_sp<-sp_dis_all%>%dplyr::filter(year==2100)
      extinct_sp<-extinct_sp%>%dplyr::filter(N_type=="EXTINCT")
      extinct_sp<-extinct_sp%>%dplyr::filter(M==dispersal)
      extinct_sp<-extinct_sp%>%dplyr::filter(TYPE==sprintf("Diversity_exposure_%d_dispersal_%d_10km_2_100km", 
                                                           exposure, dispersal))
      saveRDS(extinct_sp, sprintf("../../Objects/when_where_extinction_exposure_%d/extinct_sp_%d_10km_2_100km.rda", exposure, dispersal))
    }
  }
  
}
i=1
args = commandArgs(trailingOnly=TRUE)
g<-args[1]
if (is.na(g)){
  g<-"Mammals"
}
exposure<-as.numeric(args[2])
if (is.na(exposure)){
  exposure<-0
}



source("commonFuns/functions.r")
df<-NULL
for (dispersal in (c(0:1))){
  extinct_sp<-readRDS(sprintf("../../Objects/when_where_extinction_exposure_%d/extinct_sp_%d_10km_2_100km.rda", 
                              exposure, dispersal))
  extinct_sp<-extinct_sp%>%filter(group==g)
  for (i in c(1:nrow(extinct_sp))){
    print(paste(i, nrow(extinct_sp), g, "exposure=", exposure, "dispersal=", dispersal, sep=" - "))
    item<-extinct_sp[i,]
    target_folder<-sprintf("../../Objects/Dispersal/%s/%s", g, item$sp)
    
    #item$sp<-"Bunomys_fratrorum"
    #item$group<-"Mammals"
    init_disp_str<-sprintf("%s/initial_disp_10km_exposure_%d_dispersal_%d.rda", 
                           target_folder, exposure, dispersal)
    dispersal_log_str<-sprintf("%s/%s_%s_%d_dispersal_%d_10km.rda", target_folder, item$GCM, item$SSP, 
                                   exposure, dispersal)
    if ((!file.exists(init_disp_str))|(!file.exists(dispersal_log_str))){
      init_disp_str<-sprintf("%s/initial_disp_exposure_%d_dispersal_%d.rda", 
                             target_folder, exposure, dispersal)
      dispersal_log_str<-sprintf("%s/%s_%s_%d_dispersal_%d.rda", target_folder, item$GCM, item$SSP, 
                                     exposure, dispersal)
    }
    if (!file.exists(init_disp_str)){
      next()
    }
    st_dis<-readRDS(init_disp_str)
    colnames(st_dis)<-c("x", "y", "mask")
    #start_dis<-start_dis%>%dplyr::filter(in_out==1)
    st_dis<-st_dis%>%ungroup()%>%dplyr::distinct(x, y, mask)
    st_dis$YEAR<-2020
    
    dispersal_log<-readRDS(dispersal_log_str)
    extinct_year<-0
    if (is.null(dispersal_log)){
      extinct_year<-2021
    }
    if (length(dispersal_log)==0){
      extinct_year<-2021
    }
    if (extinct_year==0){
      dispersal_log<-rbindlist(dispersal_log)
      extinct_year<-max(dispersal_log$YEAR)+1
    }
    
    st_dis$group<-item$group
    st_dis$sp<-item$sp
    st_dis$GCM<-item$GCM
    st_dis$SSP<-item$SSP
    st_dis$extinct_year<-extinct_year
    st_dis$dispersal<-dispersal
    df<-bind_dplyr(df, st_dis)
  }
}
saveRDS(df, sprintf("../../Objects/when_where_extinction_exposure_%d/%s_10km.rda", exposure, g))
