library(data.table)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
mammal_df<-readRDS("../../Data/Mammals/mammal_df.rda")
bird_df<-readRDS("../../Data/Birds/bird_df.rda")
sp_list<-rbindlist(list(data.table(sp=unique(mammal_df$binomial), group="Mammals"), 
                        data.table(sp=unique(bird_df$SCINAME), group="Birds")))
sp_list$sp2<-gsub(" ", "_", sp_list$sp)
sp<-sp_list[1]
fitlist<-list()
fit1970list<-list()
fit100kmlist<-list()
fitebirdlist<-list()
fitseasonal1list<-list()
fitseasonal2list<-list()
for (i in 1:nrow(sp_list)){
  sp<-sp_list[i]
  print(paste(i, nrow(sp_list)))
  label<-paste(sp$sp, sp$group)
  
  f<-sprintf("../../Objects/Dispersal/%s/%s/fit.rda", sp$group, sp$sp2)
  if (file.exists(f)){
    fit<-readRDS(f)
    fit$sp<-sp$sp2
    fit$group<-sp$group
    fitlist[[label]]<-fit
  }
  
  f<-sprintf("../../Objects/Dispersal/%s/%s/fit_100km.rda", sp$group, sp$sp2)
  if (file.exists(f)){
    fit<-readRDS(f)
    fit$sp<-sp$sp2
    fit$group<-sp$group
    fit100kmlist[[label]]<-fit
  }
  
  f<-sprintf("../../Objects/Dispersal/%s/%s/fit_1970.rda", sp$group, sp$sp2)
  if (file.exists(f)){
    fit<-readRDS(f)
    fit$sp<-sp$sp2
    fit$group<-sp$group
    fit1970list[[label]]<-fit
  }
  
  f<-sprintf("../../Objects/Dispersal/%s/%s/fit_ebird.rda", sp$group, sp$sp2)
  if (file.exists(f)){
    fit<-readRDS(f)
    fit$sp<-sp$sp2
    fit$group<-sp$group
    fitebirdlist[[label]]<-fit
  }
  
  f<-sprintf("../../Objects/Dispersal/%s/%s/fit_seasonal_1.rda", sp$group, sp$sp2)
  if (file.exists(f)){
    fit<-readRDS(f)
    fit$sp<-sp$sp2
    fit$group<-sp$group
    fitseasonal1list[[label]]<-fit
  }
  
  f<-sprintf("../../Objects/Dispersal/%s/%s/fit_seasonal_2.rda", sp$group, sp$sp2)
  if (file.exists(f)){
    fit<-readRDS(f)
    fit$sp<-sp$sp2
    fit$group<-sp$group
    fitseasonal2list[[label]]<-fit
  }
  
  
}

fitlist<-rbindlist(fitlist)
saveRDS(fitlist, "../../Objects/niches/fitlist.rda")
fit1970list<-rbindlist(fit1970list)
saveRDS(fit1970list, "../../Objects/niches/fit1970list.rda")
fit100kmlist<-rbindlist(fit100kmlist)
saveRDS(fit100kmlist, "../../Objects/niches/fit100kmlist.rda")
fitebirdlist<-rbindlist(fitebirdlist, fill=T)
saveRDS(fitebirdlist, "../../Objects/niches/fitebirdlist.rda")
fitseasonal1list<-rbindlist(fitseasonal1list)
saveRDS(fitseasonal1list, "../../Objects/niches/fitseasonal1list.rda")
fitseasonal2list<-rbindlist(fitseasonal2list)
saveRDS(fitseasonal2list, "../../Objects/niches/fitseasonal2list.rda")
