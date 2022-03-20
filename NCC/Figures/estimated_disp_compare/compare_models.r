
library(ggplot2)
library(randomForest)
library(caret)
library(rgdal)
library(data.table)
library(ggpubr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

predicted_all<-readRDS("../../Objects/estimate_disp_dist/models/predicted_birds.rda")

evaluated_metrics_all<-readRDS("../../Objects/estimate_disp_dist/models/evaluated_metrics_birds.rda")
evaluated_metrics_all<-data.table(evaluated_metrics_all)
mod<-"RF"
birds_trait<-read.table("../../Data/Dispersal_distance/EltonTraits/BirdFuncDat.txt", sep="\t", head=T, stringsAsFactors = F)
ratio<-3.5

df_with_family<-data.table(read.csv("../../Data/Dispersal_distance/bird.csv", stringsAsFactors = F))

formulas<-c("log_max_dis~Diet",
            "log_max_dis~HWI",
            "log_max_dis~body_mass",
            "log_max_dis~HWI+body_mass",
            "log_max_dis~HWI+Diet",
            "log_max_dis~body_mass+Diet",
            "log_max_dis~HWI+body_mass+Diet")
#formulas<-formulas[7]
cols<-c("max_dis", "HWI", "log_body_mass", "Diet", "sp", "Order.x", "Family.name")
model_df_birds<-df_with_family[, ..cols]
model_df_birds$log_max_dis<-log(model_df_birds$max_dis)
model_df_birds$body_mass<-exp(1)^model_df_birds$log_body_mass
bird_disp<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")


for (mod in unique(evaluated_metrics_all$model)){
  if (mod=="CART"){
    mod_lab<-"rpart"
    
  }else{
    mod_lab<-mod
  }
  print(mod)
  evaluated_metrics_all_item<-evaluated_metrics_all[model==mod]
  best_model_info<-evaluated_metrics_all_item[RMSE_Full==min(RMSE_Full)]
  best_model<-readRDS(sprintf("../../Objects/estimate_disp_dist/models/%s_no_rank_birds.rda", 
                              tolower(mod_lab)))
  #best_model<-readRDS(sprintf("../../Objects/estimate_disp_dist/models/%s_no_rank.rda", "rf"))
  
  best_model<-best_model[[best_model_info$formulas]]
  cols<-c("HWI", "log_body_mass", "Diet", "Migration_3", "iucn_name")
  new_df_birds<-bird_disp[,..cols]
  new_df_birds$body_mass<-exp(1)^new_df_birds$log_body_mass
  
  new_df_birds$estimated_disp<-exp(1)^predict(best_model, new_df_birds)
  #new_df_birds[estimated_disp>dis_threshold]$estimated_disp<-new_df_birds[estimated_disp>dis_threshold]$estimated_disp*ratio
  new_df_birds$estimated_disp<-new_df_birds$estimated_disp*ratio
  
  unique(new_df_birds$Diet)
  table(new_df_birds$Migration_3)
  p<-ggplot(new_df_birds)+geom_boxplot(aes(y=estimated_disp, color=factor(Migration_3)))+theme_bw()
  p
  ggsave(p, filename=sprintf("../../Figures/Estimate_Disp/Predicted_birds%s.png", mod), width=8, height=6)
  
  
  #model_df_birds$estimated_disp<-predict(best_model, model_df_birds)
  
  
  p1<-ggplot(new_df_birds)+geom_point(aes(x=HWI, y=estimated_disp, color=factor(Diet)))+
    geom_vline(xintercept = min_HWI, color="black", linetype=2)+
    geom_vline(xintercept = max_HWI, color="black", linetype=2)+
    xlab("Hand-wing index (HWI)")+
    ylab("Estimated natal dispersal distance")+
    labs(color="Diet")+
    theme_bw()
  p1
  model_df_birds$body_mass2<-10^model_df_birds$log_body_mass
  new_df_birds$body_mass2<-10^new_df_birds$log_body_mass
  
  min_mass<-min(model_df_birds$body_mass2)
  max_mass<-max(model_df_birds$body_mass2)
  
  
  p2<-ggplot(new_df_birds)+geom_point(aes(x=body_mass2, y=estimated_disp, color=factor(Diet)))+
    geom_vline(xintercept = min_mass, color="black", linetype=2)+
    geom_vline(xintercept = max_mass, color="black", linetype=2)+
    scale_x_log10()+
    xlab("Body mass (gram, log transformed)")+
    ylab("Estimated natal dispersal distance")+
    labs(color="Diet")+
    theme_bw()
  p2
  
  p<-ggpubr::ggarrange(p2, p1, nrow=2)
  ggsave(p, filename=sprintf("../../Figures/Estimate_Disp/emperical_predicted_birds_%s.png", mod), width=8, height=8)
}
