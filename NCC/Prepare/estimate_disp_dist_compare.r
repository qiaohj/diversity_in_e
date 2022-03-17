library(ggplot2)
library(caret)
library(data.table)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

if (F){
  predicted_all<-readRDS("../../Objects/estimate_disp_dist/models/predicted_birds.rda")
  
  evaluated_metrics_all<-readRDS("../../Objects/estimate_disp_dist/models/evaluated_metrics_birds.rda")
  best_model_info_1<-evaluated_metrics_all[which(evaluated_metrics_all$RMSE_Full==min(evaluated_metrics_all$RMSE_Full)),]
  evaluated_metrics_all_2<-evaluated_metrics_all[which(evaluated_metrics_all$model!="RF"),]
  best_model_info_2<-evaluated_metrics_all_2[which(evaluated_metrics_all_2$RMSE_Full==min(evaluated_metrics_all_2$RMSE_Full)),]
  
  
  predicted_all_1<-predicted_all[model==best_model_info_1$model&formulas==best_model_info_1$formulas]
  predicted_all_2<-predicted_all[model==best_model_info_2$model&formulas==best_model_info_2$formulas]
  
  plot(predicted_all_1$pred_max_dis, predicted_all_2$pred_max_dis)
  cor(predicted_all_1$pred_max_dis, predicted_all_2$pred_max_dis)
  
  predicted_emperical_bird<-predicted_all_1
  predicted_emperical_bird$pred_max_dis_2<-predicted_all_2$pred_max_dis
  new_df_birds<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
  
  best_model_2<-readRDS(sprintf("../../Objects/estimate_disp_dist/models/%s_no_rank_birds.rda", 
                                tolower(best_model_info_2$model)))
  best_model_2<-best_model_2[[best_model_info_2$formulas]]
  
  new_df_birds$estimated_disp_2<-exp(1)^predict(best_model_2, new_df_birds)*3.5
  cor(new_df_birds$estimated_disp, new_df_birds$estimated_disp_2)
  plot(new_df_birds$estimated_disp, new_df_birds$estimated_disp_2)
  
  #mammals
  predicted_all<-readRDS("../../Objects/estimate_disp_dist/models/predicted_mammals.rda")
  
  evaluated_metrics_all<-readRDS("../../Objects/estimate_disp_dist/models/evaluated_metrics_mammals.rda")
  best_model_info_1<-evaluated_metrics_all[which(evaluated_metrics_all$RMSE_Full==min(evaluated_metrics_all$RMSE_Full)),]
  evaluated_metrics_all_2<-evaluated_metrics_all[which(evaluated_metrics_all$model!="RF"),]
  best_model_info_2<-evaluated_metrics_all_2[which(evaluated_metrics_all_2$RMSE_Full==min(evaluated_metrics_all_2$RMSE_Full)),]
  
  
  predicted_all_1<-predicted_all[model==best_model_info_1$model&formulas==best_model_info_1$formulas]
  predicted_all_2<-predicted_all[model==best_model_info_2$model&formulas==best_model_info_2$formulas]
  
  plot(predicted_all_1$pred_max_dis, predicted_all_2$pred_max_dis)
  cor(predicted_all_1$pred_max_dis, predicted_all_2$pred_max_dis)
  
  predicted_emperical_mammal<-predicted_all_1
  predicted_emperical_mammal$pred_max_dis_2<-predicted_all_2$pred_max_dis
  new_df_mammals<-readRDS("../../Objects/estimate_disp_dist/estimate_disp_dist_mammal.rda")
  
  best_model_2<-readRDS(sprintf("../../Objects/estimate_disp_dist/models/%s_no_rank_mammals.rda", 
                                tolower(best_model_info_2$model)))
  best_model_2<-best_model_2[[best_model_info_2$formulas]]
  
  new_df_mammals$estimated_disp_2<-exp(1)^predict(best_model_2, new_df_mammals)*2.5
  cor(new_df_mammals$estimated_disp, new_df_mammals$estimated_disp_2)
  plot(new_df_mammals$estimated_disp, new_df_mammals$estimated_disp_2)
  
  new_df_birds$group<-"Birds"
  cols<-c("estimated_disp", "estimated_disp_2", "group")
  new_df_birds_sub<-new_df_birds[, ..cols]
  
  new_df_mammals$group<-"Mammals"
  cols<-c("estimated_disp", "estimated_disp_2", "group")
  new_df_mammals_sub<-new_df_mammals[, ..cols]
  
  new_df<-rbindlist(list(new_df_birds_sub, new_df_mammals_sub))
  saveRDS(new_df, "../../Figures/estimated_disp_compare/all.rda")
  
  predicted_emperical_mammal$group<-"Mammals"
  cols<-c("pred_max_dis", "pred_max_dis_2", "group")
  predicted_emperical_mammal_sub<-predicted_emperical_mammal[, ..cols]
  
  predicted_emperical_bird$group<-"Birds"
  cols<-c("pred_max_dis", "pred_max_dis_2", "group")
  predicted_emperical_bird_sub<-predicted_emperical_bird[, ..cols]
  
  predicted_emperical_df<-rbindlist(list(predicted_emperical_bird_sub, predicted_emperical_mammal_sub))
  saveRDS(predicted_emperical_df, "../../Figures/estimated_disp_compare/emperical.rda")
}

emperical<-readRDS("../../Figures/estimated_disp_compare/emperical.rda")
all<-readRDS("../../Figures/estimated_disp_compare/all.rda")
colnames(emperical)<-c("RF", "NEXT", "group")
emperical$type="Emperical data"
colnames(all)<-c("RF", "NEXT", "group")
all$type="Full data"
df<-rbindlist(list(all, emperical))
p1<-ggplot(df)+geom_point(aes(x=RF, y=NEXT), size=0.4)+geom_abline(linetype=2, color="grey")+
  xlim(0, 1500)+ylim(0, 2000)+xlab("Best model")+ylab("Next best model")+
  theme_bw()+facet_grid(group~type, scale="free")
  
p1

#Group	Emperical	Full
#Birds	0.7870865	0.6556399
#Mammals	0.8316905	0.8196118

ggsave(p1, filename="../../Figures/estimated_disp_compare/estimated_disp_compare.png", width=6, height=5)
