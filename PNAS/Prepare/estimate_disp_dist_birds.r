library(randomForest)
library(caret)
library(rgdal)
library(data.table)

setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

#Sutherland
Sutherland<-read.table("../../Data/Dispersal_distance/Bird_Mammals_Dispersal_Distance/Bird.csv", head=T, sep=",")
colnames(Sutherland)<-c("sp", "body_mass", "obs_type", "dispersal_distance_median", 
                        "dispersal_distance_max", "x1", "x2", "n", "source")

Sutherland$max_str<-trimws(gsub("/", "", gsub("f", "", gsub("m", "", Sutherland$dispersal_distance_max))))
Sutherland$max_float<-as.numeric(Sutherland$max_str)

Sutherland$median_str<-trimws(gsub("/", "", gsub("f", "", gsub("m", "", Sutherland$dispersal_distance_median))))
Sutherland$median_float<-as.numeric(Sutherland$median_str)

Sutherland<-data.table(Sutherland)
Sutherland_se<-Sutherland[, .(median_dis=max(median_float, na.rm=T),
                              max_dis=max(max_float, na.rm=T)), by=list(sp)]
Sutherland_se$source<-"Sutherland, G. D., et. al., 2000"
Sutherland_se[!is.infinite(Sutherland_se$max_dis),]
Sutherland_se$ratio<-Sutherland_se$max_dis/Sutherland_se$median_dis

#summary(Sutherland_se$ratio)
#Sutherland_se<-Sutherland_se[!is.infinite(Sutherland_se$ratio),]
#Sutherland_se<-Sutherland_se[(Sutherland_se$ratio>0),]

plot(Sutherland_se$median_dis, Sutherland_se$max_dis)

birds_trait<-read.table("../../Data/Dispersal_distance/EltonTraits/Dataset_HWI.csv", sep=",", head=T, stringsAsFactors = F)
birds_trait[which(birds_trait$Diet=="fruit"), "Diet"]<-"plants"
birds_trait[which(birds_trait$Diet=="nectar"), "Diet"]<-"plants"
birds_trait[which(birds_trait$Diet=="scav"), "Diet"]<-"omnivore"
birds_trait[which(is.na(birds_trait$Diet)), "Diet"]<-"omnivore"

birds_trait[which(is.na(birds_trait$Migration_3)), "Migration_3"]<-"unknown"
birds_trait<-data.table(birds_trait)


unique(birds_trait$Diet)
iucn_bird<-read.csv("../../Data/Birds/HBW-BirdLife_List_of_Birds_v5.csv", head=T, stringsAsFactors = F)
iucn_bird<-data.table(iucn_bird)


df<-read.table("../../Data/Dispersal_distance/Bird_OPENBUGS/openbugs.csv", head=T, sep=",", stringsAsFactors = F)
df$meddisp<-as.numeric(df$meddisp)
df$mass<-as.numeric(df$mass)
df$wingspan<-as.numeric(df$wingspan)
df$shape<-(df$mass)^3/df$wingspan
df[which(df$sp=="Acrocephalus scirpaceous"), "sp"]<-"Acrocephalus scirpaceus"
df[which(df$sp=="Anas platyrhyncos"), "sp"]<-"Anas platyrhynchos"
df[which(df$sp=="Sita europea"), "sp"]<-"Sitta pygmaea"
df[which(df$sp=="Parus palustris"), "sp"]<-"Poecile palustris"
df[which(df$sp=="Parus ater"), "sp"]<-"Periparus ater"
df[which(df$sp=="Parus atricapillus"), "sp"]<-"Poecile atricapillus"
df[which(df$sp=="Parus caeruleus"), "sp"]<-"Cyanistes caeruleus"
df[which(df$sp=="Parus montanus"), "sp"]<-"Poecile montanus"
df[which(df$sp=="Carduelis choris"), "sp"]<-"Chloris chloris"
df[which(df$sp=="Carduelis cannabina"), "sp"]<-"Linaria cannabina"
df[which(df$sp=="Asio otis"), "sp"]<-"Asio otus"
df[which(df$sp=="Columba polumbus"), "sp"]<-"Columba palumbus"
df[which(df$sp=="Cygnus olur"), "sp"]<-"Cygnus olor"
df[which(df$sp=="Delichon urbica"), "sp"]<-"Delichon urbicum"
df[which(df$sp=="Emberiza citronella"), "sp"]<-"Emberiza citrinella"
df[which(df$sp=="Iridoprocne bicolor"), "sp"]<-"Tachycineta bicolor"
df[which(df$sp=="Otis asio"), "sp"]<-"Otus alius"

df_item_1<-df[, c("sp", "meddisp")]
colnames(df_item_1)[2]<-"median_dis"
df_item_1$max_dis<-NA
df_item_2<-Sutherland_se[, c("sp", "median_dis", "max_dis")]
df_item_2[which(is.infinite(df_item_2$median_dis)), "median_dis"]<-NA
df_item_2[which(is.infinite(df_item_2$max_dis)), "max_dis"]<-NA

#PARADIS
Paradis<-read.table("../../Data/Dispersal_distance/Bird-1998/Data0608.csv", head=T,
                    sep=",", stringsAsFactors = F)

Paradis$max_dis<-NA

#Others
df_item_3<-read.table("../../Data/Dispersal_distance/Bird_Others/others.csv", head=T,
                      sep=",", stringsAsFactors = F)
#df_item_3$max_dis<-NA
#df_item<-rbind(df_item_1, df_item_2)
#df_item<-rbind(rbind(df_item_1, df_item_2), Paradis)
#df_item<-rbind(df_item_1, df_item_3)
df_item<-rbind(df_item_2, df_item_3)
df_item<-df_item[which(!is.na(df_item$max_dis)),]
df_item[which(df_item$sp=="Dryobates borealis"), "sp"]<-"Leuconotopicus borealis"
df_item[which(df_item$sp=="Dryocopus pileatus"), "sp"]<-"Dryocopus pileatus"
df_item[which(df_item$sp=="Phalacrocorax pelagicus"), "sp"]<-"Phalacrocorax pelagicus"


ratio<-df_item$max_dis/df_item$median_dis
xxx<-df_item[!is.infinite(df_item$median_dis),]
plot(xxx$max_dis, xxx$median_dis)
ratio<-ratio[!is.na(ratio) & !is.infinite(ratio)]
summary(ratio)

ratio<-quantile(ratio, 0.5)
df_with_family<-merge(df_item, iucn_bird, by.x="sp", by.y="Scientific.name", all.x=T, all.y=F)
df_with_family[is.na(df_with_family$Order),]


df_with_family<-df_with_family[!is.na(df_with_family$Order),]
#df_with_family[is.na(df_with_family$median_dis), "median_dis"]<-
#  df_with_family[is.na(df_with_family$median_dis), "max_dis"]/ratio
#View(df_with_family[is.na(Family.name)])
df_with_family<-data.table(df_with_family)
df_with_family<-df_with_family[, .(max_dis = max(max_dis)), 
                               by=list(sp, Order, Family.name)]
df_with_family<-merge(df_with_family, birds_trait, by.x="sp", by.y="iucn_name", all.x=T, all.y=F)
df_with_family<-df_with_family[!is.na(df_with_family$HWI)]

df_with_family$Diet<-as.factor(df_with_family$Diet)

formulas<-c("max_dis~Diet",
            "max_dis~HWI",
            "max_dis~log_body_mass",
            "max_dis~HWI+log_body_mass",
            "max_dis~HWI+Diet",
            "max_dis~log_body_mass+Diet",
            "max_dis~HWI+log_body_mass+Diet")
#formulas<-formulas[7]
cols<-c("max_dis", "HWI", "log_body_mass", "Diet")
model_df_birds<-df_with_family[, ..cols]
dim(model_df_birds)
table(model_df_birds$Diet)

if (T){
  ########For random forest##########################
  # define training control for 10-fold Cross Validation
  evaluated_metrics_all<-NULL
  predicted_all<-NULL
  best_r2<-function(x, metric, maximize){
    #print(x)
    best(x, "Rsquared", maximize = T)
  }
  seeds_index<-1000
  setSeeds <- function(method = "cv", numbers = 1, repeats = 1, tunes = NULL, seed = 1) {
    #B is the number of resamples and integer vector of M (numbers + tune length if any)
    B <- if (method == "cv") numbers
    else if(method == "repeatedcv") numbers * repeats
    else NULL
    
    if(is.null(length)) {
      seeds <- NULL
    } else {
      set.seed(seed = seed)
      seeds <- vector(mode = "list", length = B)
      seeds <- lapply(seeds, function(x) sample.int(n = 1000000, size = numbers + ifelse(is.null(tunes), 0, tunes)))
      seeds[[length(seeds) + 1]] <- sample.int(n = 1000000, size = 1)
    }
    # return seeds
    seeds
  }
  
  
  i=7
  rf_no_rank_list<-list()
  glm_no_rank_list<-list()
  svm_no_rank_list<-list()
  brnn_no_rank_list<-list()
  rpart_no_rank_list<-list()
  evtree_no_rank_list<-list()
  
  for (i in c(1:length(formulas))){
    f<-formulas[i]
    print(f)
    # train the model
    #df_with_family$Migration_3<-as.factor(df_with_family$Migration_3)
    tunegrid <- expand.grid(.mtry=c(1:5))
    set.seed(seeds_index)
    seeds <- setSeeds("repeatedcv", numbers=10, repeats=10, tunes=5, seed=seeds_index)
    train_control <- trainControl(method="repeatedcv", number=10, repeats=10, 
                                  seeds=seeds)
    rf_no_rank <- train(as.formula(f), data=model_df_birds, 
                        trControl=train_control, method="rf", 
                        tuneGrid=tunegrid, ntree=1000)
    
    rf_no_rank_list[[f]]<-rf_no_rank
    # summarize results
    print(rf_no_rank)
    predicted<-model_df_birds
    predicted$model<-"RF"
    predicted$formulas<-f
    predicted$clade<-"Bird"
    predicted$pred<-predict(rf_no_rank, model_df_birds)
    cor<-cor(model_df_birds$max_dis, predicted$pred)
    
    evaluated_metrics<-rf_no_rank$results[which(rf_no_rank$results$RMSE==min(rf_no_rank$results$RMSE)),]
    evaluated_metrics$model<-"RF"
    evaluated_metrics$formulas<-f
    evaluated_metrics$clade<-"Bird"
    evaluated_metrics$RMSE_Full<-RMSE(predicted$pred, model_df_birds$max_dis, na.rm=T)
    evaluated_metrics$Rsquared_Full<-R2(predicted$pred, model_df_birds$max_dis, na.rm=T)
    evaluated_metrics$MAE_Full<-MAE(predicted$pred, model_df_birds$max_dis, na.rm=T)
    evaluated_metrics$cor<-cor
    evaluated_metrics<-evaluated_metrics[, c("RMSE", "Rsquared", "MAE", 
                                             "RMSE_Full", "Rsquared_Full", "MAE_Full",
                                             "cor", "model", "formulas", "clade")]
    if (is.null(evaluated_metrics_all)){
      evaluated_metrics_all<-evaluated_metrics  
      predicted_all<-predicted
    }else{
      evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)
      predicted_all<-rbind(predicted_all, predicted)
    }
    
    #######################end of random forest ################################
    
    
    ########For GLM##########################
    # train the model
    set.seed(seeds_index)
    seeds <- setSeeds("repeatedcv", numbers=10, repeats=10, tunes=5, seed=seeds_index)
    train_control <- trainControl(method="repeatedcv", number=10, repeats=10, 
                                  seeds=seeds)
    glm_no_rank <- train(as.formula(f), 
                         data=model_df_birds, trControl=train_control, 
                         method="glm")
    glm_no_rank_list[[f]]<-glm_no_rank
    # summarize results
    print(glm_no_rank)
    predicted<-model_df_birds
    predicted$model<-"GLM"
    predicted$formulas<-f
    predicted$clade<-"Bird"
    predicted$pred<-predict(glm_no_rank, model_df_birds)
    cor<-cor(model_df_birds$max_dis, predicted$pred)
    
    evaluated_metrics<-glm_no_rank$results[which(glm_no_rank$results$RMSE==min(glm_no_rank$results$RMSE)),]
    evaluated_metrics$model<-"GLM"
    evaluated_metrics$formulas<-f
    evaluated_metrics$clade<-"Bird"
    evaluated_metrics$RMSE_Full<-RMSE(predicted$pred, model_df_birds$max_dis, na.rm=T)
    evaluated_metrics$Rsquared_Full<-R2(predicted$pred, model_df_birds$max_dis, na.rm=T)
    evaluated_metrics$MAE_Full<-MAE(predicted$pred, model_df_birds$max_dis, na.rm=T)
    evaluated_metrics$cor<-cor
    evaluated_metrics<-evaluated_metrics[, c("RMSE", "Rsquared", "MAE", 
                                             "RMSE_Full", "Rsquared_Full", "MAE_Full",
                                             "cor", "model", "formulas", "clade")]
    evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)
    predicted_all<-rbind(predicted_all, predicted)
    #######################end of GLM ################################
    
    ########For Support Vector Machines with Linear Kernel ##########################
    # train the model
    set.seed(seeds_index)
    seeds <- setSeeds("repeatedcv", numbers=10, repeats=10, tunes=5, seed=seeds_index)
    train_control <- trainControl(method="repeatedcv", number=10, repeats=10, 
                                  seeds=seeds)
    svm_no_rank <- train(as.formula(f), data=model_df_birds, trControl=train_control, 
                         method="svmLinear")
    svm_no_rank_list[[f]]<-svm_no_rank
    # summarize results
    print(svm_no_rank)
    predicted<-model_df_birds
    predicted$model<-"SVM"
    predicted$formulas<-f
    predicted$clade<-"Bird"
    predicted$pred<-predict(svm_no_rank, model_df_birds)
    cor<-cor(model_df_birds$max_dis, predicted$pred)
    
    evaluated_metrics<-svm_no_rank$results[which(svm_no_rank$results$RMSE==min(svm_no_rank$results$RMSE)),]
    evaluated_metrics$model<-"SVM"
    evaluated_metrics$formulas<-f
    evaluated_metrics$clade<-"Bird"
    evaluated_metrics$RMSE_Full<-RMSE(predicted$pred, model_df_birds$max_dis, na.rm=T)
    evaluated_metrics$Rsquared_Full<-R2(predicted$pred, model_df_birds$max_dis, na.rm=T)
    evaluated_metrics$MAE_Full<-MAE(predicted$pred, model_df_birds$max_dis, na.rm=T)
    evaluated_metrics$cor<-cor
    evaluated_metrics<-evaluated_metrics[, c("RMSE", "Rsquared", "MAE", 
                                             "RMSE_Full", "Rsquared_Full", "MAE_Full",
                                             "cor", "model", "formulas", "clade")]
    evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)
    predicted_all<-rbind(predicted_all, predicted)
    #######################end of SVM ################################
    
    
    
    ########For Bayesian Regularized Neural Networks##########################
    # train the model
    set.seed(seeds_index)
    seeds <- setSeeds("repeatedcv", numbers=10, repeats=10, tunes=5, seed=seeds_index)
    train_control <- trainControl(method="repeatedcv", number=10, repeats=10, 
                                  seeds=seeds)
    brnn_no_rank <- train(as.formula(f), data=model_df_birds, trControl=train_control, 
                          method="brnn")
    brnn_no_rank_list[[f]]<-brnn_no_rank
    # summarize results
    print(brnn_no_rank)
    predicted<-model_df_birds
    predicted$model<-"BRNN"
    predicted$formulas<-f
    predicted$clade<-"Bird"
    predicted$pred<-predict(brnn_no_rank, model_df_birds)
    cor<-cor(model_df_birds$max_dis, predicted$pred)
    
    evaluated_metrics<-brnn_no_rank$results[which(brnn_no_rank$results$RMSE==min(brnn_no_rank$results$RMSE)),]
    evaluated_metrics$model<-"BRNN"
    evaluated_metrics$formulas<-f
    evaluated_metrics$clade<-"Bird"
    evaluated_metrics$RMSE_Full<-RMSE(predicted$pred, model_df_birds$max_dis, na.rm=T)
    evaluated_metrics$Rsquared_Full<-R2(predicted$pred, model_df_birds$max_dis, na.rm=T)
    evaluated_metrics$MAE_Full<-MAE(predicted$pred, model_df_birds$max_dis, na.rm=T)
    evaluated_metrics$cor<-cor
    evaluated_metrics<-evaluated_metrics[, c("RMSE", "Rsquared", "MAE", 
                                             "RMSE_Full", "Rsquared_Full", "MAE_Full",
                                             "cor", "model", "formulas", "clade")]
    evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)
    predicted_all<-rbind(predicted_all, predicted)
    #######################end of BRNN ################################
    
    ########For CART##########################
    # train the model
    set.seed(seeds_index)
    seeds <- setSeeds("repeatedcv", numbers=10, repeats=10, tunes=5, seed=seeds_index)
    train_control <- trainControl(method="repeatedcv", number=10, repeats=10, 
                                  seeds=seeds)
    rpart_no_rank <- train(as.formula(f), data=model_df_birds, trControl=train_control, 
                          method="rpart")
    rpart_no_rank_list[[f]]<-rpart_no_rank
    # summarize results
    print(rpart_no_rank)
    predicted<-model_df_birds
    predicted$model<-"CART"
    predicted$formulas<-f
    predicted$clade<-"Bird"
    predicted$pred<-predict(rpart_no_rank, model_df_birds)
    cor<-cor(model_df_birds$max_dis, predicted$pred)
    
    evaluated_metrics<-rpart_no_rank$results[which(rpart_no_rank$results$RMSE==min(rpart_no_rank$results$RMSE)),]
    evaluated_metrics$model<-"CART"
    evaluated_metrics$formulas<-f
    evaluated_metrics$clade<-"Bird"
    evaluated_metrics$RMSE_Full<-RMSE(predicted$pred, model_df_birds$max_dis, na.rm=T)
    evaluated_metrics$Rsquared_Full<-R2(predicted$pred, model_df_birds$max_dis, na.rm=T)
    evaluated_metrics$MAE_Full<-MAE(predicted$pred, model_df_birds$max_dis, na.rm=T)
    evaluated_metrics$cor<-cor
    evaluated_metrics<-evaluated_metrics[, c("RMSE", "Rsquared", "MAE", 
                                             "RMSE_Full", "Rsquared_Full", "MAE_Full",
                                             "cor", "model", "formulas", "clade")]
    evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)
    predicted_all<-rbind(predicted_all, predicted)
    #######################end of CART ################################
    
  }
  saveRDS(rf_no_rank_list, "../../Objects/estimate_disp_dist/models/rf_no_rank_birds.rda")
  
  saveRDS(glm_no_rank_list, "../../Objects/estimate_disp_dist/models/glm_no_rank_birds.rda")
  
  saveRDS(svm_no_rank_list, "../../Objects/estimate_disp_dist/models/svm_no_rank_birds.rda")
  
  saveRDS(brnn_no_rank_list, "../../Objects/estimate_disp_dist/models/brnn_no_rank_birds.rda")
  
  saveRDS(rpart_no_rank_list, "../../Objects/estimate_disp_dist/models/rpart_no_rank_birds.rda")
  
  saveRDS(evtree_no_rank_list, "../../Objects/estimate_disp_dist/models/evtree_no_rank_birds.rda")
  
  
  saveRDS(evaluated_metrics_all, "../../Objects/estimate_disp_dist/models/evaluated_metrics_birds.rda")
  
  saveRDS(predicted_all, "../../Objects/estimate_disp_dist/models/predicted_birds.rda")
  
}

library(ggplot2)
predicted_all<-readRDS("../../Objects/estimate_disp_dist/models/predicted_birds.rda")

evaluated_metrics_all<-readRDS("../../Objects/estimate_disp_dist/models/evaluated_metrics_birds.rda")

evaluated_metrics_all[which(evaluated_metrics_all$RMSE_Full==min(evaluated_metrics_all$RMSE_Full)),]
if (F){
  evaluated_metrics_all<-evaluated_metrics_all[which(!((evaluated_metrics_all$model=="RF")&
                                                         (evaluated_metrics_all$formulas=="max_dis~HWI"))),]
  evaluated_metrics_all<-evaluated_metrics_all[which(!((evaluated_metrics_all$model=="RF")&
                                                         (evaluated_metrics_all$formulas=="max_dis~log_body_mass"))),]
  evaluated_metrics_all<-evaluated_metrics_all[which(!((evaluated_metrics_all$model=="RF")&
                                                         (evaluated_metrics_all$formulas=="max_dis~HWI+log_body_mass"))),]
}
evaluated_metrics_all[which(evaluated_metrics_all$cor==max(evaluated_metrics_all$cor, na.rm = T)),]

#RMSE_Full
p<-ggplot(predicted_all)+
  geom_smooth(aes(x=max_dis, y=pred, 
                  color=factor(model), linetype=factor(formulas)),
                  alpha=0.2)+
  geom_point(aes(x=max_dis, y=pred, 
                 color=factor(model), shape=factor(formulas)))+
  theme_bw()
p
ggsave(p, filename="../../Figures/Estimate_Disp/Model_Results_birds.png", width=8, height=6)

p<-ggplot(evaluated_metrics_all)+
  geom_point(aes(x=factor(model), y=RMSE, color=factor(formulas), size=cor))+
  theme_bw()

p
ggsave(p, filename="../../Figures/Estimate_Disp/RMSE_birds.png", width=8, height=6)

p<-ggplot(evaluated_metrics_all)+
  geom_point(aes(x=factor(model), y=Rsquared, color=factor(formulas), size=cor))+
  theme_bw()

p
ggsave(p, filename="../../Figures/Estimate_Disp/Rsquared_birds.png", width=8, height=6)

p<-ggplot(evaluated_metrics_all)+
  geom_point(aes(x=factor(model), y=MAE, color=factor(formulas), size=cor))+
  theme_bw()

p
ggsave(p, filename="../../Figures/Estimate_Disp/MAE_birds.png", width=8, height=6)

plot(evaluated_metrics_all$Rsquared, evaluated_metrics_all$cor)
threshold_rsqrt<-quantile(evaluated_metrics_all$Rsquared, 0.9)
threshold_rmse<-quantile(evaluated_metrics_all$RMSE, 0.1)
best_model_info<-
  evaluated_metrics_all[
    which(evaluated_metrics_all$RMSE<=threshold_rmse),]


#best_model_info<-
#  evaluated_metrics_all[
#    which(evaluated_metrics_all$Rsquared>=threshold_rsqrt),]





best_model_info<-best_model_info[which(best_model_info$RMSE==min(best_model_info$RMSE)),]
best_model<-readRDS(sprintf("../../Objects/estimate_disp_dist/models/%s_no_rank_birds.rda", tolower(best_model_info$model)))
#best_model<-readRDS(sprintf("../../Objects/estimate_disp_dist/models/%s_no_rank.rda", "rf"))

best_model<-best_model[[best_model_info$formulas]]
cols<-c("HWI", "log_body_mass", "Diet", "Migration_3", "iucn_name")
new_df_birds<-birds_trait[,..cols]
new_df_birds[is.na(Migration_3)]
new_df_birds<-new_df_birds[(!is.na(HWI))&(!is.na(log_body_mass))&(!is.na(Diet))]
new_df_birds$estimated_disp<-predict(best_model, new_df_birds)

unique(new_df_birds$Diet)
table(new_df_birds$Migration_3)
p<-ggplot(new_df_birds)+geom_boxplot(aes(y=estimated_disp, color=factor(Migration_3)))+theme_bw()
p
ggsave(p, filename="../../Figures/Estimate_Disp/Predicted_birds.png", width=8, height=6)

new_df_birds

saveRDS(new_df_birds, "../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")

model_df_birds$estimated_disp<-predict(best_model, model_df_birds[, "HWI"])

p1<-ggplot(model_df_birds)+geom_point(aes(x=max_dis, y=estimated_disp, color=factor(Diet)))+theme_bw()
p1
p2<-ggplot(model_df_birds)+geom_point(aes(x=HWI, y=max_dis, color=factor(Diet)))+theme_bw()
p2
p<-ggpubr::ggarrange(p2, p1, nrow=2)
ggsave(p, filename="../../Figures/Estimate_Disp/emperical_predicted_birds.png", width=8, height=8)

