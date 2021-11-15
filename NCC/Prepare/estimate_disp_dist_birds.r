library(randomForest)
library(caret)
library(rgdal)
library(data.table)
library(ggpubr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")

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

dim(model_df_birds)
table(model_df_birds$Diet)

plot(log(model_df_birds$max_dis), model_df_birds$log_body_mass)
plot(log(model_df_birds$max_dis), model_df_birds$HWI)
plot(log(model_df_birds$max_dis), log(model_df_birds$HWI))
plot(model_df_birds$max_dis, model_df_birds$HWI)


model_df_birds<-model_df_birds
model<-randomForest(log_max_dis~HWI+log_body_mass+Diet, model_df_birds, ntree=1000, mtry=3)
model_df_birds$pred_values = exp(1)^predict(model, model_df_birds)
mean((model_df_birds$pred_values-model_df_birds$max_dis)^2)^0.5
p1<-ggplot(model_df_birds)+geom_point(aes(x=max_dis, y=pred_values, color=Diet))+
  xlim(min(c(model_df_birds$max_dis, model_df_birds$pred_values)), max(c(model_df_birds$max_dis, model_df_birds$pred_values)))+
  ylim(min(c(model_df_birds$max_dis, model_df_birds$pred_values)), max(c(model_df_birds$max_dis, model_df_birds$pred_values)))+
  geom_abline()+theme_bw()
p1

model<-randomForest(max_dis~HWI+log_body_mass+Diet, model_df_birds, ntree=1000, mtry=3)
model_df_birds$pred_values = predict(model, model_df_birds)
mean((model_df_birds$pred_values-model_df_birds$max_dis)^2)^0.5
p2<-ggplot(model_df_birds)+geom_point(aes(x=max_dis, y=pred_values, color=Diet))+
  xlim(min(c(model_df_birds$max_dis, model_df_birds$pred_values)), max(c(model_df_birds$max_dis, model_df_birds$pred_values)))+
  ylim(min(c(model_df_birds$max_dis, model_df_birds$pred_values)), max(c(model_df_birds$max_dis, model_df_birds$pred_values)))+
  geom_abline()+theme_bw()
p2


model<-randomForest(max_dis~HWI+body_mass+Diet, model_df_birds, ntree=1000, mtry=3)
model_df_birds$pred_values = predict(model, model_df_birds)
mean((model_df_birds$pred_values-model_df_birds$max_dis)^2)^0.5
p3<-ggplot(model_df_birds)+geom_point(aes(x=max_dis, y=pred_values, color=Diet))+
  xlim(min(c(model_df_birds$max_dis, model_df_birds$pred_values)), max(c(model_df_birds$max_dis, model_df_birds$pred_values)))+
  ylim(min(c(model_df_birds$max_dis, model_df_birds$pred_values)), max(c(model_df_birds$max_dis, model_df_birds$pred_values)))+
  geom_abline()+theme_bw()
p3

model_df<-model_df_birds[max_dis<500]
model<-randomForest(log_max_dis~HWI+body_mass+Diet, model_df, ntree=1000, mtry=3)
model_df$pred_values = exp(1)^predict(model, model_df)
model_df$time<-model_df$max_dis/model_df$pred_values
dis_threshold<-200
ratio<-mean(model_df[(max_dis>dis_threshold)&(Diet!="vertebrates")]$max_dis/
              model_df[(max_dis>dis_threshold)&(Diet!="vertebrates")]$pred_values)
model_df[(max_dis>dis_threshold)&(Diet!="vertebrates")]$pred_values<-
  model_df[(max_dis>dis_threshold)&(Diet!="vertebrates")]$pred_values*ratio
mean((model_df$pred_values-model_df$max_dis)^2)^0.5
p4<-ggplot(model_df)+geom_point(aes(x=max_dis, y=pred_values, color=Diet))+
  xlim(min(c(model_df$max_dis, model_df$pred_values)), max(c(model_df$max_dis, model_df$pred_values)))+
  ylim(min(c(model_df$max_dis, model_df$pred_values)), max(c(model_df$max_dis, model_df$pred_values)))+
  geom_abline()+theme_bw()
p4



model_df<-model_df_birds
model<-randomForest(log_max_dis~HWI+body_mass+Diet, model_df, ntree=1000, mtry=3)
model_df$pred_values = exp(1)^predict(model, model_df)
model_df$time<-model_df$max_dis/model_df$pred_values
ratio<-mean(model_df[(max_dis>dis_threshold)&(Diet!="vertebrates")]$max_dis/
              model_df[(max_dis>dis_threshold)&(Diet!="vertebrates")]$pred_values)
model_df[(max_dis>dis_threshold)&(Diet!="vertebrates")]$pred_values<-
  model_df[(max_dis>dis_threshold)&(Diet!="vertebrates")]$pred_values*ratio

mean((model_df$pred_values-model_df$max_dis)^2)^0.5
p4.2<-ggplot(model_df)+geom_point(aes(x=max_dis, y=pred_values, color=Diet))+
  xlim(min(c(model_df$max_dis, model_df$pred_values)), max(c(model_df$max_dis, model_df$pred_values)))+
  ylim(min(c(model_df$max_dis, model_df$pred_values)), max(c(model_df$max_dis, model_df$pred_values)))+
  geom_abline()+theme_bw()
p4.2

model_df<-model_df_birds[max_dis<500]
model<-randomForest(log_max_dis~HWI+body_mass+Diet, model_df, ntree=1000, mtry=3)
model_df$pred_values = exp(1)^predict(model, model_df)
model_df$time<-model_df$max_dis/model_df$pred_values

mean((model_df$pred_values-model_df$max_dis)^2)^0.5
p4.3<-ggplot(model_df)+geom_point(aes(x=max_dis, y=pred_values, color=Diet))+
  xlim(0, 100)+
  ylim(0, 100)+
  geom_abline()+theme_bw()
p4.3


ggarrange(p4, p4.2)

model_df<-model_df_birds[Diet=="vertebrates"]
model<-randomForest(log_max_dis~HWI+body_mass+Diet, model_df, ntree=1000, mtry=3)
model_df$pred_values = exp(1)^predict(model, model_df)
mean((model_df$pred_values-model_df$max_dis)^2)^0.5
p5<-ggplot(model_df)+geom_point(aes(x=max_dis, y=pred_values, color=Diet, size=body_mass))+
  xlim(min(c(model_df$max_dis, model_df$pred_values)), max(c(model_df$max_dis, model_df$pred_values)))+
  ylim(min(c(model_df$max_dis, model_df$pred_values)), max(c(model_df$max_dis, model_df$pred_values)))+
  geom_abline()+theme_bw()
p5

model_df<-model_df_birds[Diet=="invertebrates"]
model<-randomForest(log_max_dis~HWI, model_df, ntree=1000, mtry=3)
model_df$pred_values = exp(1)^predict(model, model_df)
mean((model_df$pred_values-model_df$max_dis)^2)^0.5
p6<-ggplot(model_df)+geom_point(aes(x=max_dis, y=pred_values, color=HWI, size=HWI))+
  xlim(min(c(model_df$max_dis, model_df$pred_values)), max(c(model_df$max_dis, model_df$pred_values)))+
  ylim(min(c(model_df$max_dis, model_df$pred_values)), max(c(model_df$max_dis, model_df$pred_values)))+
  geom_abline()+theme_bw()
p6

model_df<-model_df_birds[!(Order.x %in% c("PELECANIFORMES", "PICIFORMES",  "PSITTACIFORMES"))]

model<-randomForest(log_max_dis~HWI+body_mass+Diet, model_df, ntree=1000, mtry=3)
model_df$pred_values = exp(1)^predict(model, model_df)
mean((model_df$pred_values-model_df$max_dis)^2)^0.5
p7<-ggplot(model_df)+geom_point(aes(x=max_dis, y=pred_values, color=Order.x))+
  xlim(min(c(model_df$max_dis, model_df$pred_values)), max(c(model_df$max_dis, model_df$pred_values)))+
  ylim(min(c(model_df$max_dis, model_df$pred_values)), max(c(model_df$max_dis, model_df$pred_values)))+
  geom_abline()+theme_bw()
p7

library(mgcv)
library(e1071)
model_df_birds$xxx<-model_df_birds$body_mass/model_df_birds$HWI

model<-svm(log_max_dis~HWI+body_mass+Diet,data= model_df_birds, method="REML")
v<-as.numeric(exp(1)^predict(model, model_df_birds))
v[64:81]<-v[63:80]
v[63]<-0
model_df_birds$pred_values = v

mean((model_df_birds$pred_values-model_df_birds$max_dis)^2)^0.5
p8<-ggplot(model_df_birds)+geom_point(aes(x=max_dis, y=pred_values, color=Diet))+
  xlim(min(c(model_df_birds$max_dis, model_df_birds$pred_values)), max(c(model_df_birds$max_dis, model_df_birds$pred_values)))+
  ylim(min(c(model_df_birds$max_dis, model_df_birds$pred_values)), max(c(model_df_birds$max_dis, model_df_birds$pred_values)))+
  geom_abline()+theme_bw()
p8

ggarrange(p1, p2, nrow=2)


predicted_all<-readRDS("../../Objects_PNAS/estimate_disp_dist/models/predicted_birds.rda")
for (f in unique(predicted_all$formulas)){
  if (f %in% c("max_dis~Diet", "max_dis~HWI")){
    next()
  }
  for (m in unique(predicted_all$model)){
    print(paste(f, m, sep=" ///// "))
    item<-predicted_all[(formulas==f)&(model==m)]
    p<-ggplot(item)+geom_point(aes(x=max_dis, y=pred, color=Diet))+
      xlim(min(c(item$max_dis, item$pred)), max(c(item$max_dis, item$pred)))+
      ylim(min(c(item$max_dis, item$pred)), max(c(item$max_dis, item$pred)))+
      geom_abline()+theme_bw()
    print(p)
    xx <- readline(prompt="X to stop: ")
    if (tolower(xx)=='x') {
      break()
    }
  }
}

if (F){
  dis_threshold<-200
  ratio<-3.5
  formulas<-c("log_max_dis~Diet",
              "log_max_dis~HWI",
              "log_max_dis~body_mass",
              "log_max_dis~HWI+body_mass",
              "log_max_dis~HWI+Diet",
              "log_max_dis~body_mass+Diet",
              "log_max_dis~HWI+body_mass+Diet",
              "log_max_dis~log_body_mass",
              "log_max_dis~HWI+log_body_mass",
              "log_max_dis~log_body_mass+Diet",
              "log_max_dis~HWI+log_body_mass+Diet")
  
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
    tunegrid <- expand.grid(.mtry=c(3:5))
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
    predicted$pred_max_dis<-exp(1)^predicted$pred
    predicted[max_dis>200]$pred_max_dis<-predicted[max_dis>200]$pred_max_dis*ratio
    #plot(predicted$pred_max_dis, predicted$max_dis)
    cor<-cor(model_df_birds$max_dis, predicted$pred_max_dis)
    
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
    predicted$pred_max_dis<-exp(1)^predicted$pred
    predicted[max_dis>200]$pred_max_dis<-predicted[max_dis>200]$pred_max_dis*ratio
    cor<-cor(model_df_birds$max_dis, predicted$pred_max_dis)
    
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
    predicted$pred_max_dis<-exp(1)^predicted$pred
    predicted[max_dis>200]$pred_max_dis<-predicted[max_dis>200]$pred_max_dis*ratio
    cor<-cor(model_df_birds$max_dis, predicted$pred_max_dis)
    
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
    predicted$pred_max_dis<-exp(1)^predicted$pred
    predicted[max_dis>200]$pred_max_dis<-predicted[max_dis>200]$pred_max_dis*ratio
    cor<-cor(model_df_birds$max_dis, predicted$pred_max_dis)
    
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
    predicted$pred_max_dis<-exp(1)^predicted$pred
    predicted[max_dis>200]$pred_max_dis<-predicted[max_dis>200]$pred_max_dis*ratio
    cor<-cor(model_df_birds$max_dis, predicted$pred_max_dis)
    
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

write.csv(evaluated_metrics_all, "../../Figures/Estimate_Disp/evaluated_metrics_all_birds.csv")

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
  geom_smooth(aes(x=max_dis, y=pred_max_dis, 
                  color=factor(model), linetype=factor(formulas)),
                  alpha=0.2)+
  geom_point(aes(x=max_dis, y=pred_max_dis, 
                 color=factor(model), shape=factor(formulas)))+
  geom_abline()+
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



unique(model_df_birds$Diet)
best_model_info<-evaluated_metrics_all[which(evaluated_metrics_all$RMSE_Full==min(evaluated_metrics_all$RMSE_Full)),]
best_model<-readRDS(sprintf("../../Objects/estimate_disp_dist/models/%s_no_rank_birds.rda", tolower(best_model_info$model)))
#best_model<-readRDS(sprintf("../../Objects/estimate_disp_dist/models/%s_no_rank.rda", "rf"))

best_model<-best_model[[best_model_info$formulas]]
cols<-c("HWI", "log_body_mass", "Diet", "Migration_3", "iucn_name")
new_df_birds<-birds_trait[,..cols]
new_df_birds[is.na(Migration_3)]
new_df_birds<-new_df_birds[(!is.na(HWI))&(!is.na(log_body_mass))&(!is.na(Diet))]
new_df_birds[Diet=="fruit"]$Diet<-"plants"
new_df_birds[Diet=="seeds"]$Diet<-"plants"
new_df_birds[Diet=="nectar"]$Diet<-"plants"
new_df_birds[Diet=="scav"]$Diet<-"omnivore"
new_df_birds$body_mass<-exp(1)^new_df_birds$log_body_mass


new_df_birds$estimated_disp<-exp(1)^predict(best_model, new_df_birds)
#new_df_birds[estimated_disp>dis_threshold]$estimated_disp<-new_df_birds[estimated_disp>dis_threshold]$estimated_disp*ratio
new_df_birds$estimated_disp<-new_df_birds$estimated_disp*ratio

unique(new_df_birds$Diet)
table(new_df_birds$Migration_3)
p<-ggplot(new_df_birds)+geom_boxplot(aes(y=estimated_disp, color=factor(Migration_3)))+theme_bw()
p
ggsave(p, filename="../../Figures/Estimate_Disp/Predicted_birds.png", width=8, height=6)

new_df_birds

saveRDS(new_df_birds, "../../Objects/estimate_disp_dist/estimate_disp_dist_bird.rda")
write.csv(new_df_birds, "../../Objects/estimate_disp_dist/estimate_disp_dist_bird.csv", row.names=F)


model_df_birds$estimated_disp<-predict(best_model, model_df_birds)

min_HWI<-min(model_df_birds$HWI)
max_HWI<-max(model_df_birds$HWI)
p1<-ggplot(new_df_birds)+geom_point(aes(x=HWI, y=estimated_disp, color=factor(Diet)))+
  geom_vline(xintercept = min_HWI, color="black", linetype=2)+
  geom_vline(xintercept = max_HWI, color="black", linetype=2)+
  xlab("Hand-wing index (HWI)")+
  ylab("Estimated natal dispersal distance")+
  labs(color="Diet")+
  theme_bw()
p1
min_mass<-min(model_df_birds$log_body_mass)
max_mass<-max(model_df_birds$log_body_mass)


p2<-ggplot(new_df_birds)+geom_point(aes(x=body_mass, y=estimated_disp, color=factor(Diet)))+
  geom_vline(xintercept = min_mass, color="black", linetype=2)+
  geom_vline(xintercept = max_mass, color="black", linetype=2)+
  scale_x_log10()+
  xlab("Body mass (log transferred)")+
  ylab("Estimated natal dispersal distance")+
  labs(color="Diet")+
  theme_bw()
p2

p<-ggpubr::ggarrange(p2, p1, nrow=2)
ggsave(p, filename="../../Figures/Estimate_Disp/emperical_predicted_birds.png", width=8, height=8)
ggsave(p, filename="../../Figures/Estimate_Disp/emperical_predicted_birds.pdf", width=8, height=8)


