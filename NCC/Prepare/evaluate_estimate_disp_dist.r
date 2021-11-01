
library(ggplot2)
library(ggpubr)
library(caret)
library(data.table)
library(randomForest)
df_with_family_table<-read.csv("../../Data/Dispersal_distance/bird.csv")
evaluated_metrics_birds<-readRDS("../../Objects/estimate_disp_dist/models/evaluated_metrics_birds.rda")
best_model_info<-evaluated_metrics_birds[which(evaluated_metrics_birds$RMSE==min(evaluated_metrics_birds$RMSE)),]
best_model<-readRDS(sprintf("../../Objects/estimate_disp_dist/models/%s_no_rank_birds.rda", tolower(best_model_info$model)))
#best_model<-readRDS(sprintf("../../Objects/estimate_disp_dist/models/%s_no_rank.rda", "rf"))

best_model<-best_model[[best_model_info$formulas]]
df_with_family_table<-data.table(df_with_family_table)

table(df_with_family_table$sp)
plot(df_with_family_table$max_dis, df_with_family_table$estimated_disp_dist)

n_order_bird<-data.frame(table(df_with_family_table$Order.x))
write.csv(n_order_bird, "../../Objects/estimate_disp_dist/n_order_bird.csv")

n_family_bird<-data.frame(table(df_with_family_table$Family.name))
write.csv(n_family_bird, "../../Objects/estimate_disp_dist/n_family_bird.csv")
i=7
n_family_bird$cor_train<--1
n_family_bird$cor_test<--1
for (i in c(1:nrow(n_family_bird))){
  print(i)
  item<-n_family_bird[i,]
  if (item$Freq<5){
    next()
  }
  tunegrid <- expand.grid(.mtry=c(1:5))
  train<-df_with_family_table[Family.name!=item$Var1]
  test<-df_with_family_table[Family.name==item$Var1]
  train_control <- trainControl(method="repeatedcv", number=10, repeats=10)
  
  model <- train(as.formula(best_model_info$formulas), data=train, 
                      trControl=train_control, method="rf", 
                      tuneGrid=tunegrid, ntree=1000)
  predict_train<-predict(model, train)
  cor_train<-cor(train$max_dis, predict_train)
  predict_test<-predict(model, test)
  cor_test<-cor(test$max_dis, predict_test)
  n_family_bird[i, "cor_train"]<-cor_train
  n_family_bird[i, "cor_test"]<-cor_test
}
rm("cor_test")
n_family_bird_test<-n_family_bird[which(n_family_bird$cor_test>0),]
write.csv(n_family_bird_test, "../../Objects/estimate_disp_dist/n_family_bird_test.csv")



#For mammals
df_with_family_table<-readRDS("../../Data/Dispersal_distance/mammal.rda")
evaluated_metrics_mammals<-readRDS("../../Objects/estimate_disp_dist/models/evaluated_metrics_mammals.rda")
best_model_info<-evaluated_metrics_mammals[which(evaluated_metrics_mammals$RMSE==min(evaluated_metrics_mammals$RMSE)),]
best_model<-readRDS(sprintf("../../Objects/estimate_disp_dist/models/%s_no_rank_mammals.rda", 
                            tolower(best_model_info$model)))
#best_model<-readRDS(sprintf("../../Objects/estimate_disp_dist/models/%s_no_rank.rda", "rf"))

best_model<-best_model[[best_model_info$formulas]]
df_with_family_table<-data.table(df_with_family_table)

n_family_mammals<-data.frame(table(df_with_family_table$Family))
write.csv(n_family_mammals, "../../Objects/estimate_disp_dist/n_family_mammals.csv")
i=7
n_family_mammals$cor_train<--1
n_family_mammals$cor_test<--1
for (i in c(1:nrow(n_family_mammals))){
  print(i)
  item<-n_family_mammals[i,]
  if (item$Freq<5){
    next()
  }
  tunegrid <- expand.grid(.mtry=c(1:5))
  train<-df_with_family_table[Family!=item$Var1]
  test<-df_with_family_table[Family==item$Var1]
  train_control <- trainControl(method="repeatedcv", number=10, repeats=10)
    model <- train(as.formula(best_model_info$formulas), data=train, 
                   trControl=train_control, method="rf", 
                   tuneGrid=tunegrid, ntree=1000)
  predict_train<-predict(model, train)
  cor_train<-cor(train$max_dis, predict_train)
  predict_test<-predict(model, test)
  cor_test<-cor(test$max_dis, predict_test)
  n_family_mammals[i, "cor_train"]<-cor_train
  n_family_mammals[i, "cor_test"]<-cor_test
}
rm("cor_test")
n_family_mammals_test<-n_family_mammals[which(n_family_mammals$cor_test>0),]
write.csv(n_family_mammals_test, "../../Objects/estimate_disp_dist/n_family_mammals_test.csv")

