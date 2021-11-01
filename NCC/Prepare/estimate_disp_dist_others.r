

# train the model
rf_order <- train(meddisp~HWI+log_body_mass+Migration_3+Diet+Order.y, data=df_with_family, 
                  trControl=train_control, method="rf")
# summarize results
print(rf_order)
evaluated_metrics<-rf_order$results
evaluated_metrics<-evaluated_metrics[which(evaluated_metrics$RMSE==min(evaluated_metrics$RMSE)),]
evaluated_metrics$model<-"RF"
evaluated_metrics$type<-"with order"
evaluated_metrics$clade<-"Bird"
evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)


# train the model
rf_order_family <- train(meddisp~HWI+log_body_mass+Migration_3+Diet+
                           Family.name, data=df_with_family, trControl=train_control, method="rf")
# summarize results
print(rf_order_family)
evaluated_metrics<-rf_order_family$results
evaluated_metrics<-evaluated_metrics[which(evaluated_metrics$RMSE==min(evaluated_metrics$RMSE)),]
evaluated_metrics$model<-"RF"
evaluated_metrics$type<-"with order and family"
evaluated_metrics$clade<-"Bird"
evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)
evaluated_metrics_all<-evaluated_metrics_all[,-1]

# train the model
glm_order <- train(meddisp~HWI+log_body_mass+Migration_3+Diet+Order.y, data=df_with_family, 
                   trControl=train_control, method="glm")
# summarize results
print(glm_order)
evaluated_metrics<-glm_order$results
evaluated_metrics$model<-"GLM"
evaluated_metrics$type<-"with order"
evaluated_metrics$clade<-"Bird"
evaluated_metrics<-evaluated_metrics[, -1]
evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)

df_with_family$Order.y
df_with_family$Family.name
# train the model
glm_order_family <- train(meddisp~HWI+log_body_mass+Migration_3+Diet+Order.y+Family.name, 
                          data=df_with_family, trControl=train_control, method="glm")
# summarize results
print(glm_order_family)
evaluated_metrics<-glm_order_family$results
evaluated_metrics$model<-"GLM"
evaluated_metrics$type<-"with order and family"
evaluated_metrics$clade<-"Bird"
evaluated_metrics<-evaluated_metrics[, -1]
evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)
# train the model
svm_order <- train(meddisp~HWI+log_body_mass+Migration_3+Diet+Order.y, data=df_with_family, 
                   trControl=train_control, method="svmLinear")
# summarize results
print(svm_order)
evaluated_metrics<-svm_order$results
evaluated_metrics$model<-"SVM"
evaluated_metrics$type<-"with order"
evaluated_metrics$clade<-"Bird"
evaluated_metrics<-evaluated_metrics[, -1]
evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)

# train the model
svm_order_family <- train(meddisp~HWI+log_body_mass+Migration_3+Diet+Order.y+Family.name, 
                          data=df_with_family, trControl=train_control, method="svmLinear")
# summarize results
print(svm_order_family)
evaluated_metrics<-svm_order_family$results
evaluated_metrics$model<-"SVM"
evaluated_metrics$type<-"with order and family"
evaluated_metrics$clade<-"Bird"
evaluated_metrics<-evaluated_metrics[, -1]
evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)
# train the model
brnn_order <- train(meddisp~HWI+log_body_mass+Migration_3+Diet+Order.y, data=df_with_family, 
                    trControl=train_control, method="brnn")
# summarize results
print(brnn_order)
evaluated_metrics<-brnn_order$results
evaluated_metrics<-evaluated_metrics[which(evaluated_metrics$RMSE==min(evaluated_metrics$RMSE)),]
evaluated_metrics$model<-"BRNN"
evaluated_metrics$type<-"with order"
evaluated_metrics$clade<-"Bird"
evaluated_metrics<-evaluated_metrics[, -1]
evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)

# train the model
brnn_order_family <- train(meddisp~HWI+log_body_mass+Migration_3+Diet+Family.name, 
                           data=df_with_family, trControl=train_control, method="brnn")
# summarize results
print(brnn_order_family)
evaluated_metrics<-brnn_order_family$results
evaluated_metrics<-evaluated_metrics[which(evaluated_metrics$RMSE==min(evaluated_metrics$RMSE)),]
evaluated_metrics$model<-"BRNN"
evaluated_metrics$type<-"with order and family"
evaluated_metrics$clade<-"Bird"
evaluated_metrics<-evaluated_metrics[, -1]
evaluated_metrics_all<-rbind(evaluated_metrics_all, evaluated_metrics)


saveRDS(rf_order, "../../Objects/estimate_disp_dist/models/rf_order.rda")
saveRDS(rf_order_family, "../../Objects/estimate_disp_dist/models/rf_order_family.rda")

saveRDS(glm_order, "../../Objects/estimate_disp_dist/models/glm_order.rda")
saveRDS(glm_order_family, "../../Objects/estimate_disp_dist/models/glm_order_family.rda")

saveRDS(svm_order, "../../Objects/estimate_disp_dist/models/svm_order.rda")
saveRDS(svm_order_family, "../../Objects/estimate_disp_dist/models/svm_order_family.rda")

saveRDS(brnn_order, "../../Objects/estimate_disp_dist/models/brnn_order.rda")
saveRDS(brnn_order_family, "../../Objects/estimate_disp_dist/models/brnn_order_family.rda")
