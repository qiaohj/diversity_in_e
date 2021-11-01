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

df_item_2<-Sutherland_se[, c("sp", "median_dis", "max_dis")]
df_item_2[which(is.infinite(df_item_2$median_dis)), "median_dis"]<-NA
df_item_2[which(is.infinite(df_item_2$max_dis)), "max_dis"]<-NA

df_item_3<-read.table("../../Data/Dispersal_distance/Bird_Others/others.csv", head=T,
                      sep=",", stringsAsFactors = F)

df_item<-rbind(df_item_2, df_item_3)
df_item<-df_item[which(!is.na(df_item$max_dis)),]
df_item[which(df_item$sp=="Dryobates borealis"), "sp"]<-"Leuconotopicus borealis"
df_item[which(df_item$sp=="Dryocopus pileatus"), "sp"]<-"Dryocopus pileatus"
df_item[which(df_item$sp=="Phalacrocorax pelagicus"), "sp"]<-"Phalacrocorax pelagicus"



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
model_df<-df_with_family[, ..cols]
dim(model_df)
table(model_df$Diet)
f<-formulas[7]

########For random forest##########################
# define training control for 10-fold Cross Validation
evaluated_metrics_all<-NULL
predicted_all<-NULL

setSeeds <- function(method = "cv", numbers = 1, repeats = 1, tunes = NULL, seed = 1237) {
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
set.seed(100)

seeds <- setSeeds("repeatedcv", numbers=10, repeats=10, tunes=5)

train_control <- trainControl(method="repeatedcv", number=10, repeats=10, 
                              seeds=seeds)

train_control <- trainControl(method="boot", number=10, repeats=10, 
                              seeds=seeds)
tunegrid <- expand.grid(.mtry=c(1:5))

rf_no_rank <- train(as.formula(f), data=model_df, 
                    trControl=train_control, method="rf",
                    tuneGrid=tunegrid)

glm_no_rank <- train(as.formula(f), data=model_df, 
                    trControl=train_control, method="glm")

rf_no_rank$results
glm_no_rank$results

rf_no_rank$finalModel$test
pred<-rf_no_rank$finalModel$predicted
mean(rf_no_rank$finalModel$rsq)
cor(pred, rf_no_rank$finalModel$y)
caret::R2(pred, rf_no_rank$finalModel$y)
caret::RMSE(pred, rf_no_rank$finalModel$y)
