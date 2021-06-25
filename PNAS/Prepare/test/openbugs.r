library(R2OpenBUGS)

model.file <- "Y:/Sp_Richness_GCM/Script/diversity_in_e/PNAS/Prepare/disp.model.example.txt"
file.show(model.file)

df<-read.table("../../Data/Dispersal_distance/Bird_OPENBUGS/openbugs.csv", head=T, sep=",", stringsAsFactors = T)
df$meddisp<-as.numeric(df$meddisp)
df$mass<-as.numeric(df$mass)
df$wingspan<-as.numeric(df$wingspan)

meddisp<-df$meddisp
sex<-as.numeric(df$sex)
guild<-as.numeric(df$guild)
mass<-df$mass
shape<-(df$mass)^3/df$wingspan
df$shape<-shape
sp<-as.numeric(df$sp)


data <- list ("meddisp", "sex", "guild", "mass",	"shape",	"sp")

inits <- function(){
  list(prec=0.5, a=10, bs=c(1,1,NA), sd_sp=4, bg = c(NA,1,1,1), bw=c(5,2,5,2), bm = c(5,2,5,2))
}

parameters <- c("prec", "bs",	"sd_sp", "b", "m")
sim <- bugs(data, inits, parameters, model.file,
                    n.chains=3, n.iter=20000, debug=F)

v1<-log(df$meddisp)
v2<-sim$mean[["m"]]
plot(v1, v2)
print(sim)
plot(sim)
sim$
## End(Not run)

str(sim)
library(randomForest)
rf_model<-randomForest(meddisp~sex+shape+guild+wingspan, data=df)
rf_model<-randomForest(meddisp~sex+shape+guild+Order+wingspan, data=df)

pred<-predict(rf_model, df)
cor(df$meddisp, pred)

