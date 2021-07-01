library(ggplot2)
library(raster)
library(rgdal)
library(rgeos)
library(MASS)
library(cluster)
library(dplyr)

setwd("Y:/Script/diversity_in_e")

source("colors.R")
source("genCircle.R")

addEllipse <- function(mu, sigma, m = NULL, n = 100, p.interval = NULL , 
                       ci.mean = FALSE, small.sample = FALSE, 
                       do.plot = TRUE, ...){
  
  # ----------------------------------------------------------------------------
  # Some error checking
  if(small.sample & is.null(m)) message("A sample size number given by m is 
                                        required when small.sample is TRUE")
  
  if(ci.mean & is.null(m)) message("A sample size number given by m is 
                                        required when plotting confidence 
                                   ellipses of the mean with ci.mean is TRUE")
  
  
  # ----------------------------------------------------------------------------
  
  # mu is the location of the ellipse (its bivariate mean)
  # sigma describes the shape and size of the ellipse
  # p can set the predction interval for the ellipse.
  # n determines how many data points are used to draw 
  #   the ellipse. More points = smoother curves.
  
  # if ci.mean is F (default) then we are plotting quantiles of the sample, 
  # i.e. prediction ellipses, and so we set c <- 1 so that it has no effect
  # below. Else it divides the radius calculation below by sqrt(m) to include
  # the conversion from standard deviation to standard error of the mean.
  ifelse(ci.mean, 
         c.scale <- m,
         c.scale <- 1
  )
  
  # The small.sample toggles on and off (default) the small sample size 
  # correction to essentially plot the SEAc in place of the SEA. It can be 
  # used inconjuction with any prediction ellipse.
  ifelse(small.sample,
         q <- (m - 1) / (m - 2),
         q <- 1)
  
  
  # if p is NULL then plot a standard ellipse with r = 1
  # else generate a prediction ellipse that contains
  # approximately proportion p of data by scaling r
  # based on the chi-squared distribution.
  # p defaults to NULL.
  ifelse(is.null(p.interval), 
         r <- 1, 
         r <- sqrt(stats::qchisq(p.interval, df=2))
  )
  
  
  # get the eigenvalues and eigenvectors of sigma
  # if ci.mean = T then the covariance matrix is divided by the sample size
  # so as to produce confidence ellipses for the mean. Else it has no 
  # effect with c.scale = 1.
  e = eigen(sigma / c.scale)
  
  # 
  SigSqrt = e$vectors %*% diag(sqrt(e$values * q)) %*% t(e$vectors)
  
  # create a unit radius circle to transform
  cc <- genCircle(n, r)
  
  # transform the unit circle according to the covariance 
  # matrix sigma
  
  # a function to transform the points
  back.trans <- function(x) {
    return(SigSqrt %*% x + mu)
  }
  
  # apply the transformation to calculate the 
  # Maximum Likelihood estimate of the ellipse.
  
  ML.ellipse = t(apply(cc,1, back.trans))
  
  #if(grDevices::dev.cur() > 1 & do.plot) {graphics::lines(ML.ellipse, ...)}
  
  # optional return of x and y coordinates of the plotted ellipse
  return(ML.ellipse)
  
}

NDquntil <- function(nD, level) {
  n <- floor(nD * level)
  if (n > nD) 
    n <- nD
  return(n)
}
in_Ellipsoid <- stats::qchisq(0.95, 2)

args = commandArgs(trailingOnly=TRUE)
group<-args[1]

if (is.na(group)){
  group<-"Amphibians"
}

GCMs<-c("EC-Earth3-Veg", "MRI-ESM2-0", "UKESM1")
SSPs<-c("SSP119", "SSP245", "SSP585")
#VARs<-c("pr", "tasmax", "tasmin")
VARs<-c("pr", "tasmax")
start_range<-c(2000:2014)
start_layer_df<-expand.grid(GCM=GCMs, SSP=SSPs[1], VAR=VARs, Y=start_range)

var_tamplate<-"../../Raster/ENV/Annually/%s_%s_%s_%d_%s_eck4.tif"

mask<-raster("../../Raster/mask_index.tif")

start_env_layers<-readRDS("../../Objects/stacked_layers_2000_2014_df.rda")
df_list<-readRDS(sprintf("../../Objects/IUCN_List/%s.rda", group))


var_pair<-start_layer_df[which(start_layer_df$VAR=="pr"),]
var_pair1<-left_join(var_pair, start_layer_df[which(start_layer_df$VAR=="tasmax"),], by=c("GCM", "SSP", "Y"))
#var_pair2<-left_join(var_pair, start_layer_df[which(start_layer_df$VAR=="tasmin"),], by=c("GCM", "SSP", "Y"))
var_pair1<-var_pair1%>%dplyr::select(GCM, VAR.x, Y, VAR.y)
colnames(var_pair1)<-c("GCM", "PR", "Y", "TEMP")
#var_pair2<-var_pair2%>%dplyr::select(GCM, VAR.x, Y, VAR.y)
#colnames(var_pair2)<-c("GCM", "PR", "Y", "TEMP")
#var_pair<-bind_rows(var_pair1, var_pair2)
var_pair<-var_pair1
var_pair$PR_NAME<-paste(gsub("-", ".", var_pair$GCM), var_pair$PR, var_pair$Y, sep="_")
var_pair$TEMP_NAME<-paste(gsub("-", ".", var_pair$GCM), var_pair$TEMP, var_pair$Y, sep="_")
df_list<-df_list[sample(nrow(df_list), nrow(df_list)),]
i=1

if (F){
  for (i in c(1:nrow(df_list))){
    item<-df_list[i,]
    item$sp<-gsub(" ", "_", item$sp)
    if (item$area<=0){
      next()
    }
    
    
    print(paste(i, nrow(df_list), item$sp))
    occ<-readRDS(sprintf("../../Objects/IUCN_Distribution/%s/%s.rda", group, item$sp))
    if (nrow(occ)>1000){
      asdf
    }else{
      next()
    }
  }
}
example_sp<-"Scinax_x-signatus"
target_folder<-sprintf("../../Objects/Niche_Models_Mean_GCM/%s/%s", group, example_sp)
fit<-readRDS(sprintf("%s/fit.rda", target_folder))
all_v<-readRDS(sprintf("%s/occ_with_env.rda", target_folder))

mve_lines<-as.data.frame(addEllipse(fit$center, fit$cov, p.interval=0.95))
range_PR<-range(mve_lines$V1)
range_TEMP<-range(mve_lines$V2)
start_env_layers_sample<-start_env_layers[sample(nrow(start_env_layers), 5000),]
all_vs_sample<-all_v[sample(nrow(all_v), 500),]

p<-ggplot()+
  geom_point(data=start_env_layers_sample, aes(x=TEMP_MAX, y=PR), color=colors_black[4], alpha=0.4)+
  geom_point(data=all_vs_sample, aes(x=TEMP, y=PR, color=factor(in_out)))+
  scale_color_manual(breaks = c(1, 0), 
                    values=color_two)+
  geom_path(data=mve_lines,aes(x=V2, y=V1), color=colors_red[9], size=1.5)+
  geom_hline(yintercept=range_PR, linetype="dashed", color = colors_black[7], size=1.5)+
  geom_vline(xintercept=range_TEMP, linetype="dashed", color = colors_black[7], size=1.5)+
  xlim(c(30, 32.2))+
  ylim(c(0, 5000))+
  theme_bw()+
  theme(legend.position = "none",
        axis.title=element_text(size=20))+
  xlab("Annual Maximum Temperature")+ylab("Annual Precipitation")
p
ggsave(p, filename="../../Figures/Methods/niche.png", width=10, height=8)
  

