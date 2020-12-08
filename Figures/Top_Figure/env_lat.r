library(raster)
library(dplyr)
library(ggplot2)
lm_eqn <- function(p){
  m <- lm(temp ~ alt+abs_y, p);
  intercept<-as.numeric(coef(m)[1])
  alt_a<-as.numeric(coef(m)[2])
  y_a<-as.numeric(coef(m)[3])
  eq <- substitute(italic(temp) == a + b %.% italic(alt)+ c %.% italic(y)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        c = format(unname(coef(m)[3]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  return(data.frame(intercept=intercept,
                    alt_a=alt_a,
                    y_a=y_a,
                    eq=sprintf("temp = %.3f%.3f*elevation%.3f*y, r^2 = %.3f", 
                               unname(coef(m)[1]),
                               unname(coef(m)[2]),
                               unname(coef(m)[3]),
                               summary(m)$r.squared)))
  
}


temp<-raster("../../Raster/bioclim/bio1_eck4_10km.tif")
head(rasterToPoints(temp))
alt<-raster("../../Raster/ALT/alt_eck4_high_res.tif")
p<-data.frame(rasterToPoints(alt))
p$temp<-raster::extract(temp,  p[, c("x", "y")])
colnames(p)[3]<-"alt"
p<-p%>%dplyr::filter(!is.na(p$temp))
p$abs_y<-abs(p$y)

p$abs_y<-p$abs_y/1000
p$alt<-p$alt

lm.temp <- lm(data=p, temp ~ alt+abs_y)
summary(lm.temp)
grid.lines = 26
alt.pred <- seq(min(p$alt), max(p$alt), length.out = grid.lines)
y.pred <- seq(min(p$abs_y), max(p$abs_y), length.out = grid.lines)
xy <- expand.grid( alt = alt.pred, abs_y = y.pred)
temp.pred <- matrix(predict(lm.temp, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
p$fit.temp <- predict(lm.temp)
# scatter plot with regression plane
library("plot3D")
info<-lm_eqn(p)
p_sample<-p[sample(nrow(p), 1000),]


scatter3D(p_sample$alt, p_sample$abs_y, p_sample$temp, pch = 18, cex = 0.5, 
          theta = 30, phi = 20, ticktype = "detailed",
          xlab = "Elevation (m)", ylab = "Latitude (km)", zlab = "Temperature",
          surf = list(x = alt.pred, y = y.pred, z = temp.pred,  
                      facets = NA), main = info$eq)
