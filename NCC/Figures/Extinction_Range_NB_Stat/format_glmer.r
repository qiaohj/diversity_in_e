library(data.table)
library(ggplot2)
library(pivottabler)
df<-data.table(read.csv("../../Figures/NB_hist_combined/p_values_with_scaled_max_disp_dist_10km.csv", stringsAsFactors = F))
df_item<-df[exposure==" no climate resilience"&da=="with dispersal"]
for (ssp in unique(df$SSP)){
  pt <- PivotTable$new()
  pt$addData(df_item[SSP==ssp])
  pt$addColumnDataGroups("SSP", addTotal=FALSE)
  pt$addColumnDataGroups("GCM", addTotal=FALSE)
  pt$addColumnDataGroups("group", addTotal=FALSE)
  #pt$addColumnDataGroups("exposure", addTotal=FALSE)
  pt$defineCalculation(caption="AIC", calculationName="aic", summariseExpression="mean(aic)", format="%.3f")
  
  pt$defineCalculation(caption="Intercept Estimate", calculationName="Intercept_Estimate", 
                       summariseExpression="mean(Intercept_Estimate)", format="%.3f")
  pt$defineCalculation(caption="Std Error", calculationName="Intercept_Std_Error", 
                       summariseExpression="mean(Intercept_Std_Error)", format="%.3f")
  pt$defineCalculation(caption="z value", calculationName="Intercept_z_value", 
                       summariseExpression="mean(Intercept_z_value)", format="%.3f")
  pt$defineCalculation(caption="Pr(Z>|z|)", calculationName="Intercept_Pr", 
                       summariseExpression="mean(Intercept_Pr)", format="%.3f")
  
  pt$defineCalculation(caption="Niche volume Estimate", calculationName="nb_volume_Estimate", 
                       summariseExpression="mean(nb_volume_Estimate)", format="%.3f")
  pt$defineCalculation(caption="Std Error", calculationName="nb_volume_Std_Error", 
                       summariseExpression="mean(nb_volume_Std_Error)", format="%.3f")
  pt$defineCalculation(caption="z value", calculationName="nb_volume_z_value", 
                       summariseExpression="mean(nb_volume_z_value)", format="%.3f")
  pt$defineCalculation(caption="Pr(Z>|z|)", calculationName="nb_volume_Pr", 
                       summariseExpression="mean(nb_volume_Pr)", format="%.3f")
  
  pt$defineCalculation(caption="Geog range size Estimate", calculationName="N_CELL_Estimate", 
                       summariseExpression="mean(N_CELL_Estimate)", format="%.3f")
  pt$defineCalculation(caption="Std Error", calculationName="N_CELL_Std_Error", 
                       summariseExpression="mean(N_CELL_Std_Error)", format="%.3f")
  pt$defineCalculation(caption="z value", calculationName="N_CELL_z_value", 
                       summariseExpression="mean(N_CELL_z_value)", format="%.3f")
  pt$defineCalculation(caption="Pr(Z>|z|)", calculationName="N_CELL_Pr", 
                       summariseExpression="mean(N_CELL_Pr)", format="%.3f")
  
  pt$defineCalculation(caption="Mean dispersal distance Estimate", calculationName="estimated_disp_Estimate", 
                       summariseExpression="mean(estimated_disp_Estimate)", format="%.3f")
  pt$defineCalculation(caption="Std Error", calculationName="estimated_disp_Std_Error", 
                       summariseExpression="mean(estimated_disp_Std_Error)", format="%.3f")
  pt$defineCalculation(caption="z value", calculationName="estimated_disp_z_value", 
                       summariseExpression="mean(estimated_disp_z_value)", format="%.3f")
  pt$defineCalculation(caption="Pr(Z>|z|)", calculationName="estimated_disp_Pr", 
                       summariseExpression="mean(estimated_disp_Pr)", format="%.3f")
  
  pt$defineCalculation(caption="VoCC Estimate", calculationName="vocc_bio1_Estimate", 
                       summariseExpression="mean(vocc_bio1_Estimate)", format="%.3f")
  pt$defineCalculation(caption="Std Error", calculationName="vocc_bio1_Std_Error", 
                       summariseExpression="mean(vocc_bio1_Std_Error)", format="%.3f")
  pt$defineCalculation(caption="z value", calculationName="vocc_bio1_z_value", 
                       summariseExpression="mean(vocc_bio1_z_value)", format="%.3f")
  pt$defineCalculation(caption="Pr(Z>|z|)", calculationName="vocc_bio1_Pr", 
                       summariseExpression="mean(vocc_bio1_Pr)", format="%.3f")
  
  pt$addRowCalculationGroups()
  pt$renderPivot()
  pt$saveHtml(sprintf("../../Figures/NB_hist_combined/glmer_pivot/glmer_pivot_%s.html", ssp))
  
}
