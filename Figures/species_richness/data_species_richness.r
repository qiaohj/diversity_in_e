g<-"Amphibians"

df<-readRDS("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Objects/Diversity/Amphibians/EC-Earth3-Veg_SSP119_1_1/indices_df.rda")

df[["2014"]][["species.richness"]]%>%dplyr::filter(index==4611)