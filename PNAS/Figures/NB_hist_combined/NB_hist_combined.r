library(ggplot2)
library(purrr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
nb_prec<-readRDS("../../Figures_Full_species/NB_hist_combined/nb_prec.rda")
nb_prec<-nb_prec+theme(strip.background.x = element_blank(),
                       strip.text.x = element_blank())

nb_temp<-readRDS("../../Figures_Full_species/NB_hist_combined/nb_temp.rda")
nb_temp<-nb_temp+theme(strip.background.x = element_blank(),
              strip.text.x = element_blank())
nb_range_size<-readRDS("../../Figures_Full_species/NB_hist_combined/nb_range_size.rda")

nb_prec_hist<-readRDS("../../Figures_Full_species/NB_hist_combined/nb_prec_hist.rda")
nb_temp_hist<-readRDS("../../Figures_Full_species/NB_hist_combined/nb_temp_hist.rda")
nb_range_size_hist<-readRDS("../../Figures_Full_species/NB_hist_combined/nb_range_size_hist.rda")

nb_prec_hist<-nb_prec_hist+theme(axis.title.x = element_blank())
legends<-get_legend(nb_prec_hist)
nb_prec_hist<-nb_prec_hist+theme(legend.position = "none")
nb_temp_hist<-nb_temp_hist+theme(axis.title.x = element_blank(), legend.position = "none")
nb_range_size_hist<-nb_range_size_hist+theme(axis.title.x = element_blank(), legend.position = "none")
list_p<-list(nb_range_size, nb_range_size_hist, nb_temp, nb_temp_hist, nb_prec, nb_prec_hist)
pp<-ggarrange(plotlist=list_p, nrow=3, ncol=2, widths = c(2.5,1),
          common.legend = T, legend.grob = legends, legend="right")
ggsave(pp, filename="../../Figures_Full_species/NB_hist_combined/NB_hist_combined.png", width=15, height=8)
ggsave(pp, filename="../../Figures_Full_species/NB_hist_combined/NB_hist_combined.pdf", width=15, height=8)
