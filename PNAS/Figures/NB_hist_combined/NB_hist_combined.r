library(ggplot2)
library(purrr)
setwd("/media/huijieqiao/Speciation_Extin/Sp_Richness_GCM/Script/diversity_in_e")
nb_prec<-readRDS("../../Figures/NB_hist_combined/nb_prec.rda")
nb_prec<-nb_prec+theme(strip.background.x = element_blank(),
                       strip.text.x = element_blank())

nb_temp<-readRDS("../../Figures/NB_hist_combined/nb_temp.rda")
nb_temp<-nb_temp+theme(strip.background.x = element_blank(),
              strip.text.x = element_blank())
nb_range_size<-readRDS("../../Figures/NB_hist_combined/nb_range_size.rda")
disp_dist<-readRDS("../../Figures/NB_hist_combined/disp_dist.rda")
disp_dist<-disp_dist+theme(strip.background.x = element_blank(),
                       strip.text.x = element_blank())

nb_prec_hist<-readRDS("../../Figures/NB_hist_combined/nb_prec_hist.rda")
nb_temp_hist<-readRDS("../../Figures/NB_hist_combined/nb_temp_hist.rda")
nb_range_size_hist<-readRDS("../../Figures/NB_hist_combined/nb_range_size_hist.rda")
disp_dist_hist<-readRDS("../../Figures/NB_hist_combined/disp_dist_hist.rda")


nb_prec_hist<-nb_prec_hist+theme(axis.title.x = element_blank())
legends<-get_legend(nb_prec_hist)
nb_prec_hist<-nb_prec_hist+theme(legend.position = "none")
disp_dist_hist<-disp_dist_hist+theme(legend.position = "none")
disp_dist_hist<-disp_dist_hist+scale_x_log10()
nb_temp_hist<-nb_temp_hist+theme(axis.title.x = element_blank(),
                                 axis.title.y = element_blank(),
                                 legend.position = "none")
nb_prec_hist<-nb_prec_hist+theme(axis.title.x = element_blank(), 
                                 axis.title.y = element_blank(),
                                 legend.position = "none")

nb_range_size_hist<-nb_range_size_hist+theme(axis.title.x = element_blank(), legend.position = "none")
disp_dist_hist<-disp_dist_hist+theme(axis.title.x = element_blank(), 
                                     axis.title.y = element_blank(),
                                     legend.position = "none")

list_p<-list(nb_range_size, nb_range_size_hist, nb_temp, nb_temp_hist, 
             nb_prec, nb_prec_hist, disp_dist, disp_dist_hist)
pp<-ggarrange(plotlist=list_p, nrow=4, ncol=2, widths = c(2.5,1),
          common.legend = T, legend.grob = legends, legend="right")
ggsave(pp, filename="../../Figures/NB_hist_combined/NB_hist_combined.png", width=15, height=12)
ggsave(pp, filename="../../Figures/NB_hist_combined/NB_hist_combined.pdf", width=15, height=12)
