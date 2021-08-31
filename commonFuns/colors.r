#https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=9
colors_red<-c("#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a",
              "#ef3b2c", "#cb181d", "#a50f15", "#67000d")
colors_green<-c("#f7fcfd", "#e5f5f9", "#ccece6", "#99d8c9", "#66c2a4",
              "#41ae76", "#238b45", "#006d2c", "#00441b")
colors_blue<-c("#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6",
                "#4292c6", "#2171b5", "#08519c", "#08306b")
colors_purple<-c("#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8",
               "#807dba", "#6a51a3", "#54278f", "#3f007d")
colors_black<-c("#ffffff", "#f0f0f0", "#d9d9d9", "#bdbdbd", "#969696",
                 "#737373", "#525252", "#252525", "#000000")
map_background<-"#f5f5f2"
mask_color<-colors_black[3]

color_two<-c(colors_red[8], colors_blue[8])
color_two_map<-c(mask_color, colors_red[8])

color_ssp<-c("SSP119"=colors_blue[7],
             "SSP245"=colors_green[7],
             "SSP585"=colors_red[7])
color_survive<-c("extant"=colors_blue[7],
             "extinct"=colors_red[7])
color_dispersal<-c("#FE5F55",
             "#7A9E9F",
             colors_green[7])

linetype_gcm<-c("EC-Earth3-Veg"=1,
                "MRI-ESM2-0"=2,
                "UKESM1"=5)
color_causation<-c("Temperature"=colors_red[7],
                   "Precipitation"=colors_blue[7])

color_causation_ul<-c("Temperature (upper limit)"=colors_red[5],
                   "Temperature (lower limit)"=colors_red[7],
                   "Precipitation (upper limit)"=colors_blue[5],
                   "Precipitation (lower limit)"=colors_blue[7])

color_dipsersal_type<-c("neither higher elevation nor higher latitude"=colors_red[7],
                        "higher elevation and/or higher latitude"=colors_blue[7])

library(RColorBrewer)
linetype_ssp<-c("SSP119"=1,
             "SSP245"=2,
             "SSP585"=5)

color_groups<-c("Amphibians"=colors_green[7],
             "Birds"=colors_blue[7],
             "Reptiles"=colors_purple[7],
             "Mammals"=colors_red[7])

color_groups<-c("Birds"=colors_blue[7],
                "Mammals"=colors_red[7])

#display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, colorblindFriendly=T)
color_keyspots<-brewer.pal(n = 5, name = 'Dark2')

names(color_keyspots)<-c("Amazon Basin", "central North America", "Congo Basin",
                         "Qinghai-Tibet Plateau", "Southeast Asia")



color_da<-c("no dispersal"=colors_red[7],
             "with dispersal"=colors_blue[7])

map_theme<-theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_rect(fill = map_background, color = NA), 
  panel.background = element_blank(), 
  legend.background = element_rect(fill = map_background, color = NA),
  panel.border = element_blank(),
  legend.position="none"
)

linetype_label<-c(" no climate resilience, no dispersal"=1,
                  "climate resilience, no dispersal"=6,
                  " no climate resilience, with dispersal"=2,
                  "climate resilience, with dispersal"=3)
