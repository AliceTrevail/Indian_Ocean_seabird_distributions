#### Measure of residency ###
if (!require("devtools")) 
  install.packages("devtools")

devtools::install_github("sjmgarnier/viridis")
library(viridis)
library(magick)

# To understand whether movement behaviour explains patterns of species richness

### packages ####
library(rgdal)
library(tidyverse)
library(plyr)
library(ggalt)
library(scales)
library(viridisLite)
library(sf)

### Read in tracking data for RSFs, and remove available points ####
all <- read_csv("/Data/df_new_migration_dates_MCPcrop_restandardised.csv")%>%
  filter(use == "1") 

### Read in IO hex grid ####

IOhex <- readOGR("/Data/Polygons/IO hex grid.shp", layer = "IO hex grid")


#### Read in relevant map layers ####

world_map <- map_data("world")
colonies <- read.csv("/Data/Colony codes.csv", as.is = T)
mymap <- fortify(IOhex)

mymap_sf <- sf::st_as_sf(mymap, coords = c("long", "lat"))
st_crs(mymap_sf) <- 4326

mymap_sf_poly = st_sf(
  aggregate(
    mymap_sf$geometry,
    list(mymap_sf$id),
    function(g){st_cast(st_combine(g),"POLYGON")})) %>%
  mutate(hex_id = as.numeric(Group.1))

### Match tracking data to hex polygon ####

all.sp <- all
coordinates(all.sp) <- c("lon", "lat")
proj4string(all.sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

all$polygon <- (over(all.sp, IOhex))$ID
head(all)


### Days per hex ###

days.per.hex.ind <- plyr::ddply(all,
                                c("Species", "Colony", "Population", "ID_bird", "polygon", "year"),
                                summarise,
                                no.days = length(ID_burst)/2)
head(days.per.hex.ind)

mean.days.per.hex.sp <- plyr::ddply(days.per.hex.ind,
                                    c("Species"),
                                    summarise,
                                    mean.days = mean(no.days),
                                    se.days = sd(no.days)/sqrt(length(no.days)))
head(mean.days.per.hex.sp)


mean.days.per.hex.sp$Species <- fct_relevel(mean.days.per.hex.sp$Species, "WTSH", "TRSH", "BAPE", "TRPE", "RTTR", "WTTR", "SOTE", "BRNO", "LENO")
levels(mean.days.per.hex.sp$Species) <-  c("Wedge-tailed shearwater", 
                                                "Tropical shearwater", 
                                                "Barau's petrel", 
                                                "Trindade petrel", 
                                                "Red-tailed tropicbird", 
                                                "White-tailed tropicbird", 
                                                "Sooty tern", 
                                                "Brown noddy", 
                                                "Lesser noddy")

p.cumulative.box <- ggplot(mean.days.per.hex.sp, aes(x = Species, y = mean.days))+
  geom_point()+
  geom_errorbar(aes(ymin = mean.days - se.days, ymax = mean.days + se.days))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(y = "Mean total days per individual per hexagon, per calendar year")
  
ggsave("/Figures/Manuscript/Supp mat/Residency_cumulative days.png", 
       p.cumulative.box, width = 5, height = 5.5, units = c("in"))


ggplot(days.per.hex.ind, aes(x = Species, y = no.days))+
  geom_violin()+
  theme_bw()+
  scale_y_continuous(n.breaks = 10)+
  labs(y = "Total days per individual, per hexagon")



days.per.hex.sp <- plyr::ddply(days.per.hex.ind,
                               c("Species", "polygon"),
                               summarise,
                               mean.days = mean(no.days),
                               se.days = sd(no.days)/sqrt(length(no.days)),
                               total.days = sum(no.days))  %>% 
  mutate(hex_id = as.numeric(polygon))


days.per.hex.sp$Species <- fct_relevel(days.per.hex.sp$Species, "WTSH", "TRSH", "BAPE", "TRPE", "RTTR", "WTTR", "SOTE", "BRNO", "LENO")
levels(days.per.hex.sp$Species) <-  c("Wedge-tailed shearwater", 
                                           "Tropical shearwater", 
                                           "Barau's petrel", 
                                           "Trindade petrel", 
                                           "Red-tailed tropicbird", 
                                           "White-tailed tropicbird", 
                                           "Sooty tern", 
                                           "Brown noddy", 
                                           "Lesser noddy")

head(days.per.hex.sp)


cumulative_days_hex <- right_join(mymap_sf_poly, days.per.hex.sp, by = "hex_id", multiple = "all" )


plot_base <- list(geom_polygon(data = world_map, aes(x = long, y = lat, group = group), inherit.aes = FALSE, fill="gray60", colour = "gray60"),
                  coord_sf(xlim = c(25, 125), ylim = c(-45, 30), expand = F),
                  theme_bw(),
                  theme(panel.background = element_rect(fill = "white", colour = "grey60"), panel.grid = element_blank(),
                        strip.background = element_blank(), axis.text=element_text(size=7), axis.title=element_text(size=9),
                        plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")),
                  labs(x = "Longitude", y = "Latitude"))



p.cumulative.map <- ggplot() +
  geom_sf(data = cumulative_days_hex, aes(fill = mean.days), col = NA)+
  scale_fill_viridis_c(name = "Mean total\ndays per\nindividual", trans = "log", breaks_log(n = 5), option="magma", begin = 0.1)+
  facet_wrap(.~Species)+
  plot_base+
  geom_text(data = mean.days.per.hex.sp, mapping = aes(x = -Inf, y = -Inf, label = paste(formatC(mean.days, digits = 2, format = "f"), "±", formatC(se.days, digits = 2, format = "f")), fontface=2), hjust = -0.1, vjust = -0.5, cex = 3)
p.cumulative.map

ggsave("/Figures/Manuscript/Residency_cumulative_map.png", 
       p.cumulative.map, width = 8, height = 6, units = c("in"))



### mean, standard error & range to report in paper ####
mean((mean.days.per.hex.sp$mean.days))
sd(mean.days.per.hex.sp$mean.days)/sqrt(length(mean.days.per.hex.sp$mean.days))
min(days.per.hex.ind$no.days)
max(days.per.hex.ind$no.days)


