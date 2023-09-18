### Indian ocean hex grid ####

# making hex grid
library(sp)
library(rgeos)
library(rgdal)
library(geogrid)
library(raster)
library(raster)

# species distributions
library(GISTools)
library(GADMTools)
library(ggalt)
library(viridis)
library(mapproj)
library(tidyverse)
library(plyr)

world_map <- map_data("world")
colonies <- read.csv("/Data/Colony codes.csv", as.is = T)


# import IO shape file
ocean <- readOGR("/Data/Polygons/Indian Ocean.shp")
plot(ocean, col = "blue")
proj4string(ocean)

# project to equal area so that units are metres
ocean.r <- spTransform(ocean, CRS("+proj=laea +lat_0=-7 +lon_0=72 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
plot(ocean.r)

### create hexagonal grid ####

size <- 200000  # grid size (dist between centroids); units = m

# buffer around extent of shape file to entirely cover area with grid
ext <- as(extent(ocean.r) + size, "SpatialPolygons")

# create grid of points in hexagonal arrangement. Specify offset so grid is not random
hex_points <- spsample(ext, type = "hexagonal", cellsize = size,  offset = c(0, 0))
proj4string(hex_points) <- CRS("+proj=laea +lat_0=-7 +lon_0=72 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# convert to grid of hexagonal polygons
hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = size)
proj4string(hex_grid)

# re-project to lat/long
hex_l <- spTransform(hex_grid, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(ocean, col = "blue")
plot(hex_l, border = "orange", add = T)

# keep only polygons that intersect with IO shape file
hex_l2 <- hex_l[ocean,]
plot(ocean, col = "blue")
plot(hex_l2, border = "orange", add = T)

### save hexagonal grid ####
pid <- sapply(slot(hex_l2, "polygons"), function(x) slot(x, "ID"))
p.df <- data.frame( ID=1:length(hex_l2), row.names = pid)
spdf <- SpatialPolygonsDataFrame(hex_l2, p.df)


#### re-save hex grid ####

IOhex <- readOGR("/Users/at687/Documents/BIOT/Non-breeding/Data/Polygons/IO hex grid.shp", layer = "IO hex grid")
plot(IOhex, fill = "red")
setwd("/Users/at687/Documents/BIOT/Non-breeding/Data/Polygons")
raster::shapefile(IOhex, "IOhex_2.shp")


#### ************************ ####
### Observed species richness ####


##### Read in data
all <- read_csv("/Data/df_new_migration_dates_MCPcrop_restandardised.csv")%>%
  filter(use == "1")

# make spatial points df
all2 <- all
coordinates(all2) <- c("lon", "lat")
proj4string(all2) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

all$polygon <- (over(all2, hex_l2))

### number of individuals per cell
# subset data to one point per cell per species

all_subset_sp <- all[!duplicated(all[c("Species", "polygon")]),]
all_subset_sp2 <- all_subset_sp
coordinates(all_subset_sp2) <- c("lon", "lat")
proj4string(all_subset_sp2) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

poly.counts<-function (pts, polys) 
  colSums(gContains(polys, pts, byid=TRUE))

counts_sp <- poly.counts(all_subset_sp2, hex_l2)
counts_sp <- stack(counts_sp)
mymap <- fortify(hex_l2)
head(counts_sp)


#### ************************* ####
#### Species richness by month ####

all$month <- format.Date(all$date, "%b")
unique(all$month)
head(all)

months <- plyr::ddply(all, ~month, summarise, no.sp =length(unique(Species)))
months

month_hex <- as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(month_hex) <- c("values","ind","month")

for (i in unique(all$month)){
  all.m <- all[all$month == i,]
  all.m.subset <- all.m[!duplicated(all.m[c("Species", "polygon")]),]
  all.m_sp <- all.m.subset
  coordinates(all.m.subset) <- c("lon", "lat")
  proj4string(all.m.subset) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  
  counts_sp <- poly.counts(all.m.subset, hex_l2)
  counts_sp <- stack(counts_sp)
  counts_sp$month <- unique(all.m$month)
  print(head(counts_sp))
  
  month_hex <- rbind(month_hex, counts_sp)
}

month_hex$month <- factor(month_hex$month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
months$month <- factor(months$month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))


monthcols <- all[!duplicated(all[c("month", "Colony")]),c("month", "Colony")]
monthcols[c("Col_Lat", "Col_Long")] <- colonies[match(monthcols$Colony, colonies$Colony_Code), c("Col_Lat", "Col_Long")]
monthcols$month <- factor(monthcols$month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

head(monthcols)

write.csv(month_hex, "/Data/df_hex_observed_richness_bymonth.csv", row.names = F)



## convert hexagon map layers to sf objects 
mymap_sf <- sf::st_as_sf(mymap, coords = c("long", "lat"))
st_crs(mymap_sf) <- 4326

mymap_sf_poly = st_sf(
  aggregate(
    mymap_sf$geometry,
    list(mymap_sf$id),
    function(g){st_cast(st_combine(g),"POLYGON")})) %>%
  dplyr::rename(hex_id = Group.1)


