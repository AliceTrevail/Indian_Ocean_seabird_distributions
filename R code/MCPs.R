# Date change MCPs #
library(sp)
library(rgeos)
library(adehabitatHR)
library(tidyverse)
library(magick)

### Read in dataset with corrected migration dates ####

df_migrationdates <- read_csv("/Data/df_migration_dates.csv")

# filter to used points only
df_migrationdates_used <- df_migrationdates %>%
  filter(use == 1)

samplesizes <- df_migrationdates %>%
  group_by(Species) %>%
  summarise(n = length(unique(ID_bird)))

# limit tracking data to minimal columns for mcp
sp_migrationdates <- dplyr::select(df_migrationdates_used, lon, lat, Population)
coordinates(sp_migrationdates) <- c("lon", "lat")
proj4string(sp_migrationdates) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# project to laea so units are in m
sp_migrationdates_laea <- spTransform(sp_migrationdates, CRS("+proj=laea +lat_0=-7 +lon_0=72 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

# create mcp for each species
mp <- mcp(sp_migrationdates_laea, percent=100)
plot(sp_migrationdates_laea)
plot(mp, add=T, col = "red")


# expand mcps by buffer of mean uncertainty around locations
#error.x <- mean(all$x.se)
#error.y <- mean(all$y.se)  
#error <- mean(c(error.x, error.y))
error <- 90061.28/1000
error*1000
# error *1000 = 90061.28

mpb <- gBuffer(mp, byid = T, width = error*1000)


plot(sp_migrationdates_laea)
plot(mp, add=T, col = NA, border = "red")
plot(mpb, add=T, col = NA, border = "blue")

head(mpb)

mpb.latlong <- spTransform(mpb, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

df_MCPs <- fortify(mpb.latlong)

