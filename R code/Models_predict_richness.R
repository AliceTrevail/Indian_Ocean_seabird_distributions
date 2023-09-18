library(rgdal)
library(sp)
library(plyr)
library(pROC)
library(ROCR)
library(maps)
library(ggalt)
library(tidyverse)
library(sf)

##### Read in tracking data & models ######
all.sp <-read_csv("/Data/df_new_migration_dates_MCPcrop_restandardised.csv")

species <- c("BAPE", "BRNO", "LENO", "RTTR", "SOTE", "TRPE", "TRSH", "WTSH", "WTTR")

for (i in species){
  m <- readRDS(paste0("/Data/Models/Main/", i, "_main.rds"))
  assign(paste0("mod_", i), m, envir = .GlobalEnv)
}

rm(m)


#### Months per population #####

all.sp$Month <-  as.integer(format.Date(all.sp$date, "%m"))
head(all.sp)
pop.months <- all.sp %>%
  group_by(Population, Month) %>%
  summarise(Tracking = "TRUE") %>% # are birds tracked during this month?
  group_by(Population) %>%
  summarise(count = n())



#### Read in IO hex grid #####

IOhex <- readOGR("/Data/Polygons/IO hex grid.shp", layer = "IO hex grid")


####### Predict models with Birdlife & MLC colonies #######

## Read in hex env and coldist data 

hex.env.coldist <- read.csv("/Data/Hex env data/For models/hex_env_coldist.csv", as.is = T)

coldist_BL_MLC <- read.csv("/Data/BL_MLC_coldists.csv", as.is = T)
colnames(coldist_BL_MLC)

colnames(coldist_BL_MLC)[1] <- "polygon.id"
hex.env.coldist <- merge(hex.env.coldist, coldist_BL_MLC, by = "polygon.id")

head(hex.env.coldist)
str(hex.env.coldist)


## run predictions
species <- c("BAPE", "BRNO", "LENO", "RTTR", "SOTE", "TRPE", "TRSH", "WTSH", "WTTR")

predict.all <- as.data.frame(matrix(ncol = 8, nrow = 0))
colnames(predict.all) <- c("polygon.id", "Data_Year", "Data_Month", "Predict_Colony", "Predict_year", "p.fit", "p.fit.se", "Species")

for (sp in species){
  df <- as.data.frame(filter(all.sp, Species == sp))
  mod <- get(paste0("mod_", sp))
  
  df$Month <-  as.integer(format.Date(df$date, "%m"))
  
  ## subset environment to within species range & non-breeding months ##
  hex.predict <- subset(hex.env.coldist, hex.env.coldist[,sp] == 1 & Month %in% unique(df$Month))
  hex.predict$ColDist <- hex.predict[,paste0(sp, "_Coldist_BL_MLC")]
  head(hex.predict)
  
  ## standardise predict variables by use available df ##
  vars <- c("Bath_90km", "Slope_90km", "Chl_log",
            "SST_90km_1day_mean", "SSTanom_90km_1day_mean", "SSTgrad_log",
            "windSp_90km_1day_mean", "ssh_90km_1day_mean", "eke_log",
            "windWSC_90km_7day_mean", "ColDist")
  
  for (j in vars){
    hex.predict[,paste0(j, "_st")] <- (hex.predict[,j]-mean(df[,j], na.rm = T))/sd(df[,j], na.rm = T)
  }
  
  
  ### run for parameters at each colony & year [to marginalise] ##
  predict.sp <- as.data.frame(matrix(ncol = 8, nrow = 0))
  colnames(predict.sp) <- c("polygon.id", "Data_Year", "Data_Month", "Predict_Colony", "Predict_year", "p.fit", "p.fit.se", "Species")
  
  cols <- unique(df$Colony)
  yrs <- unique(df$year)
  
  for (col in cols){
    for (yr in yrs){
      
      hex.predict$Colony <- col
      hex.predict$year <- as.factor(yr)
      
      p1 <- mgcv::predict.gam(mod, newdata = hex.predict, type = "response", se.fit = T)
      hex.predict$p.fit <- as.data.frame(p1[1])[,"fit"]
      hex.predict$p.fit.se <- as.data.frame(p1[2])[,"se.fit"]
      
      predict.out <- hex.predict[,c("polygon.id", "Year", "Month", "Colony", "year", "p.fit", "p.fit.se")]
      colnames(predict.out) <- c("polygon.id", "Data_Year", "Data_Month", "Predict_Colony", "Predict_year", "p.fit", "p.fit.se")
      predict.out$Species <- sp 
      
      predict.sp <- rbind(predict.sp, predict.out)
      
      print(paste0(sp, ", ", col, ", ", yr))
    }
  }
  
  predict.all <- rbind(predict.all, predict.sp)
  
}


##### Marginalise results over parameter values for colonies and years ####

names(predict.all)
predict.all.c <- predict.all[!is.na(predict.all$p.fit),]
### here by calculating the mean
predict.marginalised <- plyr::ddply(predict.all.c,
                                    c("polygon.id", "Species", "Data_Year", "Data_Month"),
                                    summarise,
                                    p.fit.m = mean(p.fit, na.rm=T),
                                    p.fit.m.se = sd(p.fit, na.rm=T)/sqrt(length(polygon.id)))
head(predict.marginalised)


##### Calculate mean probability for each month ####

predict.month <- plyr::ddply(predict.marginalised,
                             c("polygon.id", "Species", "Data_Month"),
                             summarise,
                             p.fit.mo = mean(p.fit.m, na.rm=T),
                             p.fit.mo.se = sd(p.fit.m, na.rm=T)/sqrt(length(polygon.id)))

head(predict.month)


##### Sum presences for each year ####

# first calculate threshold of presence for each species


species <- c("BAPE", "BRNO", "LENO", "RTTR", "SOTE", "TRPE", "TRSH", "WTSH", "WTTR")

thresholds <- as.data.frame(matrix(ncol = 2, nrow = 0))
colnames(thresholds) <- c("Species", "Threshold")

sp <- "LENO"

for (sp in species){
  df <- as.data.frame(filter(all.sp, Species == sp))
  mod <- get(paste0("mod_", sp))

  p2 <- predict(mod, df, type="response")

  roccurve <- pROC::roc(df$use ~ as.numeric(p2))
  threshold = pROC::coords(roccurve, "best", ret = "threshold")[1,]

  t <- as.data.frame(bind_cols(sp, threshold))

  thresholds <- rbind(thresholds, t)
}

colnames(thresholds) <- c("Species", "threshold")


predict.month$threshold <- thresholds[match(predict.month$Species, thresholds$Species), "threshold"]
predict.month$p.use <- ifelse(predict.month$p.fit.mo > predict.month$threshold, 1, 0)
head(predict.month)


### Sum presences per polygon

predict.sum <- plyr::ddply(predict.month,
                           c("polygon.id", "Species"),
                           summarise,
                           months.present = sum(p.use, na.rm= T))

head(predict.sum)

# if present for one or more months = present
predict.sum$p.sp <- ifelse(predict.sum$months.present >=1, 1, 0)
head(predict.sum)


##### Species richness per cell #####

predict.richness <- plyr::ddply(predict.sum,
                                c("polygon.id"),
                                summarise,
                                sp.richness = sum(p.sp, na.rm= T))
head(predict.richness)


#### Plot predicted richness by month ####
# sum presences per polygon, per month #
predict.sum.month <- plyr::ddply(predict.month,
                                 c("polygon.id", "Data_Month"),
                                 summarise,
                                 sp.present = sum(p.use, na.rm= T))

predict.sum.month$month <- as.factor(predict.sum.month$Data_Month)
levels(predict.sum.month$month)
levels(predict.sum.month$month) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
head(predict.sum.month)


predict.richness.month<- predict.sum.month[predict.sum.month$polygon.id %in% hex.IO$x,]%>%
  mutate(hex_id = as.character(polygon.id-1))

predicted_richness_hex_month <- right_join(mymap_sf_poly, predict.richness.month, by = "hex_id", multiple = "all")



##### Area of high species richness ####
library(raster)

IOhex <- readOGR("/Data/Polygons/IO hex grid.shp", layer = "IO hex grid")

IOhex$area <- as.numeric(area(IOhex))
head(IOhex)

hex.IO <- read.csv("/Data/Hex env data/For models/hex_id_IO.csv", as.is = T)

predict.IO <- predict.richness[predict.richness$polygon.id %in% hex.IO$x,]%>%
  mutate(hex_id = as.character(polygon.id-1))

predict.IO$area <- IOhex@data[match(predict.IO$polygon.id, IOhex@data$ID), "area"]

predict.IO$area_km <- predict.IO$area/1e+6

sum(predict.IO$area_km[predict.IO$sp.richness == "9"])
NROW(predict.IO[predict.IO$sp.richness == "9",])
NROW(predict.IO[predict.IO$sp.richness == "9",])/NROW(predict.IO)*100

sum(predict.IO$area_km[predict.IO$sp.richness %in% c("7", "8", "9")])
NROW(predict.IO[predict.IO$sp.richness %in% c("7", "8", "9"),])
NROW(predict.IO[predict.IO$sp.richness %in% c("7", "8", "9"),])/NROW(predict.IO)*100


predict.richness.month.IO <- predict.richness.month[predict.richness.month$polygon.id %in% hex.IO$x,]%>%
  mutate(hex_id = as.character(polygon.id-1))

predict.richness.month.IO$area <- IOhex@data[match(predict.richness.month.IO$polygon.id, IOhex@data$ID), "area"]
predict.richness.month.IO$area_km <- predict.richness.month.IO$area/1e+6

head(predict.richness.month.IO)
month_area <- predict.richness.month.IO %>%
  filter(sp.present > 6) %>%
  group_by(month) %>%
  summarise(total_area = sum(area_km))



##### ************************************* #####
##### Predict GAM for tracked colonies only #####



####### Predict models #######

species <- c("BAPE", "BRNO", "LENO", "RTTR", "SOTE", "TRPE", "TRSH", "WTSH", "WTTR")


predict.all <- as.data.frame(matrix(ncol = 8, nrow = 0))
colnames(predict.all) <- c("polygon.id", "Data_Year", "Data_Month", "Predict_Colony", "Predict_year", "p.fit", "p.fit.se", "Species")

for (sp in species){
  df <- as.data.frame(filter(all.sp, Species == sp))
  mod <- get(paste0("mod_", sp))
  
  df$Month <-  as.integer(format.Date(df$date, "%m"))
  
  ## subset environment to non-breeding months and  tracked colony range ##
  
  sp.df <- df
  coordinates(sp.df) <- c("lon", "lat")
  proj4string(sp.df) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  df$polygon <- unlist(over(sp.df, IOhex))
  
  hex.predict <- subset(hex.env.coldist, Month %in% unique(df$Month) & polygon.id %in% unique(df$polygon))
  
  ## calculate mean colony distance
  poly.cdist <- plyr::ddply(df, c("Colony", "polygon"), summarise,
                            ColDist = mean(ColDist))
  poly.cdist2 <- plyr::ddply(poly.cdist, c("polygon"), summarise,
                            ColDist = min(ColDist))

  hex.predict$ColDist <- poly.cdist2[match(hex.predict$polygon.id, poly.cdist2$polygon), "ColDist"]
  
  ## standardise predict variables by use available df ##
  vars <- c("Bath_90km", "Slope_90km", "Chl_log",
            "SST_90km_1day_mean", "SSTanom_90km_1day_mean", "SSTgrad_log",
            "windSp_90km_1day_mean", "ssh_90km_1day_mean", "eke_log",
            "windWSC_90km_7day_mean", "ColDist")
  
  for (j in vars){
    hex.predict[,paste0(j, "_st")] <- (hex.predict[,j]-mean(df[,j], na.rm = T))/sd(df[,j], na.rm = T)
  }
  
  
  ### run for parameters at each colony & year [to marginalise] ##
  predict.sp <- as.data.frame(matrix(ncol = 8, nrow = 0))
  colnames(predict.sp) <- c("polygon.id", "Data_Year", "Data_Month", "Predict_Colony", "Predict_year", "p.fit", "p.fit.se", "Species")
  
  cols <- unique(df$Colony)
  yrs <- unique(df$year)
  
  for (col in cols){
    for (yr in yrs){
      
      hex.predict$Colony <- col
      hex.predict$year <- as.factor(yr)
      
      p1 <- mgcv::predict.gam(mod, newdata = hex.predict, type = "response", se.fit = T)
      hex.predict$p.fit <- as.data.frame(p1[1])[,"fit"]
      hex.predict$p.fit.se <- as.data.frame(p1[2])[,"se.fit"]
      
      predict.out <- hex.predict[,c("polygon.id", "Year", "Month", "Colony", "year", "p.fit", "p.fit.se")]
      colnames(predict.out) <- c("polygon.id", "Data_Year", "Data_Month", "Predict_Colony", "Predict_year", "p.fit", "p.fit.se")
      predict.out$Species <- sp 
      
      predict.sp <- rbind(predict.sp, predict.out)
      
      print(paste0(sp, ", ", col, ", ", yr))
    }
  }
  
  predict.all <- rbind(predict.all, predict.sp)
  
}


head(predict.all)


##### Marginalise results over parameter values for colonies and years ####

names(predict.all)
predict.all.c <- predict.all[!is.na(predict.all$p.fit),]
### here by calculating the mean
predict.marginalised <- plyr::ddply(predict.all.c,
                                    c("polygon.id", "Species", "Data_Year", "Data_Month"),
                                    summarise,
                                    p.fit.m = mean(p.fit, na.rm=T),
                                    p.fit.m.se = sd(p.fit, na.rm=T)/sqrt(length(polygon.id)))
head(predict.marginalised)


##### Calculate mean probability for each month ####

predict.month <- plyr::ddply(predict.marginalised,
                             c("polygon.id", "Species", "Data_Month"),
                             summarise,
                             p.fit.mo = mean(p.fit.m, na.rm=T),
                             p.fit.mo.se = sd(p.fit.m, na.rm=T)/sqrt(length(polygon.id)))



head(predict.month)


##### Sum presences for each year ####

predict.month$threshold <- thresholds[match(predict.month$Species, thresholds$Species), "threshold"]
predict.month$p.use <- ifelse(predict.month$p.fit.mo > predict.month$threshold, 1, 0)


### Sum presences per polygon

predict.sum <- plyr::ddply(predict.month,
                           c("polygon.id", "Species"),
                           summarise,
                           months.present = sum(p.use, na.rm= T))

head(predict.sum)

# if present for one or more months = present
predict.sum$p.sp <- ifelse(predict.sum$months.present >=1, 1, 0)
head(predict.sum)


##### Species richness per cell #####

predict.richness <- plyr::ddply(predict.sum,
                                c("polygon.id"),
                                summarise,
                                sp.richness = sum(p.sp, na.rm= T))
head(predict.richness)


predict.richness2 <- plyr::ddply(predict.sum[!predict.sum$Species %in% c("SOTE", "MAPE"),],
                                 c("polygon.id"),
                                 summarise,
                                 sp.richness = sum(p.sp, na.rm= T))
head(predict.richness2)


##### ************************** ####
##### Predict with MPA network ####

##### Read in tracking data & models ######
all.sp <-read_csv("/Data/df_new_migration_dates_MCPcrop_restandardised.csv")

species <- c("BAPE", "BRNO", "LENO", "RTTR", "SOTE", "TRPE", "TRSH", "WTSH", "WTTR")
 
for (i in species){
  m <- readRDS(paste0("/Data/Models/MPA/", i, "_MPA.rds"))
  assign(paste0("mod_", i), m, envir = .GlobalEnv)
}

rm(m)


#### Parameter estimates MPA #####


MPAcoeffs <- as.data.frame(matrix(ncol = 3, nrow = 9))
colnames(MPAcoeffs) <- c("Species", "MPA_inoutfoutside", "MPA_se")
MPAcoeffs$Species <- species

for (i in species){
  mod <- get(paste0("mod_", i))
  s <- summary(mod)
  
  c <- s$p.coeff["MPA_inoutfoutside"]
  se <- s$se["MPA_inoutfoutside"]
  
  MPAcoeffs$MPA_inoutfoutside[MPAcoeffs$Species == i] <- c
  MPAcoeffs$MPA_se[MPAcoeffs$Species == i] <- se
  
  
}

MPAcoeffs


## Read in hex env and coldist data 

hex.env.coldist <- read.csv("Hex env data/For models/hex_env_coldist.csv", as.is = T)

coldist_BL_MLC <- read.csv("/Data/BL_MLC_coldists.csv", as.is = T)
colnames(coldist_BL_MLC)

colnames(coldist_BL_MLC)[1] <- "polygon.id"
hex.env.coldist <- merge(hex.env.coldist, coldist_BL_MLC, by = "polygon.id")

head(hex.env.coldist)
str(hex.env.coldist)


#### Combine with Chagos & full MPA network ####

hex.ChagosMPA <- read.csv("/Data/Hex env data/For models/hex_MPA.csv", as.is = T)
hex.MPA <- read.csv("/Data/Hex env data/For models/hex_MPAs_IOallMPAs.csv", as.is = T)
head(hex.MPA)

hex.env.coldist$MPA_inoutf <- hex.MPA[match(hex.env.coldist$polygon.id, hex.MPA$polygon.id), "MPA_inoutf"]
hex.env.coldist$MPA_CHA_inoutf <- hex.ChagosMPA[match(hex.env.coldist$polygon.id, hex.ChagosMPA$polygon.id), "MPA_inoutf"]
head(hex.env.coldist)




####### Predict models w/ full MPA network #######

species <- c("BAPE", "BRNO", "LENO",  "RTTR", "SOTE", "TRPE", "TRSH", "WTSH", "WTTR")


predict.all <- as.data.frame(matrix(ncol = 8, nrow = 0))
colnames(predict.all) <- c("polygon.id", "Data_Year", "Data_Month", "Predict_Colony", "Predict_year", "p.fit", "p.fit.se", "Species")

for (sp in species){
  df <- as.data.frame(filter(all.sp, Species == sp))
  mod <- get(paste0("mod_", sp))
  
  df$Month <-  as.integer(format.Date(df$date, "%m"))
  
  ## subset environment to within species range & non-breeding months ##
  hex.predict <- subset(hex.env.coldist, hex.env.coldist[,sp] == 1 & Month %in% unique(df$Month))
  hex.predict$ColDist <- hex.predict[,paste0(sp, "_Coldist_BL_MLC")]
  head(hex.predict)
  
  ## standardise predict variables by use available df ##
  vars <- c("Bath_90km", "Slope_90km", "Chl_log",
            "SST_90km_1day_mean", "SSTanom_90km_1day_mean", "SSTgrad_log",
            "windSp_90km_1day_mean", "ssh_90km_1day_mean", "eke_log",
            "windWSC_90km_7day_mean", "ColDist")
  
  for (j in vars){
    hex.predict[,paste0(j, "_st")] <- (hex.predict[,j]-mean(df[,j], na.rm = T))/sd(df[,j], na.rm = T)
  }
  
  
  ### run for parameters at each colony & year [to marginalise] ##
  predict.sp <- as.data.frame(matrix(ncol = 8, nrow = 0))
  colnames(predict.sp) <- c("polygon.id", "Data_Year", "Data_Month", "Predict_Colony", "Predict_year", "p.fit", "p.fit.se", "Species")
  
  cols <- unique(df$Colony)
  yrs <- unique(df$year)
  
  for (col in cols){
    for (yr in yrs){
      
      hex.predict$Colony <- col
      hex.predict$year <- as.factor(yr)
      
      p1 <- mgcv::predict.gam(mod, newdata = hex.predict, type = "response", se.fit = T)
      hex.predict$p.fit <- as.data.frame(p1[1])[,"fit"]
      hex.predict$p.fit.se <- as.data.frame(p1[2])[,"se.fit"]
      
      predict.out <- hex.predict[,c("polygon.id", "Year", "Month", "Colony", "year", "p.fit", "p.fit.se")]
      colnames(predict.out) <- c("polygon.id", "Data_Year", "Data_Month", "Predict_Colony", "Predict_year", "p.fit", "p.fit.se")
      predict.out$Species <- sp 
      
      predict.sp <- rbind(predict.sp, predict.out)
      
      print(paste0(sp, ", ", col, ", ", yr))
    }
  }
  
  predict.all <- rbind(predict.all, predict.sp)
  
}




head(predict.all)

##### Marginalise results over parameter values for colonies and years ####

names(predict.all)
predict.all.c <- predict.all[!is.na(predict.all$p.fit),]
### here by calculating the mean
predict.marginalised <- plyr::ddply(predict.all.c,
                                    c("polygon.id", "Species", "Data_Year", "Data_Month"),
                                    summarise,
                                    p.fit.m = mean(p.fit, na.rm=T))#,
                                    #p.fit.m.se = sd(p.fit, na.rm=T)/sqrt(length(polygon.id)))
head(predict.marginalised)


##### Calculate mean probability for each month ####

predict.month <- plyr::ddply(predict.marginalised,
                             c("polygon.id", "Species", "Data_Month"),
                             summarise,
                             p.fit.mo = mean(p.fit.m, na.rm=T))#,
                             #p.fit.mo.se = sd(p.fit.m, na.rm=T)/sqrt(length(polygon.id)))

head(predict.month)


##### Sum presences for each year ####

# first calculate threshold of presence for each species

vars_narm <- c("Bath_90km", "Slope_90km", "Chl_log",
               "SST_90km_1day_mean", "SSTanom_90km_1day_mean", "SSTgrad_log",
               "windSp_90km_1day_mean", "ssh_90km_1day_mean", "eke_log",
               "windWSC_90km_7day_mean", "ColDist", "MPA_inoutf")


species <- c("BAPE", "BRNO", "LENO", "RTTR", "SOTE", "TRPE", "TRSH", "WTSH", "WTTR")

thresholds <- as.data.frame(matrix(ncol = 2, nrow = 0))
colnames(thresholds) <- c("Species", "Threshold")


for (sp in species){
  df <- as.data.frame(filter(all.sp, Species == sp))%>%
    drop_na(., any_of(vars_narm))
  mod <- get(paste0("mod_", sp))
  
  p2 <- predict(mod, type="response")
  p2_num <- as.numeric(p2)
  
  roccurve <- pROC::roc(df$use ~ as.numeric(p2))
  threshold = pROC::coords(roccurve, "best", ret = "threshold")[1,]
  
  t <- as.data.frame(bind_cols(sp, threshold))
  
  thresholds <- rbind(thresholds, t)
  print(sp)
}

colnames(thresholds) <- c("Species", "threshold")
predict.month$threshold <- thresholds[match(predict.month$Species, thresholds$Species), "threshold"]
predict.month$p.use <- ifelse(predict.month$p.fit.mo > predict.month$threshold, 1, 0)

head(predict.month)

### Sum presences per polygon

predict.sum <- plyr::ddply(predict.month,
                           c("polygon.id", "Species"),
                           summarise,
                           months.present = sum(p.use, na.rm= T))

head(predict.sum)

# if present for one or more months = present
predict.sum$p.sp <- ifelse(predict.sum$months.present >=1, 1, 0)
head(predict.sum)


##### Species richness per cell #####

predict.richness <- plyr::ddply(predict.sum,
                                c("polygon.id"),
                                summarise,
                                sp.richness = sum(p.sp, na.rm= T))
head(predict.richness)


####### ********************************** #######
####### Predict models w/ Chagos MPA only #######

species <- c("BAPE", "BRNO", "LENO",  "RTTR", "SOTE", "TRPE", "TRSH", "WTSH", "WTTR")

hex.env.coldist$MPA_inoutf <- hex.env.coldist$MPA_CHA_inoutf
head(hex.env.coldist)

predict.all <- as.data.frame(matrix(ncol = 8, nrow = 0))
colnames(predict.all) <- c("polygon.id", "Data_Year", "Data_Month", "Predict_Colony", "Predict_year", "p.fit", "p.fit.se", "Species")

for (sp in species){
  df <- as.data.frame(filter(all.sp, Species == sp))
  mod <- get(paste0("mod_", sp))
  
  df$Month <-  as.integer(format.Date(df$date, "%m"))
  
  ## subset environment to within species range & non-breeding months ##
  hex.predict <- subset(hex.env.coldist, hex.env.coldist[,sp] == 1 & Month %in% unique(df$Month))
  hex.predict$ColDist <- hex.predict[,paste0(sp, "_Coldist_BL_MLC")]
  head(hex.predict)
  
  ## standardise predict variables by use available df ##
  vars <- c("Bath_90km", "Slope_90km", "Chl_log",
            "SST_90km_1day_mean", "SSTanom_90km_1day_mean", "SSTgrad_log",
            "windSp_90km_1day_mean", "ssh_90km_1day_mean", "eke_log",
            "windWSC_90km_7day_mean", "ColDist")
  
  for (j in vars){
    hex.predict[,paste0(j, "_st")] <- (hex.predict[,j]-mean(df[,j], na.rm = T))/sd(df[,j], na.rm = T)
  }
  
  
  ### run for parameters at each colony & year [to marginalise] ##
  predict.sp <- as.data.frame(matrix(ncol = 8, nrow = 0))
  colnames(predict.sp) <- c("polygon.id", "Data_Year", "Data_Month", "Predict_Colony", "Predict_year", "p.fit", "p.fit.se", "Species")
  
  cols <- unique(df$Colony)
  yrs <- unique(df$year)
  
  for (col in cols){
    for (yr in yrs){
      
      hex.predict$Colony <- col
      hex.predict$year <- as.factor(yr)
      
      p1 <- mgcv::predict.gam(mod, newdata = hex.predict, type = "response", se.fit = T)
      hex.predict$p.fit <- as.data.frame(p1[1])[,"fit"]
      hex.predict$p.fit.se <- as.data.frame(p1[2])[,"se.fit"]
      
      predict.out <- hex.predict[,c("polygon.id", "Year", "Month", "Colony", "year", "p.fit", "p.fit.se")]
      colnames(predict.out) <- c("polygon.id", "Data_Year", "Data_Month", "Predict_Colony", "Predict_year", "p.fit", "p.fit.se")
      predict.out$Species <- sp 
      
      predict.sp <- rbind(predict.sp, predict.out)
      
      print(paste0(sp, ", ", col, ", ", yr))
    }
  }
  
  predict.all <- rbind(predict.all, predict.sp)
  
}

head(predict.all)

##### Marginalise results over parameter values for colonies and years ####

names(predict.all)
predict.all.c <- predict.all[!is.na(predict.all$p.fit),]
### here by calculating the mean
predict.marginalised <- plyr::ddply(predict.all.c,
                                    c("polygon.id", "Species", "Data_Year", "Data_Month"),
                                    summarise,
                                    p.fit.m = mean(p.fit, na.rm=T))#,
#p.fit.m.se = sd(p.fit, na.rm=T)/sqrt(length(polygon.id)))
head(predict.marginalised)


##### Calculate mean probability for each month ####

predict.month <- plyr::ddply(predict.marginalised,
                             c("polygon.id", "Species", "Data_Month"),
                             summarise,
                             p.fit.mo = mean(p.fit.m, na.rm=T))#,
#p.fit.mo.se = sd(p.fit.m, na.rm=T)/sqrt(length(polygon.id)))

head(predict.month)


##### Sum presences for each year ####

# first calculate threshold of presence for each species

vars_narm <- c("Bath_90km", "Slope_90km", "Chl_log",
               "SST_90km_1day_mean", "SSTanom_90km_1day_mean", "SSTgrad_log",
               "windSp_90km_1day_mean", "ssh_90km_1day_mean", "eke_log",
               "windWSC_90km_7day_mean", "ColDist", "MPA_inoutf")


species <- c("BAPE", "BRNO", "LENO", "RTTR", "SOTE", "TRPE", "TRSH", "WTSH", "WTTR")

thresholds <- as.data.frame(matrix(ncol = 2, nrow = 0))
colnames(thresholds) <- c("Species", "Threshold")


for (sp in species){
  df <- as.data.frame(filter(all.sp, Species == sp))%>%
    drop_na(., any_of(vars_narm))
  mod <- get(paste0("mod_", sp))
  
  
  p2 <- predict(mod, type="response")
  
  roccurve <- pROC::roc(df$use ~ as.numeric(p2))
  threshold = pROC::coords(roccurve, "best", ret = "threshold")[1,]
  
  t <- as.data.frame(bind_cols(sp, threshold))
  
  thresholds <- rbind(thresholds, t)
}

colnames(thresholds) <- c("Species", "threshold")
predict.month$threshold <- thresholds[match(predict.month$Species, thresholds$Species), "threshold"]
predict.month$p.use <- ifelse(predict.month$p.fit.mo > predict.month$threshold, 1, 0)

head(predict.month)

### Sum presences per polygon

predict.sum <- plyr::ddply(predict.month,
                           c("polygon.id", "Species"),
                           summarise,
                           months.present = sum(p.use, na.rm= T))

head(predict.sum)

# if present for one or more months = present
predict.sum$p.sp <- ifelse(predict.sum$months.present >=1, 1, 0)
head(predict.sum)


##### Species richness per cell #####

predict.richness <- plyr::ddply(predict.sum,
                                c("polygon.id"),
                                summarise,
                                sp.richness = sum(p.sp, na.rm= T))
head(predict.richness)

