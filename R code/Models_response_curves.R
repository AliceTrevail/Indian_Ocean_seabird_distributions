library(mgcv)
library(tidyverse)
options(tibble.width = Inf)
options(tibble.width = getOption("width"))


### Read in data ###

all.sp <- read_csv("/Data/df_new_migration_dates_MCPcrop.csv")%>%
  mutate(wt = ifelse(use == "1",30, 1))

### Standardise ###

fn_standardise <- function(x){
  (x - mean(x, na.rm = T))/sd(x, na.rm = T)
}
  
vars <- c("Bath_90km", "Slope_90km", "Chl_log",
          "SST_90km_1day_mean", "SSTanom_90km_1day_mean", "SSTgrad_log",
          "windSp_90km_1day_mean", "ssh_90km_1day_mean", "eke_log",
          "windWSC_90km_7day_mean", "ColDist")

vars_st <- c("Bath_90km_st", "Slope_90km_st", "Chl_log_st",
          "SST_90km_1day_mean_st", "SSTanom_90km_1day_mean_st", "SSTgrad_log_st",
          "windSp_90km_1day_mean_st","ssh_90km_1day_mean_st","eke_log_st",
          "windWSC_90km_7day_mean_st","ColDist_st")

all.sp2 <- all.sp %>%
  select(-all_of(vars_st)) %>%
  group_by(Species) %>%
  mutate(across(all_of(vars), 
                .fns = list(st = ~fn_standardise(.x))))

write_csv(all.sp2, "/Data/df_new_migration_dates_MCPcrop_restandardised.csv")


##### Models ####
# set up formula and function to run GAM
formula_RSF <- as.formula(use ~ s(Bath_90km_st, k = 4) + s(Slope_90km_st, k = 4) +
                    s(SST_90km_1day_mean_st, k = 4) + s(SSTanom_90km_1day_mean_st, k = 4) + s(SSTgrad_log_st, k = 4) +
                    s(Chl_log_st, k = 4) + s(windSp_90km_1day_mean_st, k = 4) + s(windWSC_90km_7day_mean_st, k = 4) +
                    s(eke_log_st, k = 4) + s(ssh_90km_1day_mean_st, k = 4) + ColDist_st + as.factor(year))

function_runRSF <- function(dataframe){
  m1<-gam(formula_RSF, family=binomial, weights = wt, data = dataframe)
  m1
}

formula_RSF_col <- as.formula(use ~ s(Bath_90km_st, k = 4) + s(Slope_90km_st, k = 4) +
                            s(SST_90km_1day_mean_st, k = 4) + s(SSTanom_90km_1day_mean_st, k = 4) + s(SSTgrad_log_st, k = 4) +
                            s(Chl_log_st, k = 4) + s(windSp_90km_1day_mean_st, k = 4) + s(windWSC_90km_7day_mean_st, k = 4) +
                            s(eke_log_st, k = 4) + s(ssh_90km_1day_mean_st, k = 4) + ColDist_st + as.factor(year) + Colony)


function_runRSF_col <- function(dataframe){
  m1<-gam(formula_RSF_col, family=binomial, weights = wt, data = dataframe)
  m1
}

species_multiplecols <- c("RTTR", "SOTE", "WTSH", "WTTR")
species_singlecols <- c("BAPE", "BRNO", "LENO", "TRPE", "TRSH")

for (i in species_singlecols){
  d_i <- filter(all.sp, Species == i)
  m <- function_runRSF(dataframe = d_i)
  saveRDS(m, paste0("/Data/Models/Main/", i, "_main.rds"))
  rm(m)
  print(i)
}

for (i in species_multiplecols){
  d_i <- filter(all.sp, Species == i)
  m <- function_runRSF_col(dataframe = d_i)
  saveRDS(m, paste0("/Data/Models/Main/", i, "_main.rds"))
  rm(m)
  print(i)
}


##### Extract response curves ####
species <- c("BAPE", "BRNO", "LENO",  "RTTR", "SOTE", "TRPE", "TRSH", "WTSH", "WTTR")

results <- as.data.frame(matrix(ncol =5, nrow = 0))
colnames(results) <- c("Variable", "x", "fit", "se", "Species")


par(mfrow = c(2,5))

for (i in species){
  
  mod <- readRDS(paste0("/Data/Models/Main/", i, "_main.rds"))
  df <- filter(all.sp, Species == i)
  
  pd <- plot(mod)
  
  for (j in 1:length(pd)){
    
    pd1 <- pd[[j]]
    
    pd1.2 <- data.frame(pd1$xlab, pd1$x, pd1$fit, pd1$se)
    colnames(pd1.2) <- c("Variable", "x", "fit", "se")
    pd1.2$Species <- i
    head(pd1.2)
    
    results <- rbind.data.frame(results, pd1.2)
    
    cd <- (termplot(mod, data = df, envir = environment(formula(mod)),
                    terms = c("ColDist_st"), se = T, plot = F))$ColDist_st
    
    colnames(cd) <- c("x", "fit", "se")
    cd$Variable <- "ColDist_st"
    cd$Species <- i
    
    results <- rbind.data.frame(results, cd)
    
    
  }
  
  rm(mod)
  rm(df)
  
  print(i)
}

resp.curves <- results
resp.curves$Variable <- sub("_.*", "", resp.curves$Variable)
unique(resp.curves$Variable)
str(resp.curves)

### un-standardise x axes ####
vars <- c("Bath_90km", "Slope_90km", "Chl_log",
          "SST_90km_1day_mean", "SSTanom_90km_1day_mean", "SSTgrad_log",
          "windSp_90km_1day_mean", "ssh_90km_1day_mean", "eke_log",
          "windWSC_90km_7day_mean", "ColDist")

# mean of explanatory variables from original data
env.meansd <- all.sp %>%
  group_by(Species) %>%
  summarise(across(all_of(vars), list(mean = ~mean(., na.rm = T), sd = ~sd(., na.rm = T)))) %>%
  pivot_longer(!Species, names_to = c("Variable", "Type"), values_to = c("value"), names_pattern = "(.*)_(\\w+)")%>%
  mutate(Variable = sub("_.*", "", Variable)) %>%
  pivot_wider(names_from = Type, values_from = value)

resp.curves_env <- resp.curves %>%
  rename(x_st = x) %>%
  left_join(env.meansd) %>%
  mutate(x_t = (x_st * sd)+mean) %>%  # back transform
  mutate(x = ifelse(Variable %in% c("Chl", "eke", "SSTgrad"), # inverse log
                    exp(x_t), x_t))


#### Response curves plot ####


head(resp.curves_env[resp.curves_env$Variable == "ColDist",])

resp.curves_env$upper <- resp.curves_env$fit + resp.curves_env$se
resp.curves_env$lower <- resp.curves_env$fit - resp.curves_env$se

resp.curves_env$fit_t <- exp(resp.curves_env$fit)/(1+exp(resp.curves_env$fit))
resp.curves_env$upper_t <- exp(resp.curves_env$upper)/(1+exp(resp.curves_env$upper))
resp.curves_env$lower_t <- exp(resp.curves_env$lower)/(1+exp(resp.curves_env$lower))

resp.curves_env$VariableLong <- as.factor(resp.curves_env$Variable)
levels(resp.curves_env$VariableLong)
levels(resp.curves_env$VariableLong) <- c("Bathymetry", "Chlorophyll", "Colony Distance", "Eddy kinetic energy", "Slope", "Sea surface height",
                                       "SST", "SST anomaly", "SST gradient", "Wind speed", "Wind stress curl")
resp.curves_env$VariableLong <- fct_relevel(resp.curves_env$VariableLong, "Bathymetry", "Slope", "SST", "SST anomaly", "SST gradient", "Chlorophyll", "Wind speed", "Wind stress curl", "Sea surface height", "Eddy kinetic energy")


resp.curves_env$Species <- fct_relevel(resp.curves_env$Species, "WTSH", "TRSH", "BAPE", "TRPE", "RTTR", "WTTR", "SOTE", "BRNO", "LENO")
Splabs <- c("Wedge-tailed\nshearwater", 
            "Tropical\nshearwater", 
            "Barau's petrel", 
            "Trindade petrel", 
            "Red-tailed\ntropicbird", 
            "White-tailed\ntropicbird", 
            "Sooty tern", 
            "Brown noddy", 
            "Lesser noddy")

cs <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', '#882255', '#AA4499', '#DDDDDD')


plot1 <- ggplot(resp.curves_env, aes(x = x, y = fit_t, col = Species))+
  facet_wrap(.~VariableLong, scales = "free_x")+
  geom_line(linewidth = 0.5)+
  geom_line(aes(x = x, y = upper_t), lty = "dashed", linewidth = 0.4)+
  geom_line(aes(x = x, y = lower_t), lty = "dashed", linewidth = 0.4)+
  scale_color_manual(values = cs, labels = Splabs)+
  scale_y_continuous(limits = c(0,1))+
  labs(x= "Value", y = "Probability of occurence")+
  theme_bw()+
  theme(legend.key.height = unit(1.8, "line"),  
        legend.text = element_text(margin = margin(l = 30, unit = "pt")))

ggsave(plot = plot1, filename = "/Figures/Manuscript/Supp mat/EnvResponseCurve.png", width = 8.5, height = 6, units = "in")

plot2 <- ggplot(resp.curves_env, aes(x = x, y = fit_t, col = Species))+
  facet_wrap(.~VariableLong, scales = "free_x")+
  geom_line(linewidth = 0.5)+
  geom_line(aes(x = x, y = upper_t), lty = "dashed", linewidth = 0.4)+
  geom_line(aes(x = x, y = lower_t), lty = "dashed", linewidth = 0.4)+
  scale_color_manual(values = cs, labels = Splabs)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_log10()+
  labs(x= "Value", y = "Probability of occurence")+
  theme_bw()+
  theme(legend.key.height = unit(1.8, "line"),  
        legend.text = element_text(margin = margin(l = 30, unit = "pt")))

ggsave(plot = plot2, filename = "/Figures/Manuscript/Supp mat/EnvResponseCurve_log.png", width = 8.5, height = 6, units = "in")


##### Models with MPA ####
# set up formula and function to run GAM
formula_RSF_withMPA <- as.formula(use ~ s(Bath_90km_st, k = 4) + s(Slope_90km_st, k = 4) +
                            s(SST_90km_1day_mean_st, k = 4) + s(SSTanom_90km_1day_mean_st, k = 4) + s(SSTgrad_log_st, k = 4) +
                            s(Chl_log_st, k = 4) + s(windSp_90km_1day_mean_st, k = 4) + s(windWSC_90km_7day_mean_st, k = 4) +
                            s(eke_log_st, k = 4) + s(ssh_90km_1day_mean_st, k = 4) + ColDist_st + as.factor(year) + MPA_inoutf )

function_runRSF_withMPA <- function(dataframe){
  m1<-gam(formula_RSF_withMPA, family=binomial, weights = wt, data = dataframe)
  m1
}

formula_RSF_withMPA_col <- as.formula(use ~ s(Bath_90km_st, k = 4) + s(Slope_90km_st, k = 4) +
                                s(SST_90km_1day_mean_st, k = 4) + s(SSTanom_90km_1day_mean_st, k = 4) + s(SSTgrad_log_st, k = 4) +
                                s(Chl_log_st, k = 4) + s(windSp_90km_1day_mean_st, k = 4) + s(windWSC_90km_7day_mean_st, k = 4) +
                                s(eke_log_st, k = 4) + s(ssh_90km_1day_mean_st, k = 4) + ColDist_st + as.factor(year) + MPA_inoutf + Colony)


function_runRSF_withMPA_col <- function(dataframe){
  m1<-gam(formula_RSF_withMPA_col, family=binomial, weights = wt, data = dataframe)
  m1
}

species_multiplecols <- c("RTTR", "SOTE", "WTSH", "WTTR")
species_singlecols <- c("BAPE", "BRNO", "LENO", "TRPE", "TRSH")

for (i in species_singlecols){
  d_i <- filter(all.sp, Species == i)
  m <- function_runRSF_withMPA(dataframe = d_i)
  saveRDS(m, paste0("/Data/Models/MPA/", i, "_MPA.rds"))
  rm(m)
  print(i)
}

for (i in species_multiplecols){
  d_i <- filter(all.sp, Species == i)
  m <- function_runRSF_withMPA_col(dataframe = d_i)
  saveRDS(m, paste0("/Data/Models/MPA/", i, "_MPA.rds"))
  rm(m)
  print(i)
}

