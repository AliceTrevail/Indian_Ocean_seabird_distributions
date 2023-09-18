#### Model validation #####
library(caret)
require(dplyr)
require(pROC)

# read data ####

all.sp <- read_csv("/Data/df_new_migration_dates_MCPcrop_restandardised.csv")

vars <- c("Bath_90km", "Slope_90km", "Chl_log",
          "SST_90km_1day_mean", "SSTanom_90km_1day_mean", "SSTgrad_log",
          "windSp_90km_1day_mean", "ssh_90km_1day_mean", "eke_log",
          "windWSC_90km_7day_mean", "ColDist")

all.sp_complete <- all.sp[complete.cases(all.sp[ , vars]),]
all.sp_incomplete <- all.sp[!complete.cases(all.sp[ , vars]),]


# read models ####
species <- c("BAPE", "BRNO", "LENO", "RTTR", "SOTE", "TRPE", "TRSH", "WTSH", "WTTR")

for (i in species){
  m <- readRDS(paste0("/Data/Models/Main/",i,"_main.rds"))
  assign(paste0("mod_", i), m, envir = .GlobalEnv)
}

#### create table for model diagnostics #####

species <- c("WTSH", "TRSH", "BAPE", "TRPE", "RTTR", "WTTR", "SOTE", "BRNO", "LENO")
diagnostics <- c("Species", "dev.expl", "AUC", "Threshold", "accuracy", "sensitivity", "specificity", 
                 "Pos.PP", "Neg.PP", "Precision", "Recall")

mod.valid <- as.data.frame(matrix(ncol = 11, nrow = 9))
colnames(mod.valid) <- diagnostics
mod.valid$Species <- species

#### extract model diagnostics in loop ####
i <- 1
for (i in 1:NROW(mod.valid)){
  
  sp <- mod.valid$Species[i]
  
  m <- get(paste0("mod_", sp))
  df <- filter(all.sp_complete, Species == sp)
  
  # deviance explained
  
  mod.valid$dev.expl[i] <- summary(m)$dev.expl*100
  
  
  # ROC curve, AUC and threshold 
  
  # predicted output from model, based on explanatory covariate values
  modpred <- predict(m, type="response")
  
  # calculate roc curve,
  roccurve <- pROC::roc(df$use, as.numeric(modpred))
  
  # calculate threshold in predictions where P(occurence) = presence
  t <- pROC::coords(roccurve, "best", ret = "threshold", transpose = F)[1,]
  mod.valid$Threshold[i] <- t
  
  # AUC = area under ROC curve 
  mod.valid$AUC[i] <-auc(roccurve)
  
  
  # calculate confusion matrix to understand ratios of true and false negatives and positives
  p <- as.numeric(predict(m , type="response")>t)
  df$use2 <- as.factor(df$use)
  p2 <- as.factor(p)
  cm <- confusionMatrix(p2,df$use2, positive = "1")
  
  
  # extract values from confusion matrix
  mod.valid$accuracy[i] <- cm$overall[[1]]*100
  mod.valid$sensitivity[i] <- cm$byClass[["Sensitivity"]]
  mod.valid$specificity[i] <- cm$byClass[["Specificity"]]
  mod.valid$Pos.PP[i] <- cm$byClass[["Pos Pred Value"]]
  mod.valid$Neg.PP[i] <- cm$byClass[["Neg Pred Value"]]
  mod.valid$Precision[i] <- cm$byClass[["Precision"]]
  mod.valid$Recall[i] <- cm$byClass[["Recall"]]
  
  
}



## ******************** ####
## Models including MPA ####


# read models ####
species <- c("BAPE", "BRNO", "LENO", "RTTR", "SOTE", "TRPE", "TRSH", "WTSH", "WTTR")

for (i in species){
  m <- readRDS(paste0("/Data/Models/MPA/",i,"_MPA.rds"))
  assign(paste0("mod_", i), m, envir = .GlobalEnv)
}

#### create table for model diagnostics #####

species <- c("WTSH", "TRSH", "BAPE", "TRPE", "RTTR", "WTTR", "SOTE", "BRNO", "LENO")
diagnostics <- c("Species", "dev.expl", "AUC", "Threshold", "accuracy", "sensitivity", "specificity", 
                 "Pos.PP", "Neg.PP", "Precision", "Recall")

mod.valid <- as.data.frame(matrix(ncol = 11, nrow = 9))
colnames(mod.valid) <- diagnostics
mod.valid$Species <- species

#### extract model diagnostics in loop ####
i <- 1
for (i in 1:NROW(mod.valid)){
  
  sp <- mod.valid$Species[i]
  
  m <- get(paste0("mod_", sp))
  df <- filter(all.sp_complete, Species == sp)
  
  # deviance explained
  
  mod.valid$dev.expl[i] <- summary(m)$dev.expl*100
  
  
  # ROC curve, AUC and threshold 
  
  # predicted output from model, based on explanatory covariate values
  modpred <- predict(m, type="response")
  
  # calculate roc curve,
  roccurve <- pROC::roc(df$use, as.numeric(modpred))
  
  # calculate threshold in predictions where P(occurence) = presence
  t <- pROC::coords(roccurve, "best", ret = "threshold", transpose = F)[1,]
  mod.valid$Threshold[i] <- t
  
  # AUC = area under ROC curve 
  mod.valid$AUC[i] <-auc(roccurve)
  
  
  # calculate confusion matrix to understand ratios of true and false negatives and positives
  p <- as.numeric(predict(m , type="response")>t)
  df$use2 <- as.factor(df$use)
  p2 <- as.factor(p)
  cm <- confusionMatrix(p2,df$use2, positive = "1")
  
  
  # extract values from confusion matrix
  mod.valid$accuracy[i] <- cm$overall[[1]]*100
  mod.valid$sensitivity[i] <- cm$byClass[["Sensitivity"]]
  mod.valid$specificity[i] <- cm$byClass[["Specificity"]]
  mod.valid$Pos.PP[i] <- cm$byClass[["Pos Pred Value"]]
  mod.valid$Neg.PP[i] <- cm$byClass[["Neg Pred Value"]]
  mod.valid$Precision[i] <- cm$byClass[["Precision"]]
  mod.valid$Recall[i] <- cm$byClass[["Recall"]]
  
  
}




