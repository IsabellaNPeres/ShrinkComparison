library(MASS)
library(nlshrink)
library(cvCovEst)
library(ShrinkCovMat)
library(stringr)
library(readxl)
library(dplyr)
library(RiskPortfolios)
library(tidyverse)
library(xtable)
source("utils.R")
files_names <- list.files(path = 'covShrinkage', pattern = '*.R', full.names = T)
for (name in files_names) source(name)


##########################################################
##                  Data cleaning                       ##
##########################################################

## Filter from 2000-01 to 2022-06
monthly_data <- read_xlsx("ibrx_mensal_aplicacao.xlsx", na = "-") |> 
  mutate(Data = str_replace(Data, "Jan", "01"), 
         Data = str_replace(Data, "Fev", "02"),
         Data = str_replace(Data, "Mar", "03"), 
         Data = str_replace(Data, "Abr", "04"),
         Data = str_replace(Data, "Mai", "05"), 
         Data = str_replace(Data, "Jun", "06"),
         Data = str_replace(Data, "Jul", "07"), 
         Data = str_replace(Data, "Ago", "08"),
         Data = str_replace(Data, "Set", "09"), 
         Data = str_replace(Data, "Out", "10"),
         Data = str_replace(Data, "Nov", "11"), 
         Data = str_replace(Data, "Dez", "12")) |> 
  mutate(Data = lubridate::my(Data)) |> 
  dplyr::filter(Data >= '2000-01-01', Data < '2022-07-01')

### select assets with no missing values
full_monthly_data <- monthly_data |> dplyr::select(names(which(apply(is.na(monthly_data), 2, sum) == 0)))

### extract ibov
ibov <- full_monthly_data |> dplyr::select(IBOV)

### assets without ibov
monthly_data <- full_monthly_data |> dplyr::select(-Data, -IBOV)
returns <- monthly_data

### Remove no ordinary assets or high correlated from the same company
ordinary_monthly_data <- monthly_data |> 
  dplyr::select(-BBDC4, -ELET6, -PETR4, -ITSA4, -GGBR4)

###################################################
##           Out-of-sample comparison            ##
###################################################

cov_est_names <- c("amos", "ls1", "ls2", "ls3", "ls4", "ls5", "nls_quest",
                   "nls_asym", "gis_est", "lis_est", "qis_est", "nld", 
                   "spikedF", "spikedOp", "spikedSt")

for (InS in c(60, 120)) {
  OoS <- nrow(returns) - InS
  
  # Portfolio weights
  w_mv <- list()
  
  # Portfolio returns
  Rp <- list()
  sspw <- list()
  to <- list()
  oos_table <- list()
  
  for (k in 1:15){
    w_mv[[k]] <-  matrix(NA, ncol = ncol(returns), nrow = OoS)
    Rp[[k]] <- matrix(NA, ncol = 8, nrow = OoS)
    sspw[[k]] <- matrix(NA, ncol = 8, nrow = OoS)
    to[[k]] <- matrix(NA, ncol = 8, nrow = OoS - 1)
  }
  
  for (i in 1:OoS) {
    print(sprintf("Window %i of %i", i, OoS))
    
    # In-sample returns
    ret <- as.matrix(returns[i:(i + InS - 1), ])
  
    # realized one-step-ahead returns
    realized_returns <- as.numeric(returns[i + InS, ]) 
  
    # Covariance matrix (input for portfolio allocation techniques) 
    amos <- cov(ret) 
    
    ls1 <- cov1Para(ret)  
    ls2 <- cov2Para(ret)  
    ls3 <- covCor(ret)  
    ls4 <- covDiag(ret)  
    ls5 <- shrinkcovmat.identity(t(ret))$Sigmahat
    
    
    nls_quest <- nlshrink_cov(ret) 
    nls_asym  <- nlShrinkLWEst(ret)
    gis_est       <- gis(ret)
    lis_est       <- lis(ret)
    qis_est       <- qis(ret)
    nld       <- denseLinearShrinkEst(ret)
    
    
    spikedF  <- spikedFrobeniusShrinkEst(ret, ncol(ret)/nrow(ret), num_spikes = NULL, noise = NULL)
    spikedOp <- spikedOperatorShrinkEst(ret, ncol(ret)/nrow(ret), num_spikes = NULL, noise = NULL)
    spikedSt <- spikedSteinShrinkEst(ret, ncol(ret)/nrow(ret), num_spikes = NULL, noise = NULL)
    
    cov_est <- list(amos, ls1, ls2, ls3, ls4, ls5, nls_quest, nls_asym,
                 gis_est,lis_est, qis_est, nld, spikedF, spikedOp, spikedSt)
  
    # Compute Portfolio Weights
    ## Portfolio allocation strategies implemented: Minimum Variance (MV)
    for (k in 1:15){
      w_mv[[k]][i, ] <- optimalPortfolio(Sigma = cov_est[[k]], control = list(type = 'minvol', constraint = 'lo'))
      
      # Compute Portfolio returns
      Rp[[k]][i, ] <- as.numeric(realized_returns %*% w_mv[[k]][i, ]) 
      
      # SSPW
      sspw[[k]][i, ] <- c(0, sum(w_mv[[k]][i, ]^2))
      
      # Turover
      if (i > 1) {
        to[[k]][i - 1, ] <- c(0,
                         calculate_to(w_mv[[k]][i - 1, ], w_mv[[k]][i, ], tail(ret, 1)))
      }
      
      
      oos_table[[k]] <- cbind(t(medidas(Rp[[k]][,1])), apply(to[[k]], 2, mean), apply(sspw[[k]], 2, mean))
      colnames(oos_table[[k]]) <- c("AV","SD", "SR", "ASR", "SO", "TO", "SSPW")
    }
    
  }
  
  # Save data
  for (k in 1:15){
  write.csv(Rp[[k]], paste0("Rp_", cov_est_names[[k]], "_", InS, ".csv"))
  write.csv(w_mv[[k]], paste0("w_mv_", cov_est_names[[k]], "_", InS, ".csv"))
  write.csv(oos_table[[k]], paste0("oos_", cov_est_names[[k]], "_", InS, ".csv"))
  }
  
}  

#######################################################
##                  Results reading                  ## 
#######################################################


files_names <- list.files(path = 'oos', pattern = '*.csv', full.names = T)

cov_est_names <- rep(c("amos", "ls1", "ls2", "ls3", "ls4", "ls5", "nls_quest",
                       "nls_asym", "gis_est", "lis_est", "qis_est", "nld", 
                       "spikedF", "spikedOp", "spikedSt"), each=2)

oos <- data.frame()

for (i in 1:length(files_names)) {
  oos_table <- read.csv(files_names[i])
  oos_table <- oos_table[, -1] 
  oos_table$Estimador <- cov_est_names[i]
  oos <- rbind(oos, oos_table)
}

oos$InS <- rep(c(60, 120), 15)

oos %>% filter(InS == 60) %>% select(-7) %>% 
  xtable(digits = 4)

oos %>% filter(InS == 120) %>% select(-7) %>% 
  xtable(digits = 4)

