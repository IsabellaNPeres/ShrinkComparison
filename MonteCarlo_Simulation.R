library(clusterGeneration)
library(MASS)
library(nlshrink)
library(cvCovEst)
library(ShrinkCovMat)
library(readxl)
library(dplyr)
library(psych)
library(expm)
library(mnormt)
source("/utils.R")
files_names <- list.files(path = '/covShrinkage', pattern = '*.R', full.names = T)
for (name in files_names) source(name)

####################################################################
###                       Reading Data                          ####
####################################################################
## Scenario 1: Cov matrix simular do NYSE
#cov_matriz <- as.matrix(read_excel("../Data/cov_matriz_nyse.xlsx"))
## Scenario 2: Cov matrix simular do IBRX
# cov_matriz <- as.matrix(read_excel("cov_matriz_ibrx.xlsx"))


####################################################################
###                           Setup                             ####
####################################################################
MC <- 1000

#set.seed(MC)
#cov_matriz <- genPositiveDefMat("unifcorrmat", dim = 500)$Sigma
## Scenario 1: Cov matrix simular do NYSE
cov_matrix <- as.matrix(read_excel("/Data/ClusterGeneration/cov_matriz_clusterGeneration.xlsx"))
max_cols <- 250

for (n_cols in max_cols) {
#for (n_cols in c(10, 50, 100, max_cols)) {
  for (n_obs in c(500, 1000, 2000)) {
    ####################################################################
    ###                   Monte Carlo Simulation                    ####
    ####################################################################
    true_cov <- cov_matrix[1:n_cols, 1:n_cols]
    list_cov_matrix <- list()
    for (i in 1:MC) {
      print(i)
      set.seed(i)
      #a <- rmt(n = n_obs, mean = rep(0, n_cols), S = 5/7*true_cov, df = 7, sqrt=NULL)
      a <- mvrnorm(n_obs, mu = rep(0, nrow(true_cov)), Sigma = true_cov)
      is_error <- FALSE
      tryCatch({
        list_cov_matrix$amos[[i]]     <- cov(a) 
        
        list_cov_matrix$ls1[[i]]       <- cov1Para(a)  
        list_cov_matrix$ls2[[i]]       <- cov2Para(a)  
        list_cov_matrix$ls3[[i]]       <- covCor(a)  
        list_cov_matrix$ls4[[i]]       <- covDiag(a)  
        list_cov_matrix$ls5[[i]]       <- shrinkcovmat(t(a), target = 'identity')$Sigmahat

        list_cov_matrix$nls_quest[[i]]  <- nlshrink_cov(a) 
        list_cov_matrix$nls_asym[[i]]   <- nlShrinkLWEst(a)
        list_cov_matrix$gis[[i]]        <- gis(a)
        list_cov_matrix$lis[[i]]        <- lis(a)
        list_cov_matrix$qis[[i]]        <- qis(a)
        list_cov_matrix$nld[[i]]        <- denseLinearShrinkEst(a)

        
        list_cov_matrix$spikedF[[i]]  <- spikedFrobeniusShrinkEst(a, ncol(a)/nrow(a), num_spikes = NULL, noise = NULL)
        list_cov_matrix$spikedOp[[i]] <- spikedOperatorShrinkEst(a, ncol(a)/nrow(a), num_spikes = NULL, noise = NULL)
        list_cov_matrix$spikedSt[[i]] <- spikedSteinShrinkEst(a, ncol(a)/nrow(a), num_spikes = NULL, noise = NULL)
        
      },
      error = function(e) {
        is_error <- TRUE
      })
      while (is_error == TRUE) {
        is_error <- FALSE
        a <- mvrnorm(n_obs, mu = rep(0, nrow(true_cov)), Sigma = true_cov)
        #a <- rmt(n = n_obs, mean = rep(0, n_cols), S = 5/7*true_cov, df = 7, sqrt=NULL)
        tryCatch({
          list_cov_matrix$amos[[i]]     <- cov(a) 
          
          list_cov_matrix$ls1[[i]]       <- cov1Para(a)  
          list_cov_matrix$ls2[[i]]       <- cov2Para(a)  
          list_cov_matrix$ls3[[i]]       <- covCor(a)  
          list_cov_matrix$ls4[[i]]       <- covDiag(a)  
          list_cov_matrix$ls5[[i]]       <- shrinkcovmat(a, target = 'identity')$Sigmahat
          
          list_cov_matrix$nls_quest[[i]]  <- nlshrink_cov(a) 
          list_cov_matrix$nls_asym[[i]]   <- nlShrinkLWEst(a)
          list_cov_matrix$gis[[i]]        <- gis(a)
          list_cov_matrix$lis[[i]]        <- lis(a)
          list_cov_matrix$qis[[i]]        <- qis(a)
          list_cov_matrix$nld[[i]]        <- denseLinearShrinkEst(a)
          
          list_cov_matrix$spikedF[[i]]  <- spikedFrobeniusShrinkEst(a, ncol(a)/nrow(a), num_spikes = NULL, noise = NULL)
          list_cov_matrix$spikedOp[[i]] <- spikedOperatorShrinkEst(a, ncol(a)/nrow(a), num_spikes = NULL, noise = NULL)
          list_cov_matrix$spikedSt[[i]] <- spikedSteinShrinkEst(a, ncol(a)/nrow(a), num_spikes = NULL, noise = NULL)

        },
        error = function(e) {
          is_error <- TRUE
        })
      }
    }
    name_estim <- names(list_cov_matrix)
    n_estim <- length(name_estim)
    
    
    ####################################################################
    ###                    Performance Evaluation                   ####
    ####################################################################
    rmse_var <- rmse_cov <- frobenius <- inverse_stein <- 
      min_var <- stein <- inverse_frobenius <- symmetrized_stein <-
      weighted_frobenius <- disutility <- frechet <- quadratic <- 
      inverse_quadratic <- matrix(NA, ncol = n_estim, nrow = MC)
    
    colnames(rmse_var) <- colnames(rmse_cov) <- colnames(frobenius) <- colnames(inverse_stein) <- 
      colnames(min_var) <- colnames(stein) <- colnames(inverse_frobenius) <-
      colnames(symmetrized_stein) <- colnames(weighted_frobenius) <- colnames(disutility) <- 
      colnames(frechet) <- colnames(quadratic) <- colnames(inverse_quadratic) <- name_estim
    
    ones_matrix <- matrix(rep(1, nrow(true_cov)*ncol(true_cov)), ncol = ncol(true_cov))
    
    for (i in 1:MC) {
      for (j in 1:n_estim) {
        estim_cov <- as.matrix(list_cov_matrix[[j]][[i]])
        inverse_estim_matrix <- solve(estim_cov)
        inverse_true_matrix <- solve(true_cov)
        aux <- (estim_cov - true_cov)
        
        rmse_var[i, j] <- sqrt(mean(diag(aux)^2))
        rmse_cov[i, j] <- sqrt(mean(aux[upper.tri(aux, diag = FALSE)]^2))
        frobenius[i, j] <- norm(aux, "F") / n_cols
        inverse_stein[i, j] <- (tr(inverse_estim_matrix %*% true_cov) - log(det(inverse_estim_matrix %*% true_cov))) / n_cols
        min_var[i, j] <- tr(inverse_estim_matrix %*% true_cov %*% inverse_estim_matrix) / tr(inverse_estim_matrix)^2 * n_cols
        stein[i, j] <- (tr(estim_cov %*% inverse_true_matrix) - log(det(estim_cov %*% inverse_true_matrix))) / n_cols
        inverse_frobenius[i, j] <- norm(inverse_estim_matrix - inverse_true_matrix, "F") / n_cols 
        symmetrized_stein[i, j] <- (tr(estim_cov %*% inverse_true_matrix +  inverse_estim_matrix %*% true_cov))/ n_cols 
        weighted_frobenius[i, j] <- (tr((estim_cov - true_cov) %*% (estim_cov - true_cov)%*% inverse_true_matrix)) / n_cols 
        disutility[i, j] <- (tr((inverse_estim_matrix - inverse_true_matrix) %*% (inverse_estim_matrix - inverse_true_matrix) %*% true_cov)) / n_cols 
        frechet[i, j] <- norm(sqrtm(estim_cov) - sqrtm(true_cov), "F") / n_cols 
        quadratic[i, j] <- norm(inverse_true_matrix %*% estim_cov - ones_matrix, "F") / n_cols 
        inverse_quadratic[i, j] <- norm(inverse_estim_matrix %*% true_cov - ones_matrix, "F") / n_cols 
      }
    }
    
    
    ####################################################################
    ###                        Saving Results                       ####
    ####################################################################
    write.csv(rmse_var, paste0("/Data/ClusterGeneration/Student-T/rmse_var_", n_obs, "_", n_cols, ".csv"))
    write.csv(rmse_cov, paste0("/Data/ClusterGeneration/Student-T/rmse_cov_", n_obs, "_", n_cols, ".csv"))
    write.csv(frobenius, paste0("/Data/ClusterGeneration/Student-T/frobenius_", n_obs, "_", n_cols, ".csv"))
    write.csv(inverse_stein, paste0("/Data/ClusterGeneration/Student-T/inverse_stein_", n_obs, "_", n_cols, ".csv"))
    write.csv(min_var, paste0("/Data/ClusterGeneration/Student-T/min_var_", n_obs, "_", n_cols, ".csv"))
    write.csv(stein, paste0("/Data/ClusterGeneration/Student-T/stein_", n_obs, "_", n_cols, ".csv"))
    write.csv(inverse_frobenius, paste0("/Data/ClusterGeneration/Student-T/inverse_frobenius_", n_obs, "_", n_cols, ".csv"))
    write.csv(symmetrized_stein, paste0("/Data/ClusterGeneration/Student-T/symmetrized_stein_", n_obs, "_", n_cols, ".csv"))
    write.csv(weighted_frobenius, paste0("/Data/ClusterGeneration/Student-T/weighted_frobenius_", n_obs, "_", n_cols, ".csv"))
    write.csv(disutility, paste0("/Data/ClusterGeneration/Student-T/disutility_", n_obs, "_", n_cols, ".csv"))
    write.csv(frechet, paste0("/Data/ClusterGeneration/Student-T/frechet_", n_obs, "_", n_cols, ".csv"))
    write.csv(quadratic, paste0("/Data/ClusterGeneration/Student-T/quadratic_", n_obs, "_", n_cols, ".csv"))
    write.csv(inverse_quadratic, paste0("/Data/ClusterGeneration/Student-T/inverse_quadratic_", n_obs, "_", n_cols, ".csv"))
    
  }
}

rm(list=ls())
