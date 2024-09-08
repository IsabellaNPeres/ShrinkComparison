####################################
###     Auxiliary Functions      ###
####################################


#Função pra remover conlunas 
remove_cols <- function(df, data){
  df <- df %>% filter(between(Data, data[1], data[2]))
  df <- df[, apply(df, 2, anyNA) == F]
  df_cor <- df[-1]
  cor_matriz  <- cor(df_cor) 
  cor_matriz[upper.tri(cor_matriz)] <- 0
  diag(cor_matriz) <- 0
  df_cor <- df_cor[ , apply(cor_matriz, 2, function(x) any(abs(x) < 0.95))]
  df <- cbind(df[, 1], df_cor)
  return(df)
}

#Função transformar D's em data frame
D_dataframe <- function(D, estimadores, MC, n){
  df <- data_frame(n = 1:MC)
  for (i in 1:length(D)){
    df <- cbind(df, data_frame(D = D[[i]])) 
  }
  df <- df[, -1]
  colnames(df) <- estimadores
  i <- 1:ncol(df)
  df[, i] <-  apply(df[ , i], 2, as.numeric)
  df$n <- n
  return(df)
}

#Função TryCatch
nlshrink_function <- function(dados){
  tryCatch(
    #try to do this
    {
      #some expression
      estimativa <- nlshrink_cov(dados)
      return(estimativa)
    },
    #if an error occurs, tell me the error
    error=function(e) {
      return(NA)
    },
    #if a warning occurs, tell me the warning
    warning=function(w) {
      message('A Warning Occurred')
      print(w)
    }
  )
}

#Funções para comparar estiamdores na aplicação
medidas <- function(x, rf = 0.5) {
  # Annualized Average
  AV <- mean(x)
  # Annualized SD
  SD <- sd(x)
  # Information (or Sharpe) Ratio
  SR <- (mean(x) - rf)/sd(x)
  # Adjusted Sharpe Ratio
  ASR <- SR*(1 + (moments::skewness(x)/6)*SR - ((moments::kurtosis(x) - 3)/24)*SR^2)
  # Sortino Ratio
  SO <- (mean(x) - rf)/sqrt(mean(ifelse(x - rf < 0, 0, (x - rf)^2)))
  output <- c(12*AV, sqrt(12)*SD, sqrt(12)*SR, sqrt(12)*ASR, sqrt(12)*SO)
  return(output)
}


calculate_to <- function(previous_weights, desired_weights, oos_returns) {
  num <- previous_weights * (1 + oos_returns/100)
  den <- sum(num)
  updated_weights <- num/den
  to <- sum(abs(desired_weights - updated_weights))
  return(to)
}