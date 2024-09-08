library(tidyverse)
library(stringr)
library(readxl)
library(writexl)
library(nlshrink)

source("utils.R")
setwd("/Data")

#Leitura dos dados
nyse_mensal <- read_excel("nyse_mensal.xlsx", skip = 3, col_types = c("text", rep("numeric", 1419)))
ibrx_mensal <- read_excel("ibrx_mensal.xlsx", skip = 3, col_types = c("text", rep("numeric", 97)))

#Renomeando colunas
colnames(nyse_mensal) <- c("Data", str_remove(colnames(nyse_mensal)[-1], "^(?s).*prov\n"))
colnames(ibrx_mensal) <- c("Data", str_remove(colnames(ibrx_mensal)[-1], "^(?s).*prov\n"))

#Removendo colunas com NA's
#Coluna Data no formato de data
nyse_mensal$Data <- str_glue("01-{nyse_mensal$Data}") %>% as.Date("%d-%b-%Y") 
ibrx_mensal$Data <- str_glue("01-{ibrx_mensal$Data}") %>% as.Date("%d-%b-%Y") 

lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
Sys.setlocale("LC_TIME", lct)

#Escolha da data
data <- as.Date(c('01-01-2000', '01-06-2023'), "%d-%m-%Y")

#Limpeza dos dados
nyse_mensal <- remove_cols(nyse_mensal, data)

ibrx_mensal <- remove_cols(ibrx_mensal, data)

cov_matriz_nyse <- cov(nyse_mensal[-1])
cov_matriz_ibrx <- ibrx_mensal %>% dplyr::select(-c(Data, IBOV)) %>% cov()

cov_matriz_nyse %>% as.data.frame() %>% write_xlsx("/Data/cov_matriz_nyse.xlsx")
cov_matriz_ibrx %>% as.data.frame() %>% write_xlsx("/Data/cov_matriz_ibrx.xlsx")

nls_matriz_nyse <- nlshrink_function(as.matrix(nyse_mensal[-1]))
nls_matriz_ibrx <- nlshrink_function(as.matrix(dplyr::select(ibrx_mensal, -c(Data, IBOV))))

nls_matriz_nyse %>% as.data.frame() %>% write_xlsx("/Data/nls_matriz_nyse.xlsx")
nls_matriz_ibrx %>% as.data.frame() %>% write_xlsx("/Data/nls_matriz_ibrx.xlsx")
