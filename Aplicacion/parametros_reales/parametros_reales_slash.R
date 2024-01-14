
#iniciales
setwd("ruta")
#-----------------------------------------------------------------------------
#Datos para aplicar el algoritmo MDA (Comuna Robledo, Medellín)
#-----------------------------------------------------------------------------
load(r"ruta")#base de datos real

#-----------------------------------------------------------------------------
#librería a usar
#------------------------------------------------------------------------------
library(readxl) 
library(beepr)
library(ggplot2)
library(tidyverse)

#necesario para generar los pesos del modelo Slash
#install.packages("gamlss")
#install.packages("gamlss.dist")
#install.packages("gamlss.tr")
library(gamlss)
library(gamlss.tr)
library(gamlss.dist)

#necesario para generar los datos iniciales
library(norm)
#library(norm2)
#necesario para aplicar el proceso de imputación
library(MASS)
library(mnorm)
#para medio vectorizar  a una matriz
library(fBasics)
#------------------------------------------------------------------------------
#library(heavy)
library(VGAM)

#-----------------------------------------------------------------------------
#matriz de datos (respuestas y covariables)
#---------------------------------------------------------------------------
datos <- base_de_datos #datos completos
datos <- as.matrix(datos)

#respuestas y covariables
respuestas <- datos[, 1:3] #matriz de variables respuesta
covariables <- cbind(rep(1, times = nrow(respuestas)), datos[,9], 
                    datos[, 6], datos[,8]) #matriz de covariables
colnames(covariables) <- c("X_0", "edad", "sexo", "tiempo_lechem")

#----------------------------------------------------------------------------


#------------------------------------------------------------------------------
#funciones del algoritmo MDA
#------------------------------------------------------------------------------
#1.ordenación
source("1.ordenación/1.funcion_patron_monotono.R")
source("1.ordenación/2.funcion_comienzo_patrones.R")
source("1.ordenación/3.funciones_matriz_por_patron.R")
#2.Cálculo betas_Covs (bajo los patrones y los del modelo)
source("2.Cálculo Betas-Covs/5.Estimacion_betas_cov.R")
source("2.Cálculo Betas-Covs/10.Simulacion_matriz_dist_posterior.R")
source("2.Cálculo Betas-Covs/11.simulacion_betas_posterior.R")
#3.imputación
source("3.imputación/6.imputacion_datos_destruyen_patron.R")
#4.cálculo de los pesos para los 3 modelos 
source("4. Simulación Pesos del modelo/8.simulacion_pesos_slash.R")
#5.Cálculo de los Nu's del modelo 
source("5. Simulación Nu's del modelo/13.Simulacion_nu_posterior_slash_model.R")
#6.Valores iniciales 
source("6. Betas_cov Iniciales en caso de que haya Y_mis_des/4.funcion_datos_iniciales.R")
source("6. Betas_cov Iniciales en caso de que haya Y_mis_des/0.Betas_mcov_iniciales.R")
#aplicaciones de los modelos
source("7.EM_version_MDA_slash/aplicacion_modelo_slash.R")
source("10.ambientes_simulacion_modelo_slash/slash_model_ambiente_simulacion_n_50.R")

#-----------------------------------------------------------------------------
#----------------------------datos iniciales----------------------------------
#-----------------------------------------------------------------------------
y_mp <- patron_monotono(respuestas, covariables)[[1]]
x_mp <- patron_monotono(respuestas, covariables)[[2]]
ln_ymp <- log(y_mp)
set.seed(12)#semilla para tener los mismos datos imputados iniciales
datos_init <- datos_iniciales(ln_ymp, x_mp)
set.seed(12414) #ajustamos una semilla para obtener los mismos 
#valores iniciales en los ambientes de simulación
parametros_iniciales <- betas_mvco_init(datos_init, x_mp)
M_0 <- parametros_iniciales[[1]]
B_0 <- parametros_iniciales[[2]]
nu_0 <- 7 #el nu inicial lo tomamos como propuesta de la literatura
#-----------------------------------------------------------------
#valores_iniciales <- list(B_0, M_0, nu_0)
#save(valores_iniciales, file="valores_iniciales.Rdata")

load("valores_iniciales.Rdata")
B_0 <- valores_iniciales[[1]]
M_0 <- valores_iniciales[[2]]
nu_0 <- valores_iniciales[[3]]
#-----------------------------------------------------------------------------
#----------------paramétros verdaderos modelo slash---------------------------
#-----------------------------------------------------------------------------
set.seed(112413414) #para obtener los mismos parametros en el cual
#nos vamos a basar para simular la muestra
system.time(aplicacion_slash <- sim_slash(ln_ymp, x_mp, B_0, M_0, nu_0, Nsim = 10000))["elapsed"]
save(aplicacion_slash, file="aplicacion_modelo_slash.Rdata")
#load("aplicacion_modelo_slash.Rdata")
#aquí guardamos las medianas de los parámetros de la aplicación del modelo


#función para el calculo de la mediana posterior de los coeficientes 
#del modelo slash
medianas <- function(y_mp, simulaciones){
  #obtenidas
  Nsim = length(simulaciones)
  matriz_coeficientes <- list()
  matriz_varianzas <- list()
  nu <- c()
  for (i in 1:Nsim) {
    matriz_coeficientes[[i]] <- simulaciones[[i]][[1]] #matriz de betas 
    matriz_varianzas[[i]] <- simulaciones[[i]][[2]] #matriz de varianzas y covarianzas
    nu[i] <- simulaciones[[i]][[3]]
  }
  
  A <- matrix(0, ncol = nrow(as.matrix(matriz_coeficientes[[1]]))*
                ncol(as.matrix(matriz_coeficientes[[1]])), 
              nrow = length(matriz_coeficientes))
  
  B <- matrix(0, ncol = (ncol(y_mp)*(ncol(y_mp)+1))/2,
              nrow = length(matriz_varianzas))
  colnames(A) <- c("B_01", "B_11", "B_21", "B_31",
                   "B_02", "B_12", "B_22", "B_32",
                   "B_03", "B_13", "B_23", "B_33")
  colnames(B) <- c("sigma_11", "sigma_21", "sigma_31", "sigma_22",
                   "sigma_32", "sigma_33")
  for (i in 1:length(matriz_coeficientes)) {
    coef <-  vec(matriz_coeficientes[[i]])
    A[i, ] <- coef
  }
  
  for (i in 1:length(matriz_varianzas)) {
    mcov <-  vech(matriz_varianzas[[i]])
    B[i, ] <- mcov
  }
  
  medianas_coeficientes <- apply(A, 2, median)
  median_var <- apply(B, 2, median)
  coef_estimados <- t(matrix(medianas_coeficientes, byrow = T, nrow = 3))
  var_estimadas <- matrix(c(median_var[1], median_var[2], median_var[3],
                            median_var[2], median_var[4], median_var[5],
                            median_var[3], median_var[5], median_var[6]), 
                          byrow = T, nrow = 3)
  colnames(coef_estimados) <- colnames(matriz_coeficientes[[1]])
  rownames(coef_estimados) <- rownames(matriz_coeficientes[[1]])
  
  colnames(var_estimadas) <- colnames(y_mp)
  rownames(var_estimadas) <- colnames(y_mp)
  
  nu_est <- median(nu)
  estimaciones <- list(coef_estimados, var_estimadas, nu_est)
  return(estimaciones)
}

#-----estimativas reales------------------------------------------
estim <- medianas(ln_ymp, aplicacion_slash) #cálculamos la mediana
#save(estim, file="parametros_reales_slash.Rdata")
#load("parametros_reales_slash.Rdata")
#funcion para calcular las muestras

B_1 <- estim[[1]] #paramétros reales
M_1 <- estim[[2]] #paramétros reales
nu_1 <- estim[[3]] #paramétros reales