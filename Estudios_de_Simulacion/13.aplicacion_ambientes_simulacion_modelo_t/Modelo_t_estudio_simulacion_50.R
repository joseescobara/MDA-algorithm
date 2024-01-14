#vamos a crear el ambiente de simulacion para el modelo 
#slash con n = 50. Los valores iniciales se proponen completando los datos faltantes
#por medio del algoritmo MDA propuesto por schafer(1997), luego
#utilizamos los pesos iguales a uno, suponiendo normalidad de los
#datos. Luego, corremos una iteración para extraer la matriz de 
#coeficiente y la matriz de varcov para utilizarlos como datos
#iniciales

setwd("ruta")
#-----------------------------------------------------------------------------
#Datos para aplicar el algoritmo MDA (Comuna Robledo, Medellín)
#-----------------------------------------------------------------------------
#load(r"ruta")

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
library(norm2)
#necesario para aplicar el proceso de imputación
library(mnorm)
#para medio vectorizar  a una matriz
library(fBasics)
library(mvtnorm)
#------------------------------------------------------------------------------


#-----------------------------------------------------------------------------
#matriz de datos (respuestas y covariables)
#---------------------------------------------------------------------------
#datos <- base_de_datos #datos completos
#datos <- as.matrix(datos)

#respuestas y covariables
#respuestas <- datos[, 1:3] #matriz de variables respuesta
#covariables <- cbind(rep(1, times = nrow(respuestas)), datos[,9], 
 #                    datos[, 6], datos[,8]) #matriz de covariables
#colnames(covariables) <- c("X_0", "edad", "sexo", "tiempo_lechem")

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
source("4. Simulación Pesos del modelo/7.simulacion_pesos_t.R")
#5.Cálculo de los Nu's del modelo 
source("5. Simulación Nu's del modelo/12.Simulacion_nu_posterior_t_model.R")
#6.Valores iniciales 
source("6. Betas_cov Iniciales en caso de que haya Y_mis_des/4.funcion_datos_iniciales.R")
source("6. Betas_cov Iniciales en caso de que haya Y_mis_des/0.Betas_mcov_iniciales.R")
#aplicaciones de los modelos
source("8.ECME_version_MDA_t/aplicacion_modelo_t.R")
source("11.ambientes_simulacion_modelo_t/t_model_ambiente_simulacion_n_50.R")
#----------------------------------------------------------------------------
#datos iniciales
#----------------------------------------------------------------
#y_mp <- patron_monotono(respuestas, covariables)[[1]]
#x_mp <- patron_monotono(respuestas, covariables)[[2]]
#ln_ym <- log(y_mp)
#datos_init <- datos_iniciales(ln_ym, x_mp)
#set.seed(12414) #ajustamos una semilla para obtener los mismos 
#valores iniciales en los ambientes de simulación
#parametros_iniciales <- betas_mvco_init(datos_init, x_mp)
#M_0 <- parametros_iniciales[[1]]
#B_0 <- parametros_iniciales[[2]]
#nu_0 <- 7 #el nu inicial lo tomamos como propuesta de la literatura

load("valores_iniciales.Rdata")
B_0 <- valores_iniciales[[1]]
M_0 <- valores_iniciales[[2]]
nu_0 <- valores_iniciales[[3]]
#----------------------------------------------------------------------------
#-----------------------paramétros verdaderos---------------------------------
#-----------------------------------------------------------------------------
#set.seed(112413414) #para obtener los mismos parametros en el cual
#nos vamos a basar para simular la muestra
#system.time(aplicacion_t <- sim_t(ln_ym, x_mp, B_0, M_0,Nsim = 10000))["elapsed"]
#save(aplicacion_t, file = "aplicacion_modelo_t.Rdata")
#load("aplicacion_modelo_t.Rdata")
#aquí calculamos las medianas de la aplicación del modelo.
#estim_t <- medianas(ln_mp, aplicacion_t)
#save(estim_t, file = "parametros_reales_modelo_t.Rdata")
#load("parametros_reales_modelo_t.Rdata")
#funcion para calcular las muestras

#B_1_t <- estim_t[[1]]
#M_1_t <- estim_t[[2]]
#nu_1_t <- estim_t[[3]]
#----------------------------------------------------------------------

#matriz_diseño_50_t <- matriz_diseño(n = 50)
#muestras <- muestra_t_completa_50(m_dis = matriz_diseño_50_t, 
#                            B_1 = B_1_t, M_1 = M_1_t ,nu_1_t, 6000)
#load("muestras_50_t_completa.Rdata")
#length(muestras)
#muestra_tamaño_50_t_faltantes <- muestra_faltantes_t_50(muestras, 
#matriz_diseño_50_t, B_1 = B_1_t, M_1 = M_1_t)
#length(muestra_tamaño_50_t_faltantes)
#set.seed(933124)
#muestras_extraer <- sample(1:length(muestra_tamaño_50_t_faltantes), 1000, replace = F)
#muestra_tamaño_50_t_faltantes <- muestra_tamaño_50_t_faltantes[muestras_extraer]
#save(muestras, matriz_diseño_50_t, file = "muestras_50_t_completa.Rdata")
#save(muestra_tamaño_50_t_faltantes, matriz_diseño_50_t, file = 
#"muestras_50_t_faltantes.Rdata")
load("muestras_50_t_faltantes.Rdata")
#load("estimaciones_50_modelo_t.Rdata")



#-----------------------------------------------------------------------------

lista_estimaciones_50_t <- list()

#------------------------------------------------------------------------

set.seed(121004)
semillas4 <- sample(1:1241240, length(muestra_tamaño_50_t_faltantes))
for (i in 1:length(muestra_tamaño_50_t_faltantes)){
  set.seed(semillas4[i])
  y_mp4 <- patron_monotono(muestra_tamaño_50_t_faltantes[[i]], matriz_diseño_50_t)[[1]]
  x_mp4 <- patron_monotono(muestra_tamaño_50_t_faltantes[[i]], matriz_diseño_50_t)[[2]]
  estim_50_t <- medianas_modelo_t(y_mp4, x_mp4, B_0, M_0, Nsim = 10000) 
  #aquí aplico la función medianas_modelo_t para a cada muestra simulada 
  #aplicarle las 10000 replicas MC y sacarle la mediana. 
  lista_estimaciones_50_t[[i]] <- estim_50_t
  save(lista_estimaciones_50_t, file = "estimaciones_50_modelo_t.Rdata")
  }
 

#-----------------------------------------------------------------------

simulaciones <- lista_estimaciones_50_t
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












