#vamos a crear el ambiente de simulacion para el modelo 
#slash con n = 50. Los valores iniciales se proponen completando los datos faltantes
#por medio del algoritmo MDA propuesto por schafer(1997), luego
#utilizamos los pesos iguales a uno, suponiendo normalidad de los
#datos. Luego, corremos una iteración para extraer la matriz de 
#coeficiente y la matriz de varcov para utilizarlos como datos
#iniciales
setwd("C:/Users/Jose/OneDrive - Universidad de Antioquia/Escritorio/Maestria/Tesis/Algoritmo MDA/Version_10_MDA")

#-----------------------------------------------------------------------------
#Datos para aplicar el algoritmo MDA (Comuna Robledo, Medellín)
#-----------------------------------------------------------------------------
#load(r"(C:\Users\Jose\OneDrive - Universidad de Antioquia\Escritorio\Maestria\Tesis\Algoritmo MDA\Version_9 MDA\db_1)")

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
library(mnorm)
#para medio vectorizar  a una matriz
library(fBasics)
library(mvtnorm)
library(MASS)
#------------------------------------------------------------------------------


#-----------------------------------------------------------------------------
#matriz de datos (respuestas y covariables)
#---------------------------------------------------------------------------
#datos <- base_de_datos #datos completos
#datos <- as.matrix(datos)

#respuestas y covariables
#respuestas <- datos[, 1:3] #matriz de variables respuesta
#covariables <- cbind(rep(1, times = nrow(respuestas)), datos[,9], 
#                     datos[, 6], datos[,8]) #matriz de covariables
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
source("11.ambientes_simulacion_modelo_t/t_model_ambiente_simulacion_n_100.R")

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

#load("valores_iniciales.Rdata")
#B_0 <- valores_iniciales[[1]]
#M_0 <- valores_iniciales[[2]]
#nu_0 <- valores_iniciales[[3]]
#----------------------------------------------------------------------------
#-----------------------paramétros verdaderos---------------------------------
#-----------------------------------------------------------------------------
#set.seed(112413414) #para obtener los mismos parametros en el cual
#nos vamos a basar para simular la muestra
#system.time(aplicacion_t <- sim_t(ln_ym, x_mp, B_0, M_0,Nsim = 10000))["elapsed"]
#save(aplicacion_t, file = "aplicacion_modelo_t.Rdata")
#load("aplicacion_modelo_t.Rdata")
#aquí calculamos las medianas de la aplicación del modelo.
#estim_t <- medianas(ln_ym, aplicacion_t)
#save(estim_t, file = "parametros_reales_modelo_t.Rdata")
#load("parametros_reales_modelo_t.Rdata")
#funcion para calcular las muestras

#B_1_t <- estim_t[[1]]
#M_1_t <- estim_t[[2]]
#nu_1_t <- estim_t[[3]]
#----------------------------------------------------------------------

#matriz_diseño_100_t <- matriz_diseño(n = 100)
#muestras <- muestra_t_completa_100(m_dis = matriz_diseño_100_t, 
#            B_1 = B_1_t, M_1 = M_1_t ,nu_1 =  nu_1_t, N_muestras = 13000)
#load("muestras_100_t_completa.Rdata")
#length(muestras)
#muestra_tamaño_100_t_faltante <- muestra_faltantes_t_100(muestras, 
#          matriz_diseño_100_t, B_1 = B_1_t, M_1 = M_1_t)
#length(muestra_tamaño_100_t_faltante)
#set.seed(1071287)
#muestra_extraer <- sample(1:length(muestra_tamaño_100_t_faltante), 1000)
#muestra_tamaño_100_t_faltante <- muestra_tamaño_100_t_faltante[muestra_extraer]
#save(muestras, matriz_diseño_100_t, file = "muestras_100_t_completa.Rdata")
#save(muestra_tamaño_100_t_faltante, matriz_diseño_100_t, 
#     file = "muestras_100_t_faltantes.Rdata")
load("muestras_100_t_faltantes.Rdata")
#load("estimaciones_100_modelo_t.Rdata")
#-----------------------------------------------------------------------------

lista_estimaciones_100_t <- list()
set.seed(5021444)
semillas5 <- sample(1:1241240, length(muestra_tamaño_100_t_faltante))
for (i in 1:length(muestra_tamaño_100_t_faltante)){
  set.seed(semillas5[i])
  y_mp5 <- patron_monotono(muestra_tamaño_100_t_faltante[[i]], matriz_diseño_100_t)[[1]]
  x_mp5 <- patron_monotono(muestra_tamaño_100_t_faltante[[i]], matriz_diseño_100_t)[[2]]
  estim_100_t <- medianas_modelo_t(y_mp5, x_mp5, B_0, M_0, Nsim = 10000) 
  #aquí aplico la función medianas_modelo_t para a cada muestra simulada 
  #aplicarle las 10000 replicas MC y sacarle la mediana. 
  lista_estimaciones_100_t[[i]] <- estim_100_t
  save(lista_estimaciones_100_t, file = "estimaciones_100_modelo_t.Rdata")
}

