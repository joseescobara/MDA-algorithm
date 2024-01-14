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
#------------------------------------------------------------------------------
library(heavy)

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
source("4. Simulación Pesos del modelo/8.simulacion_pesos_slash.R")
#5.Cálculo de los Nu's del modelo 
source("5. Simulación Nu's del modelo/13.Simulacion_nu_posterior_slash_model.R")
#6.Valores iniciales 
source("6. Betas_cov Iniciales en caso de que haya Y_mis_des/4.funcion_datos_iniciales.R")
source("6. Betas_cov Iniciales en caso de que haya Y_mis_des/0.Betas_mcov_iniciales.R")
#aplicaciones de los modelos
source("7.EM_version_MDA_slash/aplicacion_modelo_slash.R")
source("10.ambientes_simulacion_modelo_slash/slash_model_ambiente_simulacion_n_150.R")

#-----------------------------------------------------------------------------
#----------------------------datos iniciales----------------------------------
#-----------------------------------------------------------------------------
#y_mp <- patron_monotono(respuestas, covariables)[[1]]
#x_mp <- patron_monotono(respuestas, covariables)[[2]]
#ln_ymp <- log(y_mp)
#set.seed(12)#semilla para tener los mismos datos imputados iniciales
#datos_init <- datos_iniciales(ln_ymp, x_mp)
#set.seed(12414) #ajustamos una semilla para obtener los mismos 
#valores iniciales en los ambientes de simulación
#parametros_iniciales <- betas_mvco_init(datos_init, x_mp)
#M_0 <- parametros_iniciales[[1]]
#B_0 <- parametros_iniciales[[2]]
#nu_0 <- 7 #el nu inicial lo tomamos como propuesta de la literatura
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
#set.seed(112413414) #para obtener los mismos parametros en el cual
#nos vamos a basar para simular la muestra
#system.time(aplicacion_slash <- sim_slash(ln_ymp, x_mp, B_0, M_0, nu_0, Nsim = 10000))["elapsed"]
#save(aplicacion_slash, file="aplicacion_modelo_slash.Rdata")
#load("aplicacion_modelo_slash.Rdata")
#aquí guardamos las medianas de los parámetros de la aplicación del modelo
#estim <- medianas(ln_ymp, aplicacion_slash)
#save(estim, file="parametros_reales_slash.Rdata")
load("parametros_reales_slash.Rdata")
#funcion para calcular las muestras

B_1 <- estim[[1]] #paramétros reales
M_1 <- estim[[2]] #paramétros reales
nu_1 <- estim[[3]] #paramétros reales

#-----------------------------------------------------------------------------
#-----------simulación de las 1000 muestras para el ambiente n = 150----------
#-----------------------------------------------------------------------------

#matriz_diseño_150_slash <- matriz_diseño(150)
#muestras <- muestra_slash_completa_150(matriz_diseño_150_slash, 
#                            B_1,  M_1, nu_1, N_muestras = 80000) 50000+50
#length(muestras)
#muestras_150_slash_faltante <- muestra_faltantes_slash_150(muestras,matriz_diseño_150_slash,
#                                             B_1 = B_1, M_1 = M_1)#aquí genero
#length(muestras_150_slash_faltante)
#las muestras con datos faltantes
#set.seed(986427)
#muestra_extraer <- sample(1:length(muestras_150_slash_faltante), 1000)
#muestras_150_slash_faltante <- muestras_150_slash_faltante[muestra_extraer]
#length(muestras_150_slash_faltante)
#save(muestras, matriz_diseño_150_slash, file = "muestras_150_slash_completa.Rdata")
#save(muestras_150_slash_faltante, matriz_diseño_150_slash, 
#     file = "muestras_150_slash_faltantes.Rdata")

#load("muestras_150_slash_completa.Rdata")
load("muestras_150_slash_faltantes.Rdata")
load("estimaciones_modelo_slash_150.Rdata")
#-----------------------------------------------------------------------------
#-------------Aplicción del modelo slash a las 1000 muestras------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

#lista_estimaciones_150_slash <- list()
set.seed(100020)
semillas3 <- sample(1:1000000, length(muestras_150_slash_faltante))
for (i in 121:length(muestras_150_slash_faltante)){
  set.seed(semillas3[i])#semilla
  y_mp3 <- patron_monotono(muestras_150_slash_faltante[[i]], matriz_diseño_150_slash)[[1]]
  x_mp3 <- patron_monotono(muestras_150_slash_faltante[[i]], matriz_diseño_150_slash)[[2]]
  estim_150_slash <- medianas_modelo_slash(y_mp3, x_mp3, B_0, M_0, nu_0, 
                                           Nsim = 10000)#aquí aplico a cada 
  #muestra las 10000 replicas y le saco la mediana.
  lista_estimaciones_150_slash[[i]] <- estim_150_slash
  save(lista_estimaciones_150_slash, file = "estimaciones_modelo_slash_150.Rdata")
}
