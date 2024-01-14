setwd("ruta")
#-----------------------------------------------------------------------------
#Datos para aplicar el algoritmo MDA (Comuna Robledo, Medellín)
#-----------------------------------------------------------------------------
load(r"ruta")
load("parametros_reales_slash.Rdata")

coef.slash <- estim[[1]]
matdisp.slash <- estim[[2]]
nu.slash <- estim[[3]]
#-----------------------------------------------------------------------------
#librería a usar
#------------------------------------------------------------------------------
#install.packages("aplpack")
#install.packages("mnormt")
#install.packages("dplyr")
#install.packages("quantreg")
#install.packages("pracma")
#library(dplyr)
#library(ggplot2)
#library(aplpack)
#library(pracma)
#library(quantreg)
#library(devtools)
#library(mnormt)

#necesario para generar los pesos del modelo Slash
#install.packages("gamlss")
#install.packages("gamlss.dist")
#install.packages("gamlss.tr")
#library(gamlss)
#library(gamlss.tr)
#library(gamlss.dist)
#install_version("heavy", version="0.38.19")
#necesario para generar los datos iniciales
#library(norm)
#library(norm2)
#necesario para aplicar el proceso de imputación
library(mnorm)
#para medio vectorizar  a una matriz
library(fBasics)
library(mvtnorm)
library(heavy)
library(MASS)
#------------------------------------------------------------------------------


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

datos_obs <- data.frame(respuestas, covariables) #conjunto de datos observados
#incluyendo las covariables
#datos_1 <- datos_obs[-(which(is.na(datos_obs$per_bran) == T)), ]
#datos_completos <- datos_1[-(which(is.na(datos_1$peso_act) == T)), ]
#------------------------------------------------------------------

# Quantile function of the univariate slash distribution
qslash <- function(p, mu, sigma2, nu, m){
  # p: probability
  # mu: location parameter
  # sigma2: dispersion parameter
  # nu: tail parameter
  # m: independent draws for Monte Carlo integration
  qtl.func <- function(w) {
    set.seed(123)
    u <- rbeta(m, nu, 1)
    integrand <- Vectorize(function(u) {
      stats::pnorm(q = w, mean = mu, sd = sqrt(sigma2/u))
    })
    return(mean(integrand(u)) - p)
  }
  return(uniroot(qtl.func, c(-1e+16,1e+16))$root)
}

#------------------------calculo_cuantiles_slash------------------------------
q.stdslash_0.25 <- qslash(p = 0.25, mu = 0, sigma2 = 1, 
                          nu = nu.slash, m = 3000)
q.stdslash_0.50 <- 0
q.stdslash_0.75 <- -q.stdslash_0.25
sqrt.sigma11 <- sqrt(matdisp.slash[1,1])
sqrt.sigma22 <- sqrt(matdisp.slash[2,2])
sqrt.sigma33 <- sqrt(matdisp.slash[3,3])
mean.tlm <- round(mean(datos_obs$tiempo_lechem), 2)
mean.edad <- round(mean(datos_obs$edad), 2)
#-----------------------------------------------------------------------------

# perimetro braquial (niña)

est.perbra.0.25.girl <- exp(coef.slash[1,1] + coef.slash[2,1]* mean.edad
+ coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.25)

est.perbra.0.50.girl <- exp(coef.slash[1,1] + coef.slash[2,1]* mean.edad
+ coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.50)

est.perbra.0.75.girl <- exp(coef.slash[1,1] + coef.slash[2,1]* mean.edad
    + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.75)

quantile(datos_completos$per_bran, c(0.25, 0.50, 0.75))
#peso (niña)

est.peso.0.25.girl <- exp(coef.slash[1,2] + coef.slash[2,2]* mean.edad
                       + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * q.stdslash_0.25)

est.peso.0.50.girl <-exp(coef.slash[1,2] + coef.slash[2,2]* mean.edad
                    + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * q.stdslash_0.50)

est.peso.0.75.girl <- exp(coef.slash[1,2] + coef.slash[2,2]* mean.edad
                    + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * q.stdslash_0.75)

quantile(datos_completos$peso_act, c(0.25, 0.50, 0.75))
#talla (niña)

est.talla.0.25.girl <- exp(coef.slash[1,3] + coef.slash[2,3]* mean.edad
                     + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.25)

est.talla.0.50.girl <- exp(coef.slash[1,3] + coef.slash[2,3]* mean.edad
                     + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.50)

est.talla.0.75.girl <- exp(coef.slash[1,3] + coef.slash[2,3]* mean.edad
                     + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.75)
quantile(datos_completos$talla_act, c(0.25, 0.50, 0.75))
#-------------------------------------------------------------------------------

# perimetro braquial (niño)

est.perbra.0.25.boy <- exp(coef.slash[1,1] + coef.slash[2,1]* mean.edad + coef.slash[3,1]
     + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.25)

est.perbra.0.50.boy <- exp(coef.slash[1,1] + coef.slash[2,1]* mean.edad + coef.slash[3,1]
                       + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.50)

est.perbra.0.75.boy <- exp(coef.slash[1,1] + coef.slash[2,1]* mean.edad + coef.slash[3,1]
                       + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.75)

#peso (niño)

est.peso.0.25.boy <- exp(coef.slash[1,2] + coef.slash[2,2]* mean.edad + coef.slash[3,2]
                     + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * q.stdslash_0.25)

est.peso.0.50.boy <-exp(coef.slash[1,2] + coef.slash[2,2]* mean.edad + coef.slash[3,2]
                    + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * q.stdslash_0.50)

est.peso.0.75.boy <- exp(coef.slash[1,2] + coef.slash[2,2]* mean.edad + coef.slash[3,2]
                     + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * q.stdslash_0.75)
#talla (niño)

est.talla.0.25.boy <- exp(coef.slash[1,3] + coef.slash[2,3]* mean.edad + coef.slash[3,3]
                      + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.25)

est.talla.0.50.boy <- exp(coef.slash[1,3] + coef.slash[2,3]* mean.edad + coef.slash[3,3]
                      + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.50)

est.talla.0.75.boy <- exp(coef.slash[1,3] + coef.slash[2,3]* mean.edad + coef.slash[3,3]
                      + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.75)





