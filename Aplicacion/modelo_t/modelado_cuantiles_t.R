setwd("ruta")
#-----------------------------------------------------------------------------
#Datos para aplicar el algoritmo MDA (Comuna Robledo, Medellín)
#-----------------------------------------------------------------------------
load(r"ruta")#ruta base de datos real
load("parametros_reales_modelo_t.rdata")

coef.t <- estim_t[[1]]
matdisp.t <- estim_t[[2]]
nu.t <- estim_t[[3]]
#-----------------------------------------------------------------------------
#librería a usar
#------------------------------------------------------------------------------
#library(stats)
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
datos_1 <- datos_obs[-(which(is.na(datos_obs$per_bran) == T)), ]
datos_completos <- datos_1[-(which(is.na(datos_1$peso_act) == T)), ]
#-----------------------------------------------------------------------------
q.t.0.25 <- qt(0.25, df = nu.t)
q.t.0.50 <- qt(0.50, df = nu.t)
q.t.0.75 <- qt(0.75, df = nu.t)
sqrt.sigma11 <- sqrt(matdisp.t[1,1])
sqrt.sigma22 <- sqrt(matdisp.t[2,2])
sqrt.sigma33 <- sqrt(matdisp.t[3,3])
mean.tlm <- round(mean(datos_obs$tiempo_lechem), 2)
mean.edad <- round(mean(datos_obs$edad), 2)
#-----------------------------------------------------------------------------

# perimetro braquial (niña)

est.perbra.0.25.girl <- exp(coef.t[1,1] + coef.t[2,1]* mean.edad
                            + coef.t[4,1] * mean.tlm + sqrt.sigma11 * q.t.0.25)

est.perbra.0.50.girl <- exp(coef.t[1,1] + coef.t[2,1]* mean.edad
                            + coef.t[4,1] * mean.tlm + sqrt.sigma11 * q.t.0.50)

est.perbra.0.75.girl <- exp(coef.t[1,1] + coef.t[2,1]* mean.edad
                            + coef.t[4,1] * mean.tlm + sqrt.sigma11 * q.t.0.75)

quantile(datos_completos$per_bran, c(0.25, 0.50, 0.75))
#peso (niña)

est.peso.0.25.girl <- exp(coef.t[1,2] + coef.t[2,2]* mean.edad
                          + coef.t[4,2] * mean.tlm + sqrt.sigma22 * q.t.0.25)

est.peso.0.50.girl <-exp(coef.t[1,2] + coef.t[2,2]* mean.edad
                         + coef.t[4,2] * mean.tlm + sqrt.sigma22 * q.t.0.50)

est.peso.0.75.girl <- exp(coef.t[1,2] + coef.t[2,2]* mean.edad
                          + coef.t[4,2] * mean.tlm + sqrt.sigma22 * q.t.0.75)

quantile(datos_completos$peso_act, c(0.25, 0.50, 0.75))
#talla (niña)

est.talla.0.25.girl <- exp(coef.t[1,3] + coef.t[2,3]* mean.edad
                           + coef.t[4,3] * mean.tlm + sqrt.sigma33 * q.t.0.25)

est.talla.0.50.girl <- exp(coef.t[1,3] + coef.t[2,3]* mean.edad
                           + coef.t[4,3] * mean.tlm + sqrt.sigma33 * q.t.0.50)

est.talla.0.75.girl <- exp(coef.t[1,3] + coef.t[2,3]* mean.edad
                           + coef.t[4,3] * mean.tlm + sqrt.sigma33 * q.t.0.75)
quantile(datos_completos$talla_act, c(0.25, 0.50, 0.75))
#-------------------------------------------------------------------------------

# perimetro braquial (niño)

est.perbra.0.25.boy <- exp(coef.t[1,1] + coef.t[2,1]* mean.edad + coef.t[3,1]
                           + coef.t[4,1] * mean.tlm + sqrt.sigma11 * q.t.0.25)

est.perbra.0.50.boy <- exp(coef.t[1,1] + coef.t[2,1]* mean.edad + coef.t[3,1]
                           + coef.t[4,1] * mean.tlm + sqrt.sigma11 * q.t.0.50)

est.perbra.0.75.boy <- exp(coef.t[1,1] + coef.t[2,1]* mean.edad + coef.t[3,1]
                           + coef.t[4,1] * mean.tlm + sqrt.sigma11 * q.t.0.75)

#peso (niño)

est.peso.0.25.boy <- exp(coef.t[1,2] + coef.t[2,2]* mean.edad + coef.t[3,2]
                         + coef.t[4,2] * mean.tlm + sqrt.sigma22 * q.t.0.25)

est.peso.0.50.boy <-exp(coef.t[1,2] + coef.t[2,2]* mean.edad + coef.t[3,2]
                        + coef.t[4,2] * mean.tlm + sqrt.sigma22 * q.t.0.50)

est.peso.0.75.boy <- exp(coef.t[1,2] + coef.t[2,2]* mean.edad + coef.t[3,2]
                         + coef.t[4,2] * mean.tlm + sqrt.sigma22 * q.t.0.75)
#talla (niño)

est.talla.0.25.boy <- exp(coef.t[1,3] + coef.t[2,3]* mean.edad + coef.t[3,3]
                          + coef.t[4,3] * mean.tlm + sqrt.sigma33 * q.t.0.25)

est.talla.0.50.boy <- exp(coef.t[1,3] + coef.t[2,3]* mean.edad + coef.t[3,3]
                          + coef.t[4,3] * mean.tlm + sqrt.sigma33 * q.t.0.50)

est.talla.0.75.boy <- exp(coef.t[1,3] + coef.t[2,3]* mean.edad + coef.t[3,3]
                          + coef.t[4,3] * mean.tlm + sqrt.sigma33 * q.t.0.75)





