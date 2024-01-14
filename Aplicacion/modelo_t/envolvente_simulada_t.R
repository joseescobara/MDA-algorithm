setwd("ruta")
#-----------------------------------------------------------------------------
#Datos para aplicar el algoritmo MDA (Comuna Robledo, Medellín)
#-----------------------------------------------------------------------------
load(r"ruta")
load("parametros_reales_modelo_t.rdata")
#-----------------------------------------------------------------------------
#librería a usar
#------------------------------------------------------------------------------
#library(stats)
#------------------------------------------------------------------------------
dista_maha_indiv <- function(y_df, x_df, B_0, M_0){
  M_1 <- M_0 #matriz de varianzas y covarianzas inicial.
  B_1 <- B_0 #matriz de coeficientes iniciales. 
  distancias <- c()
  for (i in 1:nrow(y_df)) {#recorre todos los individuos
    y_i <- y_df[i:i, ]
    mis <- which(is.na( y_df[i,]) == T ) #extrae los indices de los
    #faltantes del individuo i-ésimo. 
    if(length(mis) != 0  ){
      for (j in 1:length(mis)) {
        eliminar <- which(colnames(M_1) == names(mis[j]))#me devuelve 
        #la columna correspondiente en la iteración j-ésima 
        #a la cual esa variable es faltante para el individuo i-ésimo.
        #La cual será eliminada en esa iteración de la matriz de var_cov
        y_i_obs <-  y_i[names(y_i) != names(mis[j])] #eliminar en cada
        #iteracion las variables con observaciones faltantes para ese 
        #individuo, en caso de que haya 1 o más.
        sigma_i_k <- M_1[-(eliminar:eliminar), -(eliminar:eliminar)]
        #elimina las variables faltantes para el individuo i-ésimo en
        #cada iteración en caso de que tenga una variable no observada
        b_i_k <- B_1[, - eliminar] #elimina de la matriz de coeficientes
        #los coeficientes de las variables faltantes para 
        #el individuo i-ésimo
        M_1 <- sigma_i_k #actualizo para la siguiente iteración, por si hay 
        #otro dato faltante en esa observacion
        B_1 <- b_i_k #actualizo para la siguiente iteración, por si hay 
        #otro dato faltante en esa observacion
        y_i <- y_i_obs #actualizo para la siguiente iteración, por si hay 
        #otro dato faltante en esa observacion
      }
    }
    distancias[i] <- mahalanobis(y_i, t(B_1)%*%x_df[i:i,], M_1 )
    #calcula la distancia de mahalanobis para el individuo i-ésimo
    #con su respectiva infomacion
    M_1 <- M_0 #vuelve a actualizar la matriz de var_cov para la
    #siguiente iteración del individuo i+1
    B_1 <- B_0 #vuelve a actualizar la matriz de coeficientes para la
    #siguiente iteración del individuo i+1
  }
  return(distancias)
}

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

log_datos <- log(respuestas)
#------------------------------------------------------------------------------
coef.t <- estim_t[[1]]
matdisp.t <- estim_t[[2]]
nu.t <- estim_t[[3]]
n <- nrow(log_datos)
n <- nrow(log_datos)
p <- ncol(log_datos)
#------------------------------------------------------------------------------

distancias_t <- dista_maha_indiv(log_datos, covariables, coef.t, 
                               matdisp.t)

#-------------------------------------------------------------------------------
Mahadist.logt <- sort(distancias_t)
teoquant.logt <- p * qf((1:n/(n + 1)), df1 = p, df2 = nu.t)
e.logt <- matrix(0, n, 100)
e1.logt <- numeric(n)
e2.logt <- numeric(n)
set.seed(84303)
for (i in 1:100) {
  est.sam.logt <- mnormt::rmt(n, mean = rep(0, p), S = diag(p), 
                              df = nu.t)
  e.logt[, i] <- sort(apply(est.sam.logt^2, 1, sum))
}
for (i in 1:n) {
  eo.logt <- sort(e.logt[i, ])
  e1.logt[i] <- min(eo.logt)
  e2.logt[i] <- max(eo.logt)
}
#par(mfrow = c(1, 1), pty = "s")
plot(teoquant.logt, Mahadist.logt, xaxt = "n", yaxt = "n", xlab = "", 
     ylab = "", xlim = c(0, 24.4), ylim = c(0, 24.4), lty = 1, main = "", 
     pch = 19, cex = 1)
at_1.lt <- seq(from = 2, to = 24, by = 11)
at_2.lt <- seq(from = 2, to = 24, by = 11)
axis(side = 1, at = at_1.lt, cex.axis = 2.7, line = 0.1, tick = F)
axis(side = 2, at = at_2.lt, cex.axis = 2.7, line = -0.8, tick = F)
title(xlab = "Cuantiles teóricos", line = 3, cex.lab = 2.5)
title(ylab = "Cunatiles observados", line = 2.4, cex.lab = 2.5)
lines(x = c(0, 24.4), y = c(0, 24.4), col = rgb(0.2, 0.2, 0.2), 
      type = "l", lty = 2, lwd = 3)
par(new = TRUE)
plot(teoquant.logt, sort(e1.logt), axes = F, xlab = "", ylab = "", 
     type = "l", xlim = c(0, 24.4), ylim = c(0, 24.4), lwd = 3)
par(new = TRUE)
plot(teoquant.logt, sort(e2.logt), axes = F, xlab = "", ylab = "", 
     type = "l", xlim = c(0, 24.4), ylim = c(0, 24.4), lwd = 3)

