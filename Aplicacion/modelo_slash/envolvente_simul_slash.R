setwd("ruta")
#-----------------------------------------------------------------------------
#Datos para aplicar el algoritmo MDA (Comuna Robledo, Medellín)
#-----------------------------------------------------------------------------
load(r"ruta") #ruta base de datos real
load("parametros_reales_slash.Rdata")
library(heavy)
#----------------------distancia de Mahalanobis-------------------------------
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
log_datos <- log(respuestas)

#-----------------------------------------------------------------------------
coef.slash <- estim[[1]]
matdisp.slash <- estim[[2]]
nu.slash <- estim[[3]]
n <- nrow(log_datos)
p <- ncol(log_datos)
#------------------------------------------------------------------------------

distancias_slash <- dista_maha_indiv(log_datos, covariables, coef.slash, 
                               matdisp.slash)

#_-----------------------------------------------------------------------------
# Quantile function of the distribution of the Mahalanobis 
# distances of multivariate slash random vector
qmslash <- function(alpha, p, nu) {
  # alpha: probability
  # p: number of response variables
  # nu: tail parameter
  fdamaha <- function(s) {
    return(pchisq(s, df = p) - (((2^nu) * gamma((p/2) + nu))/((s^nu) * gamma(p/2))) * 
             pchisq(s, df = (p + 2 * nu)) - alpha)
  }
  qtl <- uniroot(fdamaha, c(1e-10, 1e+11))$root
  return(qtl)
}


#=========================================================
# Quantile-quantile plots of Mahalanobis distances with 
# simulated envelope for the bivariate log-normal,      
# log-t and log-slash linear regression models          
#=========================================================
#-----------------------------------------------------------------------------
#------------------------------------------------------------------------------
Mahadist.slash <- sort(distancias_slash)
quant.mslash <- Vectorize(function(alpha) {
  qmslash(alpha = alpha, p = p, nu = nu.slash)
})
teoquant.logslash <- quant.mslash((1:n)/(n + 1))
e.logslash <- matrix(0, n, 100)
e1.logslash <- numeric(n)
e2.logslash <- numeric(n)
set.seed(179795)
for (i in 1:100) {
  est.sam.logslash <- rmslash(n, center = rep(0, p), Scatter = diag(p), 
                              df = nu.slash)
  e.logslash[, i] <- sort(apply(est.sam.logslash^2, 1, sum))
}
for (i in 1:n) {
  eo.logslash <- sort(e.logslash[i, ])
  e1.logslash[i] <- min(eo.logslash)
  e2.logslash[i] <- max(eo.logslash)
}
#par(mfrow = c(1, 1), pty = "s")
plot(teoquant.logslash, Mahadist.slash, xaxt = "n", yaxt = "n", 
     xlab = "", ylab = "", xlim = c(0, 28.6), ylim = c(0, 110.6), 
     lty = 1, main = "", pch = 19, cex = 1)
at_1.ls <- seq(from = 2, to = 28, by = 13)
at_2.ls <- seq(from = 2, to = 110, by = 54)
axis(side = 1, at = at_1.ls, cex.axis = 2.7, line = 0.1, tick = F)
axis(side = 2, at = at_2.ls, cex.axis = 2.7, line = -0.8, tick = F)
title(xlab = "Cuantiles teóricos", line = 3, cex.lab = 2.5)
title(ylab = "Cuantiles observados", line = 2.4, cex.lab = 2.5)
lines(x = c(0, 28.6), y = c(0, 28.6), col = rgb(0.2, 0.2, 0.2), type = "l", 
      lty = 2, lwd = 3)
par(new = TRUE)
plot(teoquant.logslash, sort(e1.logslash), axes = F, xlab = "", 
     ylab = "", type = "l", xlim = c(0, 28.6), ylim = c(0, 110.6), 
     lwd = 3)
par(new = TRUE)
plot(teoquant.logslash, sort(e2.logslash), axes = F, xlab = "", 
     ylab = "", type = "l", xlim = c(0, 28.6), ylim = c(0, 110.6), 
     lwd = 3)


