setwd("ruta")
#-----------------------------------------------------------------------------
#Datos para aplicar el algoritmo MDA (Comuna Robledo, Medellín)
#-----------------------------------------------------------------------------
load(r"ruta")
load("parametros_reales_modelo_t.Rdata")
load("parametros_reales_slash.Rdata")

library(mnormt)
#-----------------------------------------------------------------------------
#matriz de datos (respuestas y covariables)
#-----------------------------------------------------------------------------
datos <- base_de_datos #datos completos
datos <- as.matrix(datos)

#respuestas y covariables
respuestas <- datos[, 1:3] #matriz de variables respuesta
covariables <- cbind(rep(1, times = nrow(respuestas)), datos[,9], 
                    datos[, 6], datos[,8]) #matriz de covariables
colnames(covariables) <- c("X_0", "edad", "sexo", "tiempo_lechem")

log_datos <- log(respuestas)
#-----------------------------------------------------------------------------
#parametros modelo slash
coef.slash <- estim[[1]]
matdisp.slash <- estim[[2]]
nu.slash <- estim[[3]]

#parametros modelo t

coef.t <- estim_t[[1]]
matdisp.t <- estim_t[[2]]
nu.t <- estim_t[[3]]

#------------------------------------------------------------------------------
# PDF of the multivariate slash distribution
dmvslash <- function(x, mu, sigma, nu) {
  # x: response vector 
  # mu: location vector 
  # sigma: dispersion matrix
  # nu: tail parameter
  if (nu < 250) {
    integrand <- Vectorize(function(u) {
      mnormt::dmnorm(x * sqrt(u), 
                     mu * sqrt(u), sigma) * stats::dbeta(u, nu, 1)
    })
    dens.slash <- integrate(integrand, 0, 1)$value
  } else {
    integrand <- Vectorize(function(u) {
      mnormt::dmnorm(x * sqrt(u), mu * sqrt(u), sigma)
    })
    set.seed(123)
    u <- rbeta(3000, nu, 1)
    dens.slash <- mean(integrand(u))
  }
  return(dens.slash)
}


#-----------------------------------------------------------------------------


respect.observed.data <- function(y_df, x_df, coef, matdisp){
  lista_obserb <- list()
  M_1 <- matdisp #matriz de dispersion inicial.
  B_1 <- coef #matriz de coeficientes iniciales. 
  distancias <- c()
  for (i in 1:nrow(y_df)){#recorre todos los individuos
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
    lista_obserb[[i]] <- list(as.vector(y_i), as.vector(t(B_1)%*%x_df[i,]), M_1)
    #calcula la distancia de mahalanobis para el individuo i-ésimo
    #con su respectiva infomacion
    M_1 <- matdisp #vuelve a actualizar la matriz de disperison para la
    #siguiente iteración del individuo i+1
    B_1 <- coef #vuelve a actualizar la matriz de coeficientes para la
    #siguiente iteración del individuo i+1
  }
  return(lista_obserb)
}

#-----------------------------------------------------------------------------

bayes.factor.modelo <- function(y_df = log_datos, x_df = covariables, coef = coef.slash
,matdisp = matdisp.slash, nu = nu.slash, coef1 = coef.t, matdisp1 = matdisp.t,
nu1 = nu.t){
  
densidades.slash <- c()
densidades.t <- c()

datos.observed.slash <- respect.observed.data(y_df, x_df, coef, matdisp)
datos.observed.t <- respect.observed.data(y_df, x_df, coef1, matdisp1)


for (i in 1:nrow(y_df)){
  densidades.slash[i] <- dmvslash(datos.observed.slash[[i]][[1]], 
            datos.observed.slash[[i]][[2]], datos.observed.slash[[i]][[3]],nu)
}
  
for (j in 1:nrow(y_df)){
  densidades.t[j] <- dmt(x = datos.observed.t[[j]][[1]], 
                  mean = datos.observed.t[[j]][[2]], 
                  S = datos.observed.t[[j]][[3]],df = nu1)
}

BF <- log10(prod(densidades.t)/prod(densidades.slash))
return(BF)
}
                         

bayes.factor.modelo(y_df = log_datos, x_df = covariables, coef = coef.slash
                    ,matdisp = matdisp.slash, nu = nu.slash, coef1 = coef.t, matdisp1 = matdisp.t,
                    nu1 = nu.t)















