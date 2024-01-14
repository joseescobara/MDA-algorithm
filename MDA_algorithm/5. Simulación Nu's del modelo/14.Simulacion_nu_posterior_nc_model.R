#--------------------------------------------------------------------
#funcion para calcular las distancias de mahalanobis
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
        M_1 <- sigma_i_k
        B_1 <- b_i_k
        y_i <- y_i_obs
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
#----------------------------------------------------------------------

#funcion de densidad para el parámetro lambda del modelo
#normal contaminado

nu_modelo_normal_conta <- function(y_df = y_df, x_df = x_df, 
                                   B_0 = B_0, M_0 = M_0, gama = gama, lambda){
distan <- dista_maha_indiv(y_df, x_df, B_0, M_0) #calculamos las 
#distancias de mahalanobis
densidad <- 1 #la densidad comienza en 1
for (i in 1:length(distan)) {
  p_i_k <-length(which(is.na( y_df[i,]) == F )) #calcula la dimensión
  # de las variables observadas del individuo i_ésimo
  densidad_i <- gama*(lambda^(p_i_k/2))*exp((-lambda*distan[i])/2) +
    (1-gama)*exp((-distan[i])/2) #calcula la densidad para el 
  #individuo i_ésimo
  densidad <- densidad_i*densidad #calcula la productoria
}  
return(densidad)
}

#------------------------------------------------------------------------------
#función para muestrar lamnda
simulacion_lambda_nc <- function(y_df = y_df, x_df = x_df, 
                                 B_0 = B_0, M_0 = M_0, gama = gama){
lambda <- seq(0,1, 0.002)#creamos la rejilla  
densidad <- nu_modelo_normal_conta(y_df, x_df, B_0, M_0, gama, lambda) 
#evaluamos cada valor de la rejilla en la respectiva densidad
suma_densidades <- sum(densidad)
probabilidad <- densidad/suma_densidades
muestra_lambda <- sample(lambda, 1, replace = T, prob = probabilidad)
return(muestra_lambda)
}

#------------------------------------------------------------------------------

#funcion para extraer el gama del modelo normal contaminado


simulacion_gamma_nc <- function(a = 1, b = 1, pesos = pesos, lambda = lambda){
  a_1 <- sum((1-pesos)/(1-lambda))
  b_1 <- sum((pesos-lambda)/(1-lambda)) 
  #calculo del la segunda sumatoria
  muestra_gama <- rbeta(1, a + a_1, b + b_1)
  return(muestra_gama)
}

#------------------------------------------------------------------------------


