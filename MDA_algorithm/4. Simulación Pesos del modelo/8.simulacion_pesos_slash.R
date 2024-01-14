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
#----------------------------------------------------------------------------
#genero la gamma truncada en el interavalo c(0,1)
gen.trun(par = c(1), family = "GA", name = "0TO1", type = "right")
#----------------------------------------------------------------------------
pesos_slash <- function(y_df, x_df, B_0, M_0, nu_0){
  d <- dista_maha_indiv(y_df, x_df, B_0, M_0) #distancia de Mahalanobis para 
  #cada individuo, respecto a los valores observados que tiene
  pesos <- c()
  for (i in 1:length(d)){
    p_i_k <-length(which(is.na( y_df[i,]) == F)) # numero de 
    #variables obsevadas
    alpha <- (p_i_k/2) + nu_0
    beta1 <- (2/d[i])
    pesos[i] <- rGA0TO1(1, mu = (alpha*beta1), sigma = (1/sqrt(alpha)))
  }
  orden <- patron_monotono(y_df, x_df)[[3]]
  pesos <- pesos[orden] #ordena los pesos de acuerdo al orden de los 
  #individuos en patrón monótono
  
  return(pesos)
}



