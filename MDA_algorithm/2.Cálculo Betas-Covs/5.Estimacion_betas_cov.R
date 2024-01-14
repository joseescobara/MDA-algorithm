
#incluir matriz con los NA. Lo ajusté para que los datos entraran ordenados
individuos_de_acuerdo_patron <- function(y_mp){
  columnas_patron <- numero_columna(y_mp) #utilizo la función número_columnas
  #para verificar en cuáles columnas empizan los patrones reales.
  indivi_por_patron <- matriz_patrones(y_mp)[[2]] #utilizo la función 
  #matriz_patrones para extraer el número de individuos por patrones reales.
  r <- rep(0, times = ncol(y_mp)) #creo este vector para agregar los n_j 
  #a los p posibles posibles, j = 1, ..., p
  j <- 0 #utilizo este j para indexar a los individuos por patrón y
  #agregarlos a sus patrones reales.
  for (i in columnas_patron){ #aquí recorro las columnas con patrones reales
    j <- j +1 #
    r[i] <- indivi_por_patron[j] #reemplazo a los individuos en sus patrones
    #reales, los demás patrones quedan asignados con 0 individuos.
  }
  return(r)
}

#------------------------------------------------------------------------
#función para el cálculo de los betas estimados
#aquí la matriz Y tiene que entrar sin datos que destruyan el patron 
#monótono, y los xmp tambien organizada de acuerdo al patron. 

estimacion_betas_Bk <- function(y_mp, x_mp, pesos){
  #La función de esta parte es analizar si hay datos que destruyen el patrón 
  #monótono, en ese caso el algoritmo devuelve error
  matriz_por_patrones <- matriz_patrones(y_mp)[[1]] #aplico la función matriz
  #por patrones para extraer las matrices por patrón, y luego las analizo
  #para rectificar si tienen datos que destruyen el patrón monótono
  for (detector_datos_NA in matriz_por_patrones) {
    if(any(is.na(detector_datos_NA)) == T){
      return(c("Error, hay datos faltantes que destruyen el patrón monótono"))
      break
    }
  }
  numero_individuos <- individuos_de_acuerdo_patron(y_mp)#utilizo esta función
  #para obtener el número de individuos en los p posibles patrones
  num_respuestas <- ncol(y_mp) #numero de variables respuesta
  c <- 1 #acumulador de los p posibles patrones
  n_k <- numero_individuos[c] #numero de individuos por patrones
  estimaciones_matrices<- list() #aquí voy a guardar las estimaciones de 
  #los betas
  col <- ncol(y_mp) #numero de columas de la matriz de variables respuesta
  while (num_respuestas > 0) {
    matriz_pesos <- diag(pesos[1:n_k]) #aquí creo la matriz de pesos con los
    #individuos por patron
    betas <- solve(t(x_mp[1:n_k,, drop = F])%*%matriz_pesos%*%
                     x_mp[1:n_k,,drop = F])%*%t(x_mp[1:n_k,,drop = F])%*%matriz_pesos%*%y_mp[1:n_k,c:ncol(y_mp), drop = F] 
    #cálculo de los betas estimacion de la matriz de S_k por patrones
    B_k <- t(y_mp[1:n_k,c:ncol(y_mp), drop = F] - x_mp[1:n_k,, drop = F]%*%betas)%*%
      matriz_pesos%*%(y_mp[1:n_k,c:ncol(y_mp), drop = F] - x_mp[1:n_k,,drop = F]%*%betas)
    estimaciones_matrices[[c]] <- list(betas, B_k) #resultado de las matrices de los
    #betas y S_k estimados por patrón.
    num_respuestas <- num_respuestas - 1
    c <- c + 1 #corre los patrones
    n_k <- n_k + numero_individuos[c] #aqui me garantiza que a medida que va
    #aumentando los patrones va indluyendo los demas individuos a la estimacion
  }
  return(estimaciones_matrices)
}
























