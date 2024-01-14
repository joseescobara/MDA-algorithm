#-------------------------------------------------------------------------
#funcion que me genera las matrices (X´pesosX)^-1 , k = 1, ..., p.
#aquí en esta función la información que contiene y_mp es irrelevante. 
#Sólo es utilizada para extraer los individuos en los p posibles patrones.

generacion <- function(y_mp, x_mp,pesos){
  numero_individuos <- individuos_de_acuerdo_patron(y_mp)#individuos
  #por patrones reales
  c <- 1 #para extraer los individuos por patron, correr en el 
  #ciclo while y agregar las M matrices.
  n_k <- numero_individuos[c] #numero de individuos por patrones
  matrices <- list()
  col1 <- ncol(y_mp) #numero de variables respuesta
  while (c <= col1) {
    matriz_pesos <- diag(pesos[1:n_k])
    M <- solve(t(x_mp[1:n_k,, drop = F])%*%matriz_pesos%*%x_mp[1:n_k,, drop = F]) 
    #matriz (X´pesosX)^-1 #k = 1, ..., p. A medida que van aumentando los patrones se va agregando
    #más información.
    matrices[[c]] <- M 
    c <- c + 1 #recorre de k = 1, ..., p
    n_k <- n_k + numero_individuos[c]
  }
  return(matrices)
}
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#funcion que me genera las g_k, aquí en esta función la información que
#contiene y_mp es irrelevante. Sólo es utilizada para extraer los individuos
#en los p posibles patrones.
matrices_gk <- function(y_mp, x_mp, pesos){
  matrices_generacion <- generacion(y_mp, x_mp, pesos) #matrices 
  #(X´pesosX)^-1 , k = 1, ..., p
  g_k_matrices <- lapply(lapply(matrices_generacion, chol), t)#aplico la función chol a cada
  #uno de los elementos de la lista y luego la trasnouesta
  #g_k_matrices <- 
  return(g_k_matrices)
}
#--------------------------------------------------------------------------

#---------------------------------------------------------------------------
#función para extraer los betas estimados por patrón a partir de los datos
betas_estimados_patron <- function(y_mp, x_mp,pesos){
  r <- estimacion_betas_Bk(y_mp, x_mp, pesos) #calcula una vez
  betas <- list()
  for (i in 1:length(r)) {
    betas[[i]] <- r[[i]][[1]]
  }
  return(betas)
}

#--------------------------------------------------------------------------
#funcion para simular la primera parte de los betas
primera_parte_betas <- function(y_mp, x_mp,pesos,H){
  q <- ncol(x_mp) #número de covariables
  p <- ncol(y_mp) #número de variables respuesta
  g_k <- matrices_gk(y_mp, x_mp, pesos) #matrices g_k
  Z <- matrix(rnorm(q*p), byrow = T, nrow = q, ncol = p) #matriz aleatoria 
  #cuyos elementos son distribuidos normal estandár. 
  matriz_primera_parte <- matrix(0, nrow =q, ncol = p) #primera parte de 
  #la simulacion de los betas
  for (i in 1:p){#porque la matriz tiene p columnas
    matriz_primera_parte[, i] <- g_k[[i]]%*%Z[, i, drop = F]
  }
  matriz_primera_parte <- matriz_primera_parte%*%solve(H)
  return(matriz_primera_parte)
}

#----------------------------------------------------------------------------
#funcion que me genera la segunda parte de los betas
segunda_parte_betas <- function(y_mp, x_mp,pesos, H){
  betas <- betas_estimados_patron(y_mp, x_mp, pesos) #betas en la interacion i
  #cuando los datos tienen patrón monótono
  q <- ncol(x_mp) #numero de covariables
  p <- ncol(y_mp) #numero de variables respuestas
  matriz_segunda_parte <- matrix(0, nrow = q, ncol =p)
  for (i in 1:p){ #porque la matriz tiene p columnas
    matriz_segunda_parte[, i] <- betas[[i]]%*%H[i:p, i, drop = F]
  }
  matriz_segunda_parte <- matriz_segunda_parte%*%solve(H)
  return(matriz_segunda_parte)
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#funcion de generacion de los betas condicional a sigma
generacion_betas <- function(y_mp, x_mp, pesos,H){
  parte_1 <- primera_parte_betas(y_mp, x_mp,pesos,H)
  parte_2 <- segunda_parte_betas(y_mp, x_mp,pesos,H)
  betas_simulados <- parte_1 + parte_2
  colnames(betas_simulados) <- colnames(y_mp)
  rownames(betas_simulados) <- colnames(x_mp)
  return(betas_simulados)
}
#----------------------------------------------------------------------------



