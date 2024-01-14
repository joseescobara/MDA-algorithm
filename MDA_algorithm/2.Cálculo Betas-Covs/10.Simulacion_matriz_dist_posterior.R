#funciones para simular la distribucion posterior de la matriz
#de dispersión
#-----------------------------------------------------------------------------
#funcion que me extrae las matrices B_k para k = 1, ..., p
matrices_Bk <- function(y_mp, x_mp, pesos){
  Sks <- estimacion_betas_Bk(y_mp, x_mp, pesos) #hago uso de la función 
  #estimacion_betas_bk para extraer las matrices de sumas cuadradas 
  #y productos cruzados 
  Sk <- list()
  for (i in 1:length(Sks)) {
    Sk[[i]] <- Sks[[i]][[2]]
  }
  return(Sk)
}
#-----------------------------------------------------------------------------

#------------------------------------------------------------------------------

#función que me devuelve las matrices b_k
matrices_bk_inv <- function(y_mp, x_mp, pesos){
  #uso de la función para la estimación de las matrices de ceoficientes
  #y las matrices de productos cruzados y suma totaltes de los errores al
  #cuadrado ponderado
  matrices_BK <- matrices_Bk(y_mp, x_mp,pesos)#utilizo la función matrices_Bk
  #para extraer las matrices b_k, suponiendo que A_k = 0
  bk_inv <- lapply(matrices_BK, solve) # aplico la inversa a cada elemento
  #de las matrices_BK
  return(bk_inv)
}
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#funcion que me devuelve las matrices L_k para k = 1, ..., p
matrices_lk <- function(y_mp, x_mp, pesos){
  b_k_inv <- matrices_bk_inv(y_mp, x_mp, pesos) #hago uso de la función 
  #matrices b_k_inversas para aplicarle la factorización de choleski a cada
  #matrix
  L_k <- lapply(lapply(b_k_inv, chol), t) #aplico choleski a cada elemento de las 
  #B_K_inv y luego la trasnpuesta
  return(L_k)
}
#---------------------------------------------------------------------------

#----------------------------------------------------------------------------
#funcion que me simula los vectores u_k

vectores_uk <- function(y_mp, x_mp){
  p <- ncol(y_mp) #número de variables respuesta
  q <- ncol(x_mp) #número de covariables
  m <- p #Supuesto de la jeffreys a priori, m = p
  patron <- 1 #recorre los p posibles patrones
  nj_individuos <- individuos_de_acuerdo_patron(y_mp) #individuos de 
  #acuerdo a los p posibles patrones
  matriz_u <- matrix(0, nrow = p, ncol = p) #matriz para rellenar de acuerdo 
  #a las condiciones que debe cumplir los elementos de cada u_k
  actu_ind <- 0 #contador de individuos
  df <- 0 # grados de libertad
  
  while (patron <= p){ #aquí patrón va aumentando en 1 para recorrer los p 
    #posibles patrones
    actu_ind <- actu_ind + nj_individuos[patron] #suma de individuos a medida
    #que va aumentando el patrón para n_1 + n_2 + ... + n_k
    df <- actu_ind - patron + m - p - q + 1 #grados de libertad de la chi
    if(actu_ind > patron + (p-m-1+q)){#condicion para que la distribución 
      #de sigma_inv sea HH'
      for (fila in patron:p){ #mantengo columna fija, corren las filas.
        if(fila == patron){ #generadora de la chi i=j
          matriz_u[patron, patron] <- sqrt(rchisq(1, df)) #elmento u_jj
        }
        else #generadoras de la N(0,1) 1<=j<i<=p
        { matriz_u[fila, patron] <- rnorm(1, 0, 1) #en caso de que i != patron
        }
      }
    }
    else{
      return("No se cumple la condición específicada para que la distribución
             de sigma invesra sea igual a la distribución de HH' ")
      break
    }
    patron <- patron + 1 #aumentando los p posibles patrones
  }
  return(matriz_u)
} 
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#funcion para la simulacion de H
matriz_H <- function(y_mp, x_mp, pesos){
  u_matriz <- vectores_uk(y_mp, x_mp) #genero la matriz de los u_k's
  L <- matrices_lk(y_mp, x_mp,pesos) #genero las matrices L_k
  p <- ncol(y_mp) #numero de variables respuesta
  H <- matrix(0, nrow = p, ncol =p) #genero la matriz nula pxp, para luego ser 
  #rellenada
  for (i in 1:p){
    H[i:p, i] <- L[[i]]%*%u_matriz[i:p, i, drop = F]#relleno la parte triangular inferior
    #con L_K%*%u_k, k = 1, ..., p
  }
  return(H)
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#funcion para la simulacion de la matriz de dispersión
matriz_sigma <- function(y_mp, x_mp, pesos){
  H <- matriz_H(y_mp, x_mp,pesos) #hago uso da la función matriz_h para la 
  #respectiva simulación de H.
  mvcov <- solve(H%*%t(H)) #simulo una muestra de sigma
  colnames(mvcov) <- colnames(y_mp) #nombro las columnas de sigma
  rownames(mvcov) <- colnames(y_mp) #nombre las filas de sigma
  return(list(mvcov, H))
}
#-----------------------------------------------------------------------------













