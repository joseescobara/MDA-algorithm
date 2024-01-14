
#----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#funcion para la imputacion de datos que destruyen el patrón monótono

#-----------------------------------------------------------------------------
imputacion <- function(y_mp, x_mp, pesos, B_0, M_0){
  datos_obs <- matriz_patrones(y_mp)[[1]] #devuelve la lista de los datos
  #por patrones
  filas_correspo <- nrow(y_mp) #filas del total de datos
  columnas_patrones <- numero_columna(y_mp) #columnas donde empiezan los 
  #patrones
  i <- 1 #contador de comienzo de patrones y extractor de matrices de acuerdo
  #a cada patrón para luego ser analizadas
  c <- length(datos_obs) #total de patrones o matrices por patrones
  conta <- 0 #contador de las filas
  datos_imputados <- list() #datos finales
  
  while (i <= c) {
    n <- datos_obs[[i]] #datos de los individuos en el patron i
    col <- columnas_patrones[i] #columnas donde empizan los patrones
    #para ajustarle el modelo adecuado a cada individuo de acuerdo a su patron
    k <- nrow(n) #número de filas de la matriz en los patrones
    m <- ncol(n) #número de columnas de la matriz en los patrones
    #aquí mantengo las filas y corro por columnas.
    for (j in 1:k){ #corredor de individuos
      conta <- conta + 1 #este me mantiene la secuencia de las filas en orden
      #para que no se me actualicen cada vez que j comienza a recorrer las filas
      for (f in 1:m) {#corren las columnas
        if(is.na(n[j,f]) == T){
          imp <- mvrnorm(1, t(B_0[, col:ncol(y_mp), drop = F])%*%x_mp[conta:conta,], 
                        M_0[col: ncol(y_mp), col: ncol(y_mp)]/pesos[conta])
          #aquí es suponiendo que conozco mv y betas
          n[j,f] <- imp[f]
        }
      }
    }
    datos_imputados[[i]] <- n #agrega las matrices con sus datos imputados 
    #respectivamente 
    i <- i + 1
  }
  return(datos_imputados)  
}
#------------------------------------------------------------------------------
#funcion para obtener la matriz final imputada
datos_imputados <- function(y_mp, x_mp, pesos, B_0, M_0){
  patrones_m1 <- imputacion(y_mp, x_mp, pesos, B_0, M_0) #lista de las 
  #matrices imputadas por la funcion impotacion
  indiv_m1 <- matriz_patrones(y_mp)[[2]] #vector de individuos por
  #patron
  
  if(length(indiv_m1) > 1){
  y_mp1 <- y_mp
  y_mp1[1:indiv_m1[1],
            colnames(patrones_m1[[1]])]<- patrones_m1[[1]]
  #en la operación anterior, reemplaza la primera parte de la matriz por parte
  #de la matriz imputada
  acum1 <- indiv_m1[1]
  
  for(j in 2:(length(indiv_m1))){ #recorre los patrones desde dos
    y_mp1[(acum1 + 1):(acum1 + indiv_m1[(j)]),
              colnames(patrones_m1[[j]])] <- patrones_m1[[j]]
    
    acum1 <- acum1 + indiv_m1[(j)]
  }
  }
  else{y_mp1 <- y_mp}
  return(y_mp1)
  }
#-----------------------------------------------------------------------------  



