#-funcion para ordenar la matriz de covariables respecto a la nueva ordenacion 
#de la matriz y en patron monotono-------------------------------------

ordena_covariables <- function(y_df, x_df){
  x_df <- as.matrix(x_df) #cambio el formato de la matriz de covariables 
  #a formato matrix
  x <- as.matrix(y_df)#convierte la matriz de datos en matriz
  n <- nrow(x) #filas de la matriz de datos
  p <- ncol(x) #columnas de la matriz de datos
  storage.mode(x) <- "double" #los convierte en datos reales
  r <- 1 * is.na(x) #matriz de indicadores de datos faltantes, marcando
  #1 los que son faltantes y 0 los que están observados
  nmis <- as.integer(apply(r, 2, sum)) #suma de datos faltantes por columnas
  nmis_fila <- as.integer(apply(r, 1, sum)) #una alternativa para visualizar
  #algunas cosas, para ordenar por las filas que menos tienen hasta las que 
  #más tienen
  k_1 <- y_df[order(nmis_fila, decreasing = F), ] #ordena las filas de
  #acuerdo a la proporción de datos faltantes en orden creciente
  ord_fil <- c(1:nrow(x_df))#filas de la matriz de covariables 
  x_df <- cbind(x_df,ord_fil) #crea una nueva matriz de covariables, 
  #agregandole la columna de filas
  x_df1 <- x_df[order(nmis_fila, decreasing = F), ] #ordena la matriz de
  #covariables respecto al orden de las filas de la matriz de variablles
  #respuesta, desde la que menos datos faltantes tiene a la que más datos 
  #faltantes tiene
  k_2 <- k_1[,order(nmis, decreasing = T)] #ordena las columnas de la 
  #matriz de variables respuestas por la proporción de datos faltantes 
  #en orden decreciente
  faltantes <- is.na(k_2[, 1]) #toma la primera columna de de la matriz
  #de variables respuesta ordenada de acuerdo a la proporción de datos 
  #faltantes por culumnas en orden decreciente y obtenemos el vector 
  #lógico de True para datos faltantes y False para datos no faltantes 
  k_3 <- k_2[order(faltantes),] #es la nueva matriz ordenada de acuerdo
  #a la primera columna, dado que, todas las filas que tengan datos
  #faltantes en la primera entrada quedaran debajo, mientras que las otras
  #quedaran arriba
  x_df2 <- x_df1[order(faltantes),] #segunda ordenación de la matriz
  #de covariables de acuerdo a la nueva ordenación de la nueva matriz
  #de variables respuesta
  matriz_auxiliar <- 1*is.na(k_3) #creo una matriz auxiliar de indicadores de 
  #datos faltantes para ser utilizada posteriormente para el detectar cambio
  #de patron y encontrar el patrón monótono
  return(list(k_3, matriz_auxiliar, x_df2))
}


#---------------funcion para organizar el patrón monótono filas------------

conteos_filas <- function(matriz_auxiliar){
  n_filas <- nrow(matriz_auxiliar)#número de filas de la matriz indicadora
  #de datos faltantes
  filas <- seq(1, n_filas) #secuencia de 1 a el total del número de filas
  n_col <- ncol(matriz_auxiliar)#columnas de la matriz auxiliar
  p <- 1 #acumulador de filas
  faltantes_filas <- numeric() #numero de faltantes en cada fila
  z <- 0 #contador de faltantes por filas en cada iteración
    for (i in 1:n_filas) {#iterador de filas (filas fijas)
      for (j in 1:n_col) {#iterador de columnas (recorre todas las columnas)
        if( matriz_auxiliar[i, j] == 1){
          z <- z + 1
        }
        else{break}
      }
      faltantes_filas[p] <- z #agregador del número de faltantes por fila
      p <- p + 1 #acumulador de filas
      z <- 0 #actualización del contador de faltantes para una nueva iteración
    }
  return(list(filas, faltantes_filas))  
}

#-----------------funcion para obtener el patrón monótono------------------

patron_monotono <- function(y_df, x_df){
  matriz_auxiliar <- ordena_covariables(y_df, x_df)[[2]] #hace uso de la
  #función ordena covariables para extraer la matriz indicadora de datos faltantes
  ordenadores <- conteos_filas(matriz_auxiliar)[[2]] #aplica la función conteos fila
  orden_faltantes <- order(ordenadores) #orden de las filas respecto al
  #número de datos faltantes en orden creciente
  y_df1 <- ordena_covariables(y_df, x_df)[[1]] #hago uso de la función
  #ordena covariables
  x_df1 <- ordena_covariables(y_df, x_df)[[3]] #hago uso de la función 
  #ordena covariables
  y_mp <- y_df1[orden_faltantes,  ] #variables respuestas ordenadas de acuerdo
  #a la proporción de datos faltantes consecutivos por fila extraida de la 
  #función conteos fila.
  x_mp <- x_df1[orden_faltantes, ] #variables respuestas ordenadas de acuerdo
  #a la proporción de datos faltantes consecutivos por fila extraida de la 
  #función conteos fila.
  return(list(y_mp, x_mp[,1:(ncol(x_mp)-1)], x_mp[,ncol(x_mp)]))
}

