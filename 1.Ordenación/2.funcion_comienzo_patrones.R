#funcion que me devuelve el comienzo de cada patrón

#-----------------------------------------------------------------------
#funcion begin_patrones
begin_patrones <- function(y_mp){
  matriz_auxiliar <- t(1*apply(y_mp, 1, is.na)) #creo una matriz indicadora 
  #de datos faltantes
  faltantes_filas <- conteos_filas(matriz_auxiliar)[[2]] #cuento los faltantes
  #por fila hasta que encuentre uno diferente a un faltante
  comienzo_patrones <- c(1) #el primer patrón siempre comienza en la fila 1
  j <- 2 #dado que en el índice 1 es igual a 1, este comienza en el segundo 
  #índice para el siguiente patrón
  for (i in 2:length(faltantes_filas)){#i comienza en 2
    if(faltantes_filas[i] != faltantes_filas[i-1]){
      comienzo_patrones[j] <- i
      j <- j + 1
    }
  }
  return(comienzo_patrones)
}
#---------------------------------------------------------------------------





















