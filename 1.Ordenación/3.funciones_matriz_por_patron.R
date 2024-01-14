
#------------------------------------------------------------------------------
#en este caso modifiqué la función a que solo me reciba la matriz en patrón 
#monótono
find_pattern <- function(y_mp){
  begin <- begin_patrones(y_mp) #hago uso de la función begin_patrones
  #para extraer las filas donde comienzan los patrones
  h <- list() #lista donde guardo las filas donde comienza cada patrón
  for(i in 1:length(begin)){
    h[[i]] <- y_mp[begin[i]:begin[i], , drop = F] #extraigo las filas
    #donde empieza cada patron
  }
  return(h)
}
#-----------------------------------------------------------------------------
numero_columna <- function(y_mp){
  filas_comienzo <- find_pattern(y_mp) #utiliza la funcion find_pattern
  #para utilizar las filas donde comienzan los patrones
  columnas <- numeric() #aquí guardaré las columnas de los patrones
  for(i in 1:length(filas_comienzo)){ #recorre los elementos de la lista
      k <- filas_comienzo[[i]] #extraigo cada elemento de la lista obtenida
      #por medio de la función find_pattern para analizarlos.
      columnas[i] <- which(is.na(k)==F)[1] #aplico la función which para
      #que me extraiga los elementos de ese vector que cumplen la condicón
    }
  return(columnas)
}
#-----------------------------------------------------------------------------
#le hice el ajuste de que solo recibiera los datos en patrón monótono
matriz_patrones <- function(y_mp){
  datos_mp <- data.frame(y_mp, check.names = F) #convierto los datos en un data frame
  begin <- begin_patrones(datos_mp) #usa la funcion find_patter para detectar
  #las filas donde comienza cada patrón
  m_patrones <- list() #lista que va a contoner los individuos desagregados por
  #patron
  patrones <- numero_columna(datos_mp) #hace uso de la función numero_columna
  #para obtener las columnas donde comienza cada patrón
  numero_individuos_por_patron <- numeric() #numero de individuos por patrón
  t <- length(begin) #numero de patrones reales en la matriz de datos
  r <- 0 #acumulador de filas que se van utilizando
  while (t > 0) {
    #extraigo las matrices por patrón de atrás hacia adelante
    K <- datos_mp[begin[t]:(nrow(y_mp)-r),patrones[t]:ncol(y_mp), drop = F]
    K <- as.matrix(K) #convierto las matrices por patrónes en formato de 
    #matriz
    num_filas <- nrow(K) #filas que voy utilizando
    numero_individuos_por_patron[t] <- num_filas #cuento el número el número 
    #de individuos por patrón a medida que voy sacando las amtrices por patrón
    m_patrones[[t]] <- K #voy guardando las matrices extraidas de acorde al 
    #patrón
    t <- t - 1 #va disminuyendo para sacar las filas donde comienzan los 
    #patrones de atrás hacia adelante
    r <- r + num_filas #va acumulando las filas que voy utilizando
  }
  salida <- list(m_patrones, numero_individuos_por_patron)
  return(salida)
}
#-----------------------------------------------------------------------------













