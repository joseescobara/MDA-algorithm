matriz_diseño <- function(n = 50){
  set.seed(1)#semilla
  edad <- rgamma(n, 2.16585, rate = 1.308587)#simulo la variable edad 
  #con una gamma
  set.seed(2)#semilla
  sexo <- rbinom(n, 1,  0.6127168)#simulamos la variable sexo con una
  #distribucion bernoulli
  set.seed(3)#semilla
  tiempo_lechem <- rgamma(n, 1.814255, rate = 0.1975243)#simulo la variable
  #leche materna con una gamma
  x_0 <- rep(1, times = n)
  Matriz_Diseño <- cbind(x_0, edad, sexo, tiempo_lechem)
  return(Matriz_Diseño)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

muestra_slash_completa <- function(m_dis, B_1, M_1, nu_1, N_muestras = 1000){
  lista_muestras <- list()
  set.seed(12365)
  semillas <- sample(1:132513, N_muestras)
  for (i in 1:N_muestras) {
    set.seed(semillas[i])
    lista_muestras[[i]] <- m_dis%*%B_1 + rmslash(nrow(m_dis), 
                            center = c(0, 0, 0), Scatter = M_1, df = nu_1)
  }
  
  return(lista_muestras)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#vamos a crear una función para detectar los datos no atípicos 


detecta_atipicos_0 <- function(muestras, m_diseño, B_1, M_1){
  lista_filas <- list() #aquí creo una lista donde voy a guardar las filas
  #de las muestras simuladas que están cercanas a su media
  muestras_requisito <- c() #aquí guardo las filas con esos requisitos
  conta <- 1
  for (i in 1:length(muestras)){
    d_m <- dista_maha_indiv(muestras[[i]], m_diseño, B_1, M_1) #le calculo las 
    #distancias de mahalanobis a cada muestra para detectar cuáles valores está 
    #cercanos a su media. 
    menores_1 <- which(d_m < 1) #aquí valido cuáles son esos individuos
    if(length(menores_1) !=0){
    lista_filas[[conta]] <- menores_1
    muestras_requisito[conta] <- i #aquí guardo las muestras que cumplen ese 
    #requisito
    conta <- conta + 1
      }
  }
  resultado <- list(lista_filas, muestras_requisito)
  return(resultado)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
detecta_atipicos_1 <- function(muestras, m_diseño, B_1, M_1){
  lista_filas <- list() #aquí creo una lista donde voy a guardar las filas
  #de las muestras simuladas que están cercanas a su media
  muestras_requisito <- c() #aquí guardo las filas con esos requisitos
  conta <- 1
  for (i in 1:length(muestras)){
    d_m <- dista_maha_indiv(muestras[[i]], m_diseño, B_1, M_1) #le calculo las 
    #distancias de mahalanobis a cada muestra para detectar cuáles valores está 
    #cercanos a su media. 
    menores_1 <- which(d_m <= 2) #aquí valido cuáles son esos individuos
    if(length(menores_1) !=0){
      lista_filas[[conta]] <- menores_1
      muestras_requisito[conta] <- i #aquí guardo las muestras que cumplen ese 
      #requisito
      conta <- conta + 1
    }
  }
  resultado <- list(lista_filas, muestras_requisito)
  return(resultado)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

muestra_faltantes_slash <- function(muestras, m_diseño, B_1, M_1){
  N_muestras <- length(muestras) #numero de muestras a analizar
  lista_muestra_faltantes <- list()#aquí voy a guardar las muestras resultantes
  filas_1_analizar <- detecta_atipicos_0(muestras, m_diseño, B_1, M_1)[[1]]
  muestras_1_analizar <- detecta_atipicos_0(muestras, m_diseño, B_1, M_1)[[2]]
  
  muestras_faltantes <- setdiff(N_muestras, muestras_1_analizar)#analizar las muestras
  #faltantes (los índices)
  n <- nrow(m_diseño)
  conta <- 1
  if(length(muestras_1_analizar) !=0 ){
    set.seed(78013)
    corredor <- length(muestras_1_analizar)
    semillas <- sample(1:9821461, corredor)
    for (i in 1:corredor){
      set.seed(semillas[i])
      filas <- filas_1_analizar[[i]]
      p <- runif(length(filas))
      muestra_filas1 <- sample(filas, n*(63/173), prob = p, replace = T)
      muestra_filas2 <- sample(filas, n*(8/173), prob = p, replace = T)
      muestra_filas3 <- sample(filas, n*(6/173), prob = p, replace = T)
      muestras_insertar <- muestras[[muestras_1_analizar[[i]]]]
      muestras_insertar[muestra_filas1, 1] <- NA
      muestras_insertar[muestra_filas2, 2] <- NA
      muestras_insertar[muestra_filas3, 3] <- NA
    
      suma_NA <- apply(is.na(muestras_insertar), 1, sum)
      detector <- which(suma_NA > 2)
      if(length(detector) == 0){
        lista_muestra_faltantes[[conta]] <- muestras_insertar
        conta <- conta + 1
      }
      
      }
  }
  
if(length(muestras_faltantes) != 0){
  muestras_2 <- muestras[[muestras_faltantes]]
  filas_restantes <- detecta_atipicos_1(muestras_2, m_diseño, B_1, M_1)[[1]]
  #filas que cumplen con la condicion
  muestras_2_analizar <- detecta_atipicos_1(muestras_2, m_diseño, B_1, M_1)[[2]]
  #inidices de muestras que voy a analizar
  set.seed(7180113)
  corredor1 <- length(muestras_2)
  semillas1 <- sample(1:98211631, corredor1)
  for (i in 1:corredor1){
    set.seed(semillas1[i])
    filas <- filas_restantes[[i]]
    p <- runif(length(filas))
    muestra1_filas1 <- sample(filas, n*(63/173), prob = p, replace = T)
    muestra1_filas2 <- sample(filas, n*(8/173), prob = p, replace = T)
    muestra1_filas3 <- sample(filas, n*(6/173), prob = p, replace = T)
    muestras_insertar <- muestras[[muestras_2[[i]]]]
    muestras_insertar[muestra1_filas1, 1] <- NA
    muestras_insertar[muestra1_filas2, 2] <- NA
    muestras_insertar[muestra1_filas3, 3] <- NA
    suma_NA1 <- apply(is.na(muestras_insertar), 1, sum)
    detector1 <- which(suma_NA1 > 2)
    prop_faltantes <- apply(is.na(muestras_insertar), 2, sum)#para asegurar que 
    #la segunda columna siempre tiene mayor numero de datos faltantes
    if(length(detector1) == 0 & prop_faltantes[2] >= prop_faltantes[3]){
      lista_muestra_faltantes[[conta]] <- muestras_insertar
      conta <- conta + 1
    }
  }
}
return(lista_muestra_faltantes)
}









