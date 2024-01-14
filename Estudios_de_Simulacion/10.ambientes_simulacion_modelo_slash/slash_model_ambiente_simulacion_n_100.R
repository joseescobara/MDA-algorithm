
matriz_diseño_100 <- function(n = 100){
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

muestra_slash_completa_100 <- function(m_dis, B_1,M_1, nu_1, N_muestras = 1000){
  lista_muestras <- list()
  set.seed(12365)
  semillas <- sample(1:4032513, N_muestras)
  conta <- 1
  for (i in 1:N_muestras){
    set.seed(semillas[i])
    
    muestra_ten <- m_dis%*%B_1 + rmslash(nrow(m_dis), 
                                         center = c(0, 0, 0), Scatter = M_1, df = nu_1)
    d_comp <- sum(dista_maha_indiv(muestra_ten, m_dis, B_1, M_1)) #asegurar 
    #atípicos
    if(d_comp >= 639 & d_comp <= 689){
      lista_muestras[[conta]] <- muestra_ten
      conta <- conta + 1
    }
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
    menores_1 <- which(d_m > 2 & d_m < 15) #aquí valido cuáles son esos individuos
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

muestra_faltantes_slash_100 <- function(muestras, m_diseño, B_1, M_1){
  N_muestras <- length(muestras) #numero de muestras a analizar
  lista_muestra_faltantes <- list()#aquí voy a guardar las muestras resultantes
  filas_1_analizar <- detecta_atipicos_0(muestras, m_diseño, B_1, M_1)[[1]]
  muestras_1_analizar <- detecta_atipicos_0(muestras, m_diseño, B_1, M_1)[[2]]
  n <- nrow(m_diseño)
  conta <- 1
  if(length(muestras_1_analizar) !=0 ){
    set.seed(78013)
    corredor <- length(muestras_1_analizar)
    semillas <- sample(1:9821461, corredor)
    for (i in 1:corredor){
      set.seed(semillas[i])
      muestras_insertar <- muestras[[muestras_1_analizar[[i]]]] #muestra ith a
      #analizar
      filas <- filas_1_analizar[[i]] #filas de la muestra ith a analizar
      
      #patron 1 c(1,0,0)
      p <- runif(length(filas)) #probabilidad para extraer cada fila
      muestra_filas1 <- sample(filas, ceiling(n*(62/173)), prob = p, replace = F)
      muestras_insertar[muestra_filas1, 1] <- NA #incluyo patron 1
      #incluyo patron 2 (0,1,1)
      filas2 <- setdiff(filas, muestra_filas1) 
      p1 <- runif(length(filas2))
      muestra_filas2 <- sample(filas2, floor(n*(6/173)), prob = p1, replace = F)
      muestras_insertar[muestra_filas2, c(2,3)] <- NA 
      #para incluir el patron 3 (1,1,0)
      filas3 <- setdiff(filas2, muestra_filas2) 
      p2 <- runif(length(filas3))
      muestra_filas3 <- sample(filas3, ceiling(n*(1/173)), prob = p2, replace = F)
      muestras_insertar[muestra_filas3, c(1,2)] <- NA
      #inclusión patron 4 #(0,1,0)
      filas4 <- setdiff(filas3, muestra_filas3)
      p3 <- runif(length(filas4))
      muestra_filas4 <- sample(filas4, ceiling(n*(1/173)), prob = p3, replace = F)
      muestras_insertar[muestra_filas4, 2] <- NA
      #---------
      prop_faltantes <- apply(is.na(muestras_insertar), 2, sum)
      d_comp <- sum(dista_maha_indiv(muestras_insertar, m_diseño, B_1, M_1))
      if(prop_faltantes[1] >= prop_faltantes[2] & prop_faltantes[2] >= 
         prop_faltantes[3] & d_comp >=  469 & d_comp <= 689){
        lista_muestra_faltantes[[conta]] <- muestras_insertar
        conta <- conta + 1
      }
      
    }
  }
  return(lista_muestra_faltantes)
}
