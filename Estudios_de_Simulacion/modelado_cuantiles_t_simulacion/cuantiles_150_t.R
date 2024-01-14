#modelamiento de cuantiles muestra de tamaño 150 t
setwd("ruta")
load("estimaciones_150_modelo_t.Rdata")
load("muestra_tamaño_150_t_completa.Rdata")
X_dis <- matriz_diseño_150_t
library(stats)

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------


#-----------------------calculo_cuantiles_slash_25-----------------------

calculo_cuantiles_niñas_25 <- function(lista_estimaciones, X_dis){
  mean.edad <- mean(X_dis[, 2])
  mean.tlm <- mean(X_dis[, 4])
  
  cuantiles_per_braq <- c()
  cuantiles_peso <- c()
  cuantiles_talla <- c()
  
  for (i in 1:length(lista_estimaciones)){
    coef.t <- lista_estimaciones[[i]][[1]]
    matdisp.t <- lista_estimaciones[[i]][[2]]
    nu.t <- lista_estimaciones[[i]][[3]]
    q.t_0.25 <- qt(p = 0.25, df = nu.t)
    
    sqrt.sigma11 <- sqrt(matdisp.t[1,1])
    sqrt.sigma22 <- sqrt(matdisp.t[2,2])
    sqrt.sigma33 <- sqrt(matdisp.t[3,3])
    
    cuantiles_per_braq[i]<- exp(coef.t[1,1] + coef.t[2,1]* mean.edad
                                + coef.t[4,1] * mean.tlm + 
                                  sqrt.sigma11 * q.t_0.25)
    
    cuantiles_peso[i]<- exp(coef.t[1,2] + coef.t[2,2]* mean.edad
                            + coef.t[4,2] * mean.tlm + 
                              sqrt.sigma22 * q.t_0.25)
    
    cuantiles_talla[i]<- exp(coef.t[1,3] + coef.t[2,3]* mean.edad
                             + coef.t[4,3] * mean.tlm + 
                               sqrt.sigma33 * q.t_0.25)
    
  }
  
  return(list(cuantiles_per_braq, cuantiles_peso, cuantiles_talla))
  
}


cuantil_25_niña <- calculo_cuantiles_niñas_25(lista_estimaciones_150_t, matriz_diseño_150_t)
median(cuantil_25_niña[[1]])
median(cuantil_25_niña[[2]])
median(cuantil_25_niña[[3]])


m1 <- median(cuantil_25_niña[[1]])
mades1 <- abs(cuantil_25_niña[[1]] - m1)
round(median(mades1), 3)
round(m1, 3)

m2 <- median(cuantil_25_niña[[2]])
mades2 <- abs(cuantil_25_niña[[2]] - m2)
round(median(mades2), 3)
round(m2, 3)


m3 <- median(cuantil_25_niña[[3]])
mades3 <- abs(cuantil_25_niña[[3]] - m3)
round(median(mades3), 3)
m3

#------------------------------------------------------------------------------

#--------------------calculo_cuantiles_slash_50---------------------------------

calculo_cuantiles_niñas_50 <- function(lista_estimaciones, X_dis){
  mean.edad <- mean(X_dis[, 2])
  mean.tlm <- mean(X_dis[, 4])
  
  cuantiles_per_braq <- c()
  cuantiles_peso <- c()
  cuantiles_talla <- c()
  
  for (i in 1:length(lista_estimaciones)){
    coef.t <- lista_estimaciones[[i]][[1]]
    matdisp.t <- lista_estimaciones[[i]][[2]]
    nu.t <- lista_estimaciones[[i]][[3]]
    q.t_0.50 <- 0
    
    sqrt.sigma11 <- sqrt(matdisp.t[1,1])
    sqrt.sigma22 <- sqrt(matdisp.t[2,2])
    sqrt.sigma33 <- sqrt(matdisp.t[3,3])
    
    cuantiles_per_braq[i]<- exp(coef.t[1,1] + coef.t[2,1]* mean.edad
                                + coef.t[4,1] * mean.tlm + 
                                  sqrt.sigma11 * q.t_0.50)
    
    cuantiles_peso[i]<- exp(coef.t[1,2] + coef.t[2,2]* mean.edad
                            + coef.t[4,2] * mean.tlm + 
                              sqrt.sigma22 * q.t_0.50)
    
    cuantiles_talla[i]<- exp(coef.t[1,3] + coef.t[2,3]* mean.edad
                             + coef.t[4,3] * mean.tlm + 
                               sqrt.sigma33 * q.t_0.50)
    
  }
  return(list(cuantiles_per_braq, cuantiles_peso, cuantiles_talla))
}


cuantil_50_niñas <- calculo_cuantiles_niñas_50(lista_estimaciones_150_t, matriz_diseño_150_t)

median(cuantil_50_niñas[[1]])
median(cuantil_50_niñas[[2]])
median(cuantil_50_niñas[[3]])



m4 <- median(cuantil_50_niñas[[1]])
mades4 <- abs(cuantil_50_niñas[[1]] - m4)
round(median(mades4), 4)
m4

m5 <- median(cuantil_50_niñas[[2]])
mades5 <- abs(cuantil_50_niñas[[2]] - m5)
round(median(mades5), 4)
m5


m6 <- median(cuantil_50_niñas[[3]])
mades6 <- abs(cuantil_50_niñas[[3]] - m6)
round(median(mades6), 4)
m6

#-------------------------------------------------------------------------------


#--------------------calculo_cuantiles_slash_75---------------------------------

calculo_cuantiles_niñas_75 <- function(lista_estimaciones, X_dis){
  mean.edad <- mean(X_dis[, 2])
  mean.tlm <- mean(X_dis[, 4])
  
  cuantiles_per_braq <- c()
  cuantiles_peso <- c()
  cuantiles_talla <- c()
  
  for (i in 1:length(lista_estimaciones)){
    coef.t <- lista_estimaciones[[i]][[1]]
    matdisp.t <- lista_estimaciones[[i]][[2]]
    nu.t <- lista_estimaciones[[i]][[3]]
    q.t_0.75 <- qt(p = 0.75, df = nu.t)
    
    sqrt.sigma11 <- sqrt(matdisp.t[1,1])
    sqrt.sigma22 <- sqrt(matdisp.t[2,2])
    sqrt.sigma33 <- sqrt(matdisp.t[3,3])
    
    cuantiles_per_braq[i]<- exp(coef.t[1,1] + coef.t[2,1]* mean.edad
                                + coef.t[4,1] * mean.tlm + 
                                  sqrt.sigma11 * q.t_0.75)
    
    cuantiles_peso[i]<- exp(coef.t[1,2] + coef.t[2,2]* mean.edad
                            + coef.t[4,2] * mean.tlm + 
                              sqrt.sigma22 * q.t_0.75)
    
    cuantiles_talla[i]<- exp(coef.t[1,3] + coef.t[2,3]* mean.edad
                             + coef.t[4,3] * mean.tlm + 
                               sqrt.sigma33 * q.t_0.75)
    
  }
  return(list(cuantiles_per_braq, cuantiles_peso, cuantiles_talla))
}


cuantil_75_niñas <- calculo_cuantiles_niñas_75(lista_estimaciones_150_t, 
                                               matriz_diseño_150_t)

median(cuantil_75_niñas[[1]])
median(cuantil_75_niñas[[2]])
median(cuantil_75_niñas[[3]])

m7 <- median(cuantil_75_niñas[[1]])
mades7 <- abs(cuantil_75_niñas[[1]] - m7)
round(median(mades7), 3)
round(m7, 3)

m8 <- median(cuantil_75_niñas[[2]])
mades8 <- abs(cuantil_75_niñas[[2]] - m8)
round(median(mades8), 3)
round(m8,3)


m9 <- median(cuantil_75_niñas[[3]])
mades9 <- abs(cuantil_75_niñas[[3]] - m9)
round(median(mades9), 4)
m9


#------------------------------------------------------------------------------
#-----------------------------calculo_cuantiles_slash_25_niños------------------------------

calculo_cuantiles_niños_25 <- function(lista_estimaciones, X_dis){
  mean.edad <- mean(X_dis[, 2])
  mean.tlm <- mean(X_dis[, 4])
  
  cuantiles_per_braq <- c()
  cuantiles_peso <- c()
  cuantiles_talla <- c()
  
  for (i in 1:length(lista_estimaciones)){
    coef.t <- lista_estimaciones[[i]][[1]]
    matdisp.t <- lista_estimaciones[[i]][[2]]
    nu.t <- lista_estimaciones[[i]][[3]]
    q.t_0.25 <- qt(p = 0.25, df = nu.t)
    
    sqrt.sigma11 <- sqrt(matdisp.t[1,1])
    sqrt.sigma22 <- sqrt(matdisp.t[2,2])
    sqrt.sigma33 <- sqrt(matdisp.t[3,3])
    
    cuantiles_per_braq[i]<- exp(coef.t[1,1] + coef.t[2,1]* mean.edad +
                                  coef.t[3,1] + coef.t[4,1] * mean.tlm + 
                                  sqrt.sigma11 * q.t_0.25)
    
    cuantiles_peso[i]<- exp(coef.t[1,2] + coef.t[2,2]* mean.edad + 
                              coef.t[3,2] + coef.t[4,2] * mean.tlm + 
                              sqrt.sigma22 * q.t_0.25)
    
    cuantiles_talla[i]<- exp(coef.t[1,3] + coef.t[2,3]* mean.edad + 
                               coef.t[3,3] + coef.t[4,3] * mean.tlm + 
                               sqrt.sigma33 * q.t_0.25)
    
  }
  
  return(list(cuantiles_per_braq, cuantiles_peso, cuantiles_talla))
  
}


cuantil_25_niños <- calculo_cuantiles_niños_25(lista_estimaciones_150_t, 
                                               matriz_diseño_150_t)
median(cuantil_25_niños[[1]])
median(cuantil_25_niños[[2]])
median(cuantil_25_niños[[3]])


m1 <- median(cuantil_25_niños[[1]])
mades1 <- abs(cuantil_25_niños[[1]] - m1)
round(median(mades1), 4)
m1

m2 <- median(cuantil_25_niños[[2]])
mades2 <- abs(cuantil_25_niños[[2]] - m2)
round(median(mades2), 4)
m2


m3 <- median(cuantil_25_niños[[3]])
mades3 <- abs(cuantil_25_niños[[3]] - m3)
round(median(mades3), 4)
m3

#------------------------------------------------------------------------------

#--------------------calculo_cuantiles_slash_50_niños---------------------------------

calculo_cuantiles_niños_50 <- function(lista_estimaciones, X_dis){
  mean.edad <- mean(X_dis[, 2])
  mean.tlm <- mean(X_dis[, 4])
  
  cuantiles_per_braq <- c()
  cuantiles_peso <- c()
  cuantiles_talla <- c()
  
  for (i in 1:length(lista_estimaciones)){
    coef.t <- lista_estimaciones[[i]][[1]]
    matdisp.t <- lista_estimaciones[[i]][[2]]
    nu.t <- lista_estimaciones[[i]][[3]]
    q.t_0.50 <- 0
    
    sqrt.sigma11 <- sqrt(matdisp.t[1,1])
    sqrt.sigma22 <- sqrt(matdisp.t[2,2])
    sqrt.sigma33 <- sqrt(matdisp.t[3,3])
    
    cuantiles_per_braq[i]<- exp(coef.t[1,1] + coef.t[2,1]* mean.edad
                                +coef.t[3,1] + coef.t[4,1] * mean.tlm + 
                                  sqrt.sigma11 * q.t_0.50)
    
    cuantiles_peso[i]<- exp(coef.t[1,2] + coef.t[2,2]* mean.edad
                            +coef.t[3,2]+ coef.t[4,2] * mean.tlm + 
                              sqrt.sigma22 * q.t_0.50)
    
    cuantiles_talla[i]<- exp(coef.t[1,3] + coef.t[2,3]* mean.edad
                             + coef.t[3,3]  + coef.t[4,3] * mean.tlm + 
                               sqrt.sigma33 * q.t_0.50)
    
  }
  return(list(cuantiles_per_braq, cuantiles_peso, cuantiles_talla))
}


cuantil_50_niños <- calculo_cuantiles_niños_50(lista_estimaciones_150_t, 
                                               matriz_diseño_150_t)

median(cuantil_50_niños[[1]])
median(cuantil_50_niños[[2]])
median(cuantil_50_niños[[3]])



m4 <- median(cuantil_50_niños[[1]])
mades4 <- abs(cuantil_50_niños[[1]] - m4)
round(median(mades4), 3)
m4

m5 <- median(cuantil_50_niños[[2]])
mades5 <- abs(cuantil_50_niños[[2]] - m5)
round(median(mades5), 3)
m5


m6 <- median(cuantil_50_niños[[3]])
mades6 <- abs(cuantil_50_niños[[3]] - m6)
round(median(mades6), 3)
m6

#-------------------------------------------------------------------------------


#--------------------calculo_cuantiles_slash_75_niños---------------------------------

calculo_cuantiles_niños_75 <- function(lista_estimaciones, X_dis){
  mean.edad <- mean(X_dis[, 2])
  mean.tlm <- mean(X_dis[, 4])
  
  cuantiles_per_braq <- c()
  cuantiles_peso <- c()
  cuantiles_talla <- c()
  
  for (i in 1:length(lista_estimaciones)){
    coef.t <- lista_estimaciones[[i]][[1]]
    matdisp.t <- lista_estimaciones[[i]][[2]]
    nu.t <- lista_estimaciones[[i]][[3]]
    q.t_0.25 <- qt(p = 0.25, df = nu.t)
    q.t_0.75 <- -q.t_0.25
    
    sqrt.sigma11 <- sqrt(matdisp.t[1,1])
    sqrt.sigma22 <- sqrt(matdisp.t[2,2])
    sqrt.sigma33 <- sqrt(matdisp.t[3,3])
    
    cuantiles_per_braq[i]<- exp(coef.t[1,1] + coef.t[2,1]* mean.edad
                                +coef.t[3,1]+ coef.t[4,1] * mean.tlm + 
                                  sqrt.sigma11 * q.t_0.75)
    
    cuantiles_peso[i]<- exp(coef.t[1,2] + coef.t[2,2]* mean.edad
                            +coef.t[3,2]+ coef.t[4,2] * mean.tlm + 
                              sqrt.sigma22 * q.t_0.75)
    
    cuantiles_talla[i]<- exp(coef.t[1,3] + coef.t[2,3]* mean.edad
                             +coef.t[3,3]+ coef.t[4,3] * mean.tlm + 
                               sqrt.sigma33 * q.t_0.75)
    
  }
  return(list(cuantiles_per_braq, cuantiles_peso, cuantiles_talla))
}


cuantil_75_niños <- calculo_cuantiles_niños_75(lista_estimaciones_150_t, 
                                               matriz_diseño_150_t)

median(cuantil_75_niños[[1]])
median(cuantil_75_niños[[2]])
median(cuantil_75_niños[[3]])

m7 <- median(cuantil_75_niños[[1]])
mades7 <- abs(cuantil_75_niños[[1]] - m7)
round(median(mades7), 3)
m7

m8 <- median(cuantil_75_niños[[2]])
mades8 <- abs(cuantil_75_niños[[2]] - m8)
round(median(mades8), 3)
m8


m9 <- median(cuantil_75_niños[[3]])
mades9 <- abs(cuantil_75_niños[[3]] - m9)
round(median(mades9), 3)
m9
