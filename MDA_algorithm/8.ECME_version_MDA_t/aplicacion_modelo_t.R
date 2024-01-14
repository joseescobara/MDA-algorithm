#ECME version del MDA
ECME_version_MDA_t <- function(y_mp, x_mp, B_0, M_0){
  #1. I Step---------------------------------------------------------
  nu <- muestreador_nu(y_mp, x_mp, B_0, M_0) #simulo el nu
  pesos <- imputa_pesos(y_mp, x_mp, B_0, M_0, nu)#simulamos 
  #los pesos de su distribucion posterior
  nuevos_datos <- datos_imputados(y_mp, x_mp, pesos, B_0, M_0) #imputamos
  #los datos que destruyen el patrón monótono
  #2. P Step-----------------------------------------------------------
  dis_post_mvcov<- matriz_sigma(nuevos_datos, x_mp , pesos) #simulamos la 
  #matriz de varianza y covarianza de su distribucion posterior
  M_1 <- dis_post_mvcov[[1]] #muestra de la matriz de mvcov de su dist 
  #posterior
  H <- dis_post_mvcov[[2]]
  B_1 <- generacion_betas(nuevos_datos, x_mp, pesos, H) #simulo la matriz
  #de coeficientes del modelo dado la simulacion de la matriz de varianzas
  #y covarianzas
  simulaciones <- list(B_1, M_1, nu)
  return(simulaciones)
}

#------------------------------------------------------------------------
#esta función me hace Nsim simulaciones montercarlo de las distribuciones
#posteriores de los paramétros del modelo. 
sim_t <- function(y_mp, x_mp, B_0, M_0, Nsim = 10000){
  sim_t_MC <- list() #creo la lista donde voy a guardar las simulacones
  for (i in 1:Nsim) {
    sim_parametros <- ECME_version_MDA_t(y_mp, x_mp, B_0, M_0)
    #aplico la iteracion i-esima para muestrar los paramétros de su 
    #distribucion posterior
    B_0 <- sim_parametros[[1]] #actualizo los paramteros abternidos en la
    #iteracion i
    M_0 <- sim_parametros[[2]] #actualizo los paramteros abternidos en la
    #iteracion i
    nu_0 <- sim_parametros[[3]] #actualizo los paramteros abternidos en la
    #iteracion i
    sim_t_MC[[i]] <- list(B_0, M_0, nu_0)
  }
  
  sim_t_MC <- sim_t_MC[(Nsim*0.10+1):Nsim] 
  return(sim_t_MC)
}


#---------------------------------------------------------------------
#------------------------------------------------------------------------
#funcion para calcular la mediana de cada uno de los paramétros del modelo
#bajo las Nsim muestras montercalos
medianas_modelo_t <- function(y_mp, x_mp, B_0, M_0, Nsim = 10000){
  simulaciones <- sim_t(y_mp, x_mp, B_0, M_0,  Nsim) #aquí aplico
  #la función anterior, para luego sacarle la mediana a las muestras
  #obtenidas
  matriz_coeficientes <- list()
  matriz_varianzas <- list()
  nu <- c()
  for (i in 1:length(simulaciones)) {
    matriz_coeficientes[[i]] <- simulaciones[[i]][[1]] #matriz de betas 
    matriz_varianzas[[i]] <- simulaciones[[i]][[2]] #matriz de varianzas y covarianzas
    nu[i] <- simulaciones[[i]][[3]]
  }
  
  A <- matrix(0, ncol = nrow(as.matrix(matriz_coeficientes[[1]]))*
                ncol(as.matrix(matriz_coeficientes[[1]])), 
              nrow = length(matriz_coeficientes))
  
  B <- matrix(0, ncol = (ncol(y_mp)*(ncol(y_mp)+1))/2,
              nrow = length(matriz_varianzas))
  colnames(A) <- c("B_01", "B_11", "B_21", "B_31",
                   "B_02", "B_12", "B_22", "B_32",
                   "B_03", "B_13", "B_23", "B_33")
  colnames(B) <- c("sigma_11", "sigma_21", "sigma_31", "sigma_22",
                   "sigma_32", "sigma_33")
  for (i in 1:length(matriz_coeficientes)) {
    coef <-  vec(matriz_coeficientes[[i]])
    A[i, ] <- coef
  }
  
  for (i in 1:length(matriz_varianzas)) {
    mcov <-  vech(matriz_varianzas[[i]])
    B[i, ] <- mcov
  }
  
  medianas_coeficientes <- apply(A, 2, median)
  median_var <- apply(B, 2, median)
  coef_estimados <- t(matrix(medianas_coeficientes, byrow = T, nrow = 3))
  var_estimadas <- matrix(c(median_var[1], median_var[2], median_var[3],
                            median_var[2], median_var[4], median_var[5],
                            median_var[3], median_var[5], median_var[6]), 
                          byrow = T, nrow = 3)
  colnames(coef_estimados) <- colnames(matriz_coeficientes[[1]])
  rownames(coef_estimados) <- rownames(matriz_coeficientes[[1]])
  
  colnames(var_estimadas) <- colnames(y_mp)
  rownames(var_estimadas) <- colnames(y_mp)
  
  nu_est <- median(nu)
  estimaciones <- list(coef_estimados, var_estimadas, nu_est)
  return(estimaciones)
}


#------------------------------------------------------------------------------

medianas <- function(y_mp, simulaciones){
  #obtenidas
  Nsim = length(simulaciones)
  matriz_coeficientes <- list()
  matriz_varianzas <- list()
  nu <- c()
  for (i in 1:Nsim) {
    matriz_coeficientes[[i]] <- simulaciones[[i]][[1]] #matriz de betas 
    matriz_varianzas[[i]] <- simulaciones[[i]][[2]] #matriz de varianzas y covarianzas
    nu[i] <- simulaciones[[i]][[3]]
  }
  
  A <- matrix(0, ncol = nrow(as.matrix(matriz_coeficientes[[1]]))*
                ncol(as.matrix(matriz_coeficientes[[1]])), 
              nrow = length(matriz_coeficientes))
  
  B <- matrix(0, ncol = (ncol(y_mp)*(ncol(y_mp)+1))/2,
              nrow = length(matriz_varianzas))
  colnames(A) <- c("B_01", "B_11", "B_21", "B_31",
                   "B_02", "B_12", "B_22", "B_32",
                   "B_03", "B_13", "B_23", "B_33")
  colnames(B) <- c("sigma_11", "sigma_21", "sigma_31", "sigma_22",
                   "sigma_32", "sigma_33")
  for (i in 1:length(matriz_coeficientes)) {
    coef <-  vec(matriz_coeficientes[[i]])
    A[i, ] <- coef
  }
  
  for (i in 1:length(matriz_varianzas)) {
    mcov <-  vech(matriz_varianzas[[i]])
    B[i, ] <- mcov
  }
  
  medianas_coeficientes <- apply(A, 2, median)
  median_var <- apply(B, 2, median)
  coef_estimados <- t(matrix(medianas_coeficientes, byrow = T, nrow = 3))
  var_estimadas <- matrix(c(median_var[1], median_var[2], median_var[3],
                            median_var[2], median_var[4], median_var[5],
                            median_var[3], median_var[5], median_var[6]), 
                          byrow = T, nrow = 3)
  colnames(coef_estimados) <- colnames(matriz_coeficientes[[1]])
  rownames(coef_estimados) <- rownames(matriz_coeficientes[[1]])
  
  colnames(var_estimadas) <- colnames(y_mp)
  rownames(var_estimadas) <- colnames(y_mp)
  
  nu_est <- median(nu)
  estimaciones <- list(coef_estimados, var_estimadas, nu_est)
  return(estimaciones)
}















