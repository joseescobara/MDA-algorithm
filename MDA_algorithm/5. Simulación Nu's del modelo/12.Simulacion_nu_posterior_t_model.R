
#-----------------------------------------------------------------------
#------------------------------------------------------------------------------
#funcion posterior del hipeparámetro
posterior_nu <- function(y_df, x_df, B_0, M_0, nu){
  n <- nrow(y_df)
  dista <- dista_maha_indiv(y_df, x_df, B_0, M_0)
  k_1 <- 0 #para calcular las sumatorias
  k_2 <- 0 #para calcular las sumatorias
  for (i in 1:n){
    p_i_k <-length(which(is.na( y_df[i,]) == F ))
    k_1 <- k_1 + log(gamma((nu + p_i_k)/2)) #primera sumatoria de la densidad
    k_2 <- k_2 + (nu + p_i_k)/2 * log(nu + dista[i]) #segunda sumatoria
    #de la densidad
  }
  prior <- 1/nu^2 #nu >= 1
  densidad <- exp(-log(prior) + k_1 -n*log(gamma(nu/2)) + (n*nu)/2 * log(nu)
                  -k_2)
  return(densidad)
}

#-----------------------------------------------------------------------------
#función para muestrear nu de su dist posterior
muestreador_nu <- function(y_df, x_df, B_0, M_0){
  nu <- seq(2.5, 15, 0.0125)
  densidades <- posterior_nu(y_df, x_df, B_0, M_0,nu)
  probabilidades <- densidades/sum(densidades)
  muestra <- sample(nu, 1, prob = probabilidades)
  return(muestra)
}


