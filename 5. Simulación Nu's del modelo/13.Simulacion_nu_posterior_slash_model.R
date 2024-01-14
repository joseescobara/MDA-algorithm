#simulacion del nu de su distribucion posterior
nu_modelo_slash <- function(a = 6, b = 2, pesos){
  n <- length(pesos) #numero de individuos
  ln_pesos <- log(pesos) #logaritmo de los pesos
  muestra <- rgamma(1, shape = a + n, rate = b - sum(ln_pesos))
  return(muestra)
}

