#funcion para calcular los betas iniciales
betas_mvco_init <- function(datos_init, x_mp){
pesos <- rep(1, times = nrow(datos_init))
M <- matriz_sigma(datos_init, x_mp, pesos)
M_0 <- M[[1]]
H <- M[[2]]
B_0 <- generacion_betas(datos_init, x_mp, pesos, H)
return(list(M_0, B_0))
}




