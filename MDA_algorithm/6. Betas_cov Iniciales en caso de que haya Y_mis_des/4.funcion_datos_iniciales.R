
#-----------------------------------------------------------------------------
#está función será utilizada solamente para proponer los datos iniciales
#del algortimo.
#-----------------------------------------------------------------------------
datos_iniciales <- function(y_df, x_df){
  y_df1 <- as.matrix(y_df)
  y_mp <- patron_monotono(y_df1, x_df)[[1]]
  #utilizo las siguientes funciones para imputar los datos faltantes
  s <- prelim.norm(y_df1)   #do preliminary manipulations
  thetahat <- em.norm(s, showits = F)   #find the mle
  rngseed(1234567)   #set random number generator seed
  ximp <- imp.norm(s,thetahat,y_df1)  #impute missing data under the MLE
  k <- patron_monotono(y_df1, x_df)[[3]] #extraigo las
  #filas del ordenamiento final en patron monotono
  ximp2 <- ximp[k,colnames(patron_monotono(y_df1, x_df)[[1]])] #ordeno la
  #matriz imputada de acuerdo a la ordenación final por columnas y por filas
  patrones_m1 <- matriz_patrones(y_mp) #utilizo la función matriz_patrones
  #para extraer las matrices de acuerdo a su patrón, para luego ser utilizada
  #para obtener la matriz final con los datos que destruyen el patrón monótono
  #imputados
  indiv_m1 <- matriz_patrones(y_mp)[[2]] #extraigo por medio de la función
  #matriz patrones, los individuos por patrones reales.
  y_df1_ini <- patron_monotono(y_df1, x_df)[[1]] #extrae la matriz de variables
  #respuestas ordenada en patrón monótono, incluyendo los datos observados y 
  #los datos que destruyen el patrón monótono. 
  y_df1_ini[1:indiv_m1[1],
            colnames(patrones_m1[[1]][[1]])]<-ximp2[1:indiv_m1[1],
                                                    colnames(patrones_m1[[1]][[1]])]
  #en la operación anterior, reemplaza la primera parte de la matriz por parte
  #de la matriz imputada
  acum1 <- indiv_m1[1] 
  for(j in 2:(length(indiv_m1))){ #recorre los patrones desde dos
    y_df1_ini[(acum1 + 1):(acum1 + indiv_m1[(j)]),
              colnames(patrones_m1[[1]][[(j)]])]<-ximp2[(acum1+1) :(acum1 +indiv_m1[(j)]),
                                                        colnames(patrones_m1[[1]][[(j)]])]
    acum1 <- acum1 + indiv_m1[(j)]
  }
  return(y_df1_ini)
}
#-----------------------------------------------------------------------------
