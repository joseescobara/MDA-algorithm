setwd("ruta")
load("parametros_reales_slash.Rdata")
load("aplicacion_modelo_slash.Rdata")

coef.slash <- estim[[1]]

#---------------------------muestras_coef--------------------------------------
#------------------------------------------------------------------------------
simulaciones <- aplicacion_slash
Nsim <- length(simulaciones)
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

colnames(A) <- c("B_01", "B_11", "B_21", "B_31",
                 "B_02", "B_12", "B_22", "B_32",
                 "B_03", "B_13", "B_23", "B_33")
for (i in 1:length(matriz_coeficientes)) {
  coef <-  vec(matriz_coeficientes[[i]])
  A[i, ] <- coef
}

#------------------------------------------------------------------------------
#arm circunference coefficients
intercept.armcircunferences <- coef.slash[1,1] 
age.armcircunferece <- coef.slash[2,1]
gender.armcircunference <- coef.slash[3,1]
tleche.armcircunference <- coef.slash[4,1]

#wigth coefficients

intercept.weight <- coef.slash[1,2] 
age.weight <- coef.slash[2,2]
gender.weight <- coef.slash[3,2]
tleche.weight <- coef.slash[4,2]

#length coefficients

intercept.length <- coef.slash[1,3]
age.length <- coef.slash[2,3]
gender.length <- coef.slash[3,3]
tleche.length <- coef.slash[4,3]

#--------------------------------standard error (SE)---------------------------

#arm circunference SE
SE.intercept.armcircunferences <- sum((A[, 1]-median(A[, 1]))^{2})/length(A[, 1]) 
SE.age.armcircunferece <- sum((A[, 2] - median(A[, 2]))^{2})/length(A[, 2])
SE.gender.armcircunference <- sum((A[, 3] - median(A[, 3]))^{2})/length(A[, 3])
SE.tleche.armcircunference <- sum((A[, 4] - median(A[, 4]))^{2})/length(A[, 4])

#wigth SE

SE.intercept.weight <- sum((A[, 5] - median(A[, 5]))^{2})/length(A[, 5]) 
SE.age.weight <- sum((A[, 6] - median(A[, 6]))^{2})/length(A[, 6])
SE.gender.weight <- sum((A[, 7] - median(A[, 7]))^{2})/length(A[, 7])
SE.tleche.weight <- sum((A[, 8] - median(A[, 8]))^{2})/length(A[, 8])

#length SE

SE.intercept.length <- sum((A[, 9] - median(A[, 9]))^{2})/length(A[, 9])
SE.age.length <- sum((A[, 10] - median(A[, 10]))^{2})/length(A[, 10])
SE.gender.length <- sum((A[, 11] - median(A[, 11]))^{2})/length(A[, 11])
SE.tleche.length <- sum((A[, 12] - median(A[, 12]))^{2})/length(A[, 12])

#--------------------------------probability interval--------------------------

#arm circunference probability interval

low.lim.intercept.armcircunference <- quantile(A[, 1], 0.025)
up.lim.intercept.armcircunference <- quantile(A[, 1], 0.975)


low.lim.age.armcircunference <- quantile(A[, 2], 0.025)
up.lim.age.armcircunference <- quantile(A[, 2], 0.975)

low.lim.gender.armcircunference <- quantile(A[, 3], 0.025)
up.lim.gender.armcircunference <- quantile(A[, 3], 0.975)

low.lim.tleche.armcircunference <- quantile(A[, 4], 0.025)
up.lim.tleche.armcircunference <- quantile(A[, 4], 0.975)

#weight probability interval

low.lim.intercept.weight <- quantile(A[, 5], 0.025)
up.lim.intercept.weight <- quantile(A[, 5], 0.975)

low.lim.age.weight <- quantile(A[, 6], 0.025)
up.lim.age.weight <- quantile(A[, 6], 0.975)

low.lim.gender.weight <- quantile(A[, 7], 0.025)
up.lim.gender.weight <- quantile(A[, 7], 0.975)

low.lim.tleche.weight <- quantile(A[, 8], 0.025)
up.lim.tleche.weight <- quantile(A[, 8], 0.975)

#length probability interval

low.lim.intercept.length <- quantile(A[, 9], 0.025)
up.lim.intercept.length <- quantile(A[, 9], 0.975)

low.lim.age.length <- quantile(A[, 10], 0.025)
up.lim.age.length <- quantile(A[, 10], 0.975)

low.lim.gender.length <- quantile(A[, 11], 0.025)
up.lim.gender.length <- quantile(A[, 11], 0.975)

low.lim.tleche.length <- quantile(A[, 12], 0.025)
up.lim.tleche.length <- quantile(A[, 12], 0.975)

#=======================================================
# Table 4 (Estimate, asymptotic standard error (SE),  
# lower and upper bounds of asymptotic 95% confidence  
# interval and p-value of the Wald test associated to 
# bivariate log-slash linear regression model)        
#=======================================================
estim.armcircunference <- c(intercept.armcircunferences, 
                  age.armcircunferece,
                  gender.armcircunference,
                  tleche.armcircunference)

SE.armcircunference <- c(SE.intercept.armcircunferences,
               SE.age.armcircunferece,
               SE.gender.armcircunference,
               SE.tleche.armcircunference)

low.armcircunference <- c(low.lim.intercept.armcircunference,
                low.lim.age.armcircunference,
                low.lim.gender.armcircunference,
                low.lim.tleche.armcircunference)

up.armcircunference <- c(up.lim.intercept.armcircunference,
               up.lim.age.armcircunference,
               up.lim.gender.armcircunference,
               up.lim.tleche.armcircunference)

expl.vbles <- c("intercept", "age", "gender", "tleche")
table_estim.armcircunference <- data.frame(estim.armcircunference,
                                 SE.armcircunference,
                                 low.armcircunference,
                                 up.armcircunference,
                                 row.names = expl.vbles)
table_estim.armcircunference <- signif(round(table_estim.armcircunference, 5), 6)
colnames(table_estim.armcircunference) <- c("Estimate", "SE", "Lower", "Upper")



#-----------------------------------------------------------------------------
estim.weight <- c(intercept.weight, 
                  age.weight,
                  gender.weight,
                  tleche.weight)

SE.weight <- c(SE.intercept.weight,
               SE.age.weight,
               SE.gender.weight,
               SE.tleche.weight)

low.weight <- c(low.lim.intercept.weight,
                low.lim.age.weight,
                low.lim.gender.weight,
                low.lim.tleche.weight)

up.weight <- c(up.lim.intercept.weight,
               up.lim.age.weight,
               up.lim.gender.weight,
               up.lim.tleche.weight)

expl.vbles <- c("intercept", "age", "gender", "tleche")
table_estim.weight <- data.frame(estim.weight,
                                 SE.weight,
                                 low.weight,
                                 up.weight,
                                 row.names = expl.vbles)
table_estim.weight <- signif(round(table_estim.weight, 5), 6)
colnames(table_estim.weight) <- c("Estimate", "SE", "Lower", "Upper")
#------------------------------------------------------------------------------
estim.length <- c(intercept.length,
                  age.length,
                  gender.length,
                  tleche.length)

SE.length <- c(SE.intercept.length,
               SE.age.length,
               SE.gender.length,
               SE.tleche.length)

low.length <- c(low.lim.intercept.length,
                low.lim.age.length,
                low.lim.gender.length,
                low.lim.tleche.length)

up.length <- c(up.lim.intercept.length,
               up.lim.age.length,
               up.lim.gender.length,
               up.lim.tleche.length)

expl.vbles <- c("intercept", "age","gender", "tleche")
table_estim.length <- data.frame(estim.length, 
                                 SE.length, 
                                 low.length,
                                 up.length, 
                                 row.names = expl.vbles)
table_estim.length <- signif(round(table_estim.length, 5), 6)
colnames(table_estim.length) <- c("Estimate", "SE", "Lower", "Upper")

table_estim <- list(table_estim.armcircunference,
                    table_estim.weight, table_estim.length)
names(table_estim) <- c("armcircuf.","weight", "length")
table4 <- table_estim
capture.output(print(table4, print.gap=3), file = "C:/Users/Jose/OneDrive - Universidad de Antioquia/Escritorio/Maestria/Tesis/Algoritmo MDA/Aplicacion_Quantile/Datos Reales/modelo_slash_application/Table_4.txt")
