#distribución posterior de los parámetros del modelo slash
setwd("ruta")
load(r"ruta")#base de datos real
load("aplicacion_modelo_t.Rdata")


load("estimaciones_50_modelo_t.Rdata")
load("estimaciones_100_modelo_t.Rdata")
load("estimaciones_150_modelo_t.Rdata")


library(fBasics)


datos <- base_de_datos #datos completos
datos <- as.matrix(datos)

#respuestas y covariables
respuestas <- datos[, 1:3] #matriz de variables respuesta
covariables <- cbind(rep(1, times = nrow(respuestas)), datos[,9], 
                     datos[, 6], datos[,8]) #matriz de covariables
colnames(covariables) <- c("X_0", "edad", "sexo", "tiempo_lechem")


#----------------------preparacion posterior real----------------------------------
#extraemos la cadena de cada uno de los coeficientes
simulaciones <- aplicacion_t
y_mp <- respuestas
#obtenidas
Nsim <- length(simulaciones)
matriz_coeficientes <- list()
for (i in 1:Nsim) {
  matriz_coeficientes[[i]] <- simulaciones[[i]][[1]] #matriz de betas 
  
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
#----------------------preparacion posterior 50----------------------------------
#extraemos la cadena de cada uno de los coeficientes
simulaciones1 <- lista_estimaciones_50_t
#obtenidas
Nsim1 <- length(simulaciones1)
matriz_coeficientes1 <- list()
for (i in 1:Nsim1) {
  matriz_coeficientes1[[i]] <- simulaciones1[[i]][[1]] #matriz de betas 
  
}
A1 <- matrix(0, ncol = nrow(as.matrix(matriz_coeficientes1[[1]]))*
               ncol(as.matrix(matriz_coeficientes1[[1]])), 
             nrow = length(matriz_coeficientes1))
colnames(A1) <- c("B_01", "B_11", "B_21", "B_31",
                  "B_02", "B_12", "B_22", "B_32",
                  "B_03", "B_13", "B_23", "B_33")
for (i in 1:length(matriz_coeficientes1)) {
  coef1 <-  vec(matriz_coeficientes1[[i]])
  A1[i, ] <- coef1
}

#------------------------------------------------------------------------------
#----------------------preparacion posterior 100----------------------------------
#extraemos la cadena de cada uno de los coeficientes
simulaciones2 <- lista_estimaciones_100_t
#obtenidas
Nsim2 <- length(simulaciones2)
matriz_coeficientes2 <- list()
for (i in 1:Nsim2) {
  matriz_coeficientes2[[i]] <- simulaciones2[[i]][[1]] #matriz de betas 
  
}
A2 <- matrix(0, ncol = nrow(as.matrix(matriz_coeficientes2[[1]]))*
               ncol(as.matrix(matriz_coeficientes2[[1]])), 
             nrow = length(matriz_coeficientes2))
colnames(A2) <- c("B_01", "B_11", "B_21", "B_31",
                  "B_02", "B_12", "B_22", "B_32",
                  "B_03", "B_13", "B_23", "B_33")
for (i in 1:length(matriz_coeficientes2)) {
  coef2 <-  vec(matriz_coeficientes2[[i]])
  A2[i, ] <- coef2
}

#------------------------------------------------------------------------------
#----------------------preparacion posterior 150----------------------------------
#extraemos la cadena de cada uno de los coeficientes
simulaciones3 <- lista_estimaciones_150_t
#obtenidas
Nsim3 <- length(simulaciones3)
matriz_coeficientes3 <- list()
for (i in 1:Nsim3) {
  matriz_coeficientes3[[i]] <- simulaciones3[[i]][[1]] #matriz de betas 
  
}
A3 <- matrix(0, ncol = nrow(as.matrix(matriz_coeficientes3[[1]]))*
               ncol(as.matrix(matriz_coeficientes3[[1]])), 
             nrow = length(matriz_coeficientes3))
colnames(A3) <- c("B_01", "B_11", "B_21", "B_31",
                  "B_02", "B_12", "B_22", "B_32",
                  "B_03", "B_13", "B_23", "B_33")
for (i in 1:length(matriz_coeficientes3)) {
  coef3 <-  vec(matriz_coeficientes3[[i]])
  A3[i, ] <- coef3
}


#-------------------------------------------B_01------------------------------

library(ggplot2)
B_01 <- A[, 1]
B_01_50 <- A1[, 1]
B_01_100 <- A2[, 1]
B_01_150 <- A3[, 1]

#-----------------------------------------------------------------------------
df <- data.frame(
  valor = c(A[, 1], A1[, 1]),
  densidad = rep(c("verdadero", "muestra 50"), times = c(length(B_01), length(B_01_50)))
)

# Graficar las curvas de densidad
ggplot(df, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[0][1])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-----------------------------------------------------------------------------
df1 <- data.frame(
  valor = c(A[, 1], A2[, 1]),
  densidad = rep(c("verdadero", "muestra 100"), times = c(length(B_01), length(B_01_100)))
)

# Graficar las curvas de densidad
ggplot(df1, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[0][1])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#----------------------------------------------------------------------------
df2 <- data.frame(
  valor = c(A[, 1], A3[, 1]),
  densidad = rep(c("verdadero", "muestra 150"), times = c(length(B_01), length(B_01_150)))
)

# Graficar las curvas de densidad
ggplot(df2, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[0][1])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-------------------------------------------B_11------------------------------

B_11 <- A[, 2]
B_11_50 <- A1[, 2]
B_11_100 <- A2[, 2]
B_11_150 <- A3[, 2]



#-----------------------------------------------------------------------------
df4 <- data.frame(
  valor = c(B_11, B_11_50),
  densidad = rep(c("verdadero", "muestra 50"), times = c(length(B_11), length(B_11_50)))
)

# Graficar las curvas de densidad
ggplot(df4, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[1][1])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-----------------------------------------------------------------------------
df5 <- data.frame(
  valor = c(B_11, B_11_100),
  densidad = rep(c("verdadero", "muestra 100"), times = c(length(B_11), length(B_11_100)))
)

# Graficar las curvas de densidad
ggplot(df5, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[1][1])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#----------------------------------------------------------------------------
df6 <- data.frame(
  valor = c(B_11, B_11_150),
  densidad = rep(c("verdadero", "muestra 150"), times = c(length(B_11), length(B_11_150)))
)

# Graficar las curvas de densidad
ggplot(df6, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[1][1])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))


#-------------------------------------------B_21------------------------------

B_21 <- A[, 3]
B_21_50 <- A1[, 3]
B_21_100 <- A2[, 3]
B_21_150 <- A3[, 3]



#-----------------------------------------------------------------------------
df7 <- data.frame(
  valor = c(B_21, B_21_50),
  densidad = rep(c("verdadero", "muestra 50"), times = c(length(B_21), length(B_21_50)))
)

# Graficar las curvas de densidad
ggplot(df7, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[2][1])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-----------------------------------------------------------------------------
df8 <- data.frame(
  valor = c(B_21, B_21_100),
  densidad = rep(c("verdadero", "muestra 100"), times = c(length(B_21), length(B_21_100)))
)

# Graficar las curvas de densidad
ggplot(df8, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[2][1])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#------------------------------------------------------------------------------
df9 <- data.frame(
  valor = c(B_21, B_21_150),
  densidad = rep(c("verdadero", "muestra 150"), times = c(length(B_21), length(B_21_150)))
)

# Graficar las curvas de densidad
ggplot(df9, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[2][1])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-------------------------------------------B_31------------------------------

B_31 <- A[, 4]
B_31_50 <- A1[, 4]
B_31_100 <- A2[, 4]
B_31_150 <- A3[, 4]



#-----------------------------------------------------------------------------
df10 <- data.frame(
  valor = c(B_31, B_31_50),
  densidad = rep(c("verdadero", "muestra 50"), times = c(length(B_31), length(B_31_50)))
)

# Graficar las curvas de densidad
ggplot(df10, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[3][1])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-----------------------------------------------------------------------------
df11 <- data.frame(
  valor = c(B_31, B_31_100),
  densidad = rep(c("verdadero", "muestra 100"), times = c(length(B_31), length(B_31_100)))
)

# Graficar las curvas de densidad
ggplot(df11, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[3][1])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#----------------------------------------------------------------------------
df12 <- data.frame(
  valor = c(B_31, B_31_150),
  densidad = rep(c("verdadero", "muestra 150"), times = c(length(B_31), length(B_31_150)))
)

# Graficar las curvas de densidad
ggplot(df12, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[3][1])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-------------------------------------------B_02------------------------------

B_02 <- A[, 5]
B_02_50 <- A1[, 5]
B_02_100 <- A2[, 5]
B_02_150 <- A3[, 5]



#-----------------------------------------------------------------------------
df13 <- data.frame(
  valor = c(B_02, B_02_50),
  densidad = rep(c("verdadero", "muestra 50"), times = c(length(B_02), length(B_02_50)))
)

# Graficar las curvas de densidad
ggplot(df13, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[0][2])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-----------------------------------------------------------------------------
df14 <- data.frame(
  valor = c(B_02, B_02_100),
  densidad = rep(c("verdadero", "muestra 100"), times = c(length(B_02), length(B_02_100)))
)

# Graficar las curvas de densidad
ggplot(df14, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[0][2])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#----------------------------------------------------------------------------
df15 <- data.frame(
  valor = c(B_02, B_02_150),
  densidad = rep(c("verdadero", "muestra 150"), times = c(length(B_02), length(B_02_150)))
)

# Graficar las curvas de densidad
ggplot(df15, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[0][2])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))


#-------------------------------------------B_12------------------------------

B_12 <- A[, 6]
B_12_50 <- A1[, 6]
B_12_100 <- A2[, 6]
B_12_150 <- A3[, 6]



#-----------------------------------------------------------------------------
df16 <- data.frame(
  valor = c(B_12, B_12_50),
  densidad = rep(c("verdadero", "muestra 50"), times = c(length(B_12), length(B_12_50)))
)

# Graficar las curvas de densidad
ggplot(df16, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[1][2])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-----------------------------------------------------------------------------
df17 <- data.frame(
  valor = c(B_12, B_12_100),
  densidad = rep(c("verdadero", "muestra 100"), times = c(length(B_12), length(B_12_100)))
)

# Graficar las curvas de densidad
ggplot(df17, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[1][2])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#----------------------------------------------------------------------------
df18 <- data.frame(
  valor = c(B_12, B_12_150),
  densidad = rep(c("verdadero", "muestra 150"), times = c(length(B_12), length(B_12_150)))
)

# Graficar las curvas de densidad
ggplot(df18, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[1][2])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-------------------------------------------B_22------------------------------

B_22 <- A[, 7]
B_22_50 <- A1[, 7]
B_22_100 <- A2[, 7]
B_22_150 <- A3[, 7]



#-----------------------------------------------------------------------------
df19 <- data.frame(
  valor = c(B_22, B_22_50),
  densidad = rep(c("verdadero", "muestra 50"), times = c(length(B_22), length(B_22_50)))
)

# Graficar las curvas de densidad
ggplot(df19, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[2][2])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-----------------------------------------------------------------------------
df20 <- data.frame(
  valor = c(B_22, B_22_100),
  densidad = rep(c("verdadero", "muestra 100"), times = c(length(B_22), length(B_22_100)))
)

# Graficar las curvas de densidad
ggplot(df20, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[2][2])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#----------------------------------------------------------------------------
df21 <- data.frame(
  valor = c(B_22, B_22_150),
  densidad = rep(c("verdadero", "muestra 150"), times = c(length(B_22), length(B_22_150)))
)

# Graficar las curvas de densidad
ggplot(df21, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[2][2])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))


#-------------------------------------------B_32------------------------------

B_32 <- A[, 8]
B_32_50 <- A1[, 8]
B_32_100 <- A2[, 8]
B_32_150 <- A3[, 8]



#-----------------------------------------------------------------------------
df22 <- data.frame(
  valor = c(B_32, B_32_50),
  densidad = rep(c("verdadero", "muestra 50"), times = c(length(B_32), length(B_32_50)))
)

# Graficar las curvas de densidad
ggplot(df22, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[3][2])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-----------------------------------------------------------------------------
df23 <- data.frame(
  valor = c(B_32, B_32_100),
  densidad = rep(c("verdadero", "muestra 100"), times = c(length(B_32), length(B_32_100)))
)

# Graficar las curvas de densidad
ggplot(df23, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[3][2])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#----------------------------------------------------------------------------
df24 <- data.frame(
  valor = c(B_32, B_32_150),
  densidad = rep(c("verdadero", "muestra 150"), times = c(length(B_32), length(B_32_150)))
)

# Graficar las curvas de densidad
ggplot(df24, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[3][2])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-------------------------------------------B_03------------------------------

B_03 <- A[, 9]
B_03_50 <- A1[, 9]
B_03_100 <- A2[, 9]
B_03_150 <- A3[, 9]



#-----------------------------------------------------------------------------
df25 <- data.frame(
  valor = c(B_03, B_03_50),
  densidad = rep(c("verdadero", "muestra 50"), times = c(length(B_03), length(B_03_50)))
)

# Graficar las curvas de densidad
ggplot(df25, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[0][3])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-----------------------------------------------------------------------------
df26 <- data.frame(
  valor = c(B_03, B_03_100),
  densidad = rep(c("verdadero", "muestra 100"), times = c(length(B_03), length(B_03_100)))
)

# Graficar las curvas de densidad
ggplot(df26, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[0][3])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#----------------------------------------------------------------------------
df27 <- data.frame(
  valor = c(B_03, B_03_150),
  densidad = rep(c("verdadero", "muestra 150"), times = c(length(B_03), length(B_03_150)))
)

# Graficar las curvas de densidad
ggplot(df27, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[0][3])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-------------------------------------------B_13------------------------------

B_13 <- A[, 10]
B_13_50 <- A1[, 10]
B_13_100 <- A2[, 10]
B_13_150 <- A3[, 10]



#-----------------------------------------------------------------------------
df28 <- data.frame(
  valor = c(B_13, B_13_50),
  densidad = rep(c("verdadero", "muestra 50"), times = c(length(B_13), length(B_13_50)))
)

# Graficar las curvas de densidad
ggplot(df28, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[1][3])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-----------------------------------------------------------------------------
df29 <- data.frame(
  valor = c(B_13, B_13_100),
  densidad = rep(c("verdadero", "muestra 100"), times = c(length(B_13), length(B_13_100)))
)

# Graficar las curvas de densidad
ggplot(df29, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[1][3])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#----------------------------------------------------------------------------
df30 <- data.frame(
  valor = c(B_13, B_13_150),
  densidad = rep(c("verdadero", "muestra 150"), times = c(length(B_13), length(B_13_150)))
)

# Graficar las curvas de densidad
ggplot(df30, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[1][3])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-------------------------------------------B_23------------------------------

B_23 <- A[, 11]
B_23_50 <- A1[, 11]
B_23_100 <- A2[, 11]
B_23_150 <- A3[, 11]



#-----------------------------------------------------------------------------
df31 <- data.frame(
  valor = c(B_23, B_23_50),
  densidad = rep(c("verdadero", "muestra 50"), times = c(length(B_23), length(B_23_50)))
)

# Graficar las curvas de densidad
ggplot(df31, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[2][3])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-----------------------------------------------------------------------------
df32 <- data.frame(
  valor = c(B_23, B_23_100),
  densidad = rep(c("verdadero", "muestra 100"), times = c(length(B_23), length(B_23_100)))
)

# Graficar las curvas de densidad
ggplot(df32, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[2][3])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#----------------------------------------------------------------------------
df33 <- data.frame(
  valor = c(B_23, B_23_150),
  densidad = rep(c("verdadero", "muestra 150"), times = c(length(B_23), length(B_23_150)))
)

# Graficar las curvas de densidad
ggplot(df33, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[2][3])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-------------------------------------------B_33------------------------------

B_33 <- A[, 12]
B_33_50 <- A1[, 12]
B_33_100 <- A2[, 12]
B_33_150 <- A3[, 12]



#-----------------------------------------------------------------------------
df34 <- data.frame(
  valor = c(B_33, B_33_50),
  densidad = rep(c("verdadero", "muestra 50"), times = c(length(B_33), length(B_33_50)))
)

# Graficar las curvas de densidad
ggplot(df34, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[3][3])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#-----------------------------------------------------------------------------
df35 <- data.frame(
  valor = c(B_33, B_33_100),
  densidad = rep(c("verdadero", "muestra 100"), times = c(length(B_33), length(B_33_100)))
)

# Graficar las curvas de densidad
ggplot(df35, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[3][3])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))

#----------------------------------------------------------------------------
df36 <- data.frame(
  valor = c(B_33, B_33_150),
  densidad = rep(c("verdadero", "muestra 150"), times = c(length(B_33), length(B_33_150)))
)

# Graficar las curvas de densidad
ggplot(df36, aes(x = valor, fill = densidad, linetype = densidad)) +
  geom_density(alpha = 0.0) +
  theme_minimal() +
  scale_linetype_manual(values = c("dashed", "solid" )) +
  labs(
    x = expression(paste(beta[3][3])),
    y = "Densidad") +
  theme(legend.position = "none")+
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1))
